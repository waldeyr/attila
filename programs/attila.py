#!/usr/bin/env python3
"""
attila.py - Native Python 3 orchestrator and pipeline for immunoglobulin sequence analysis.
Replaces autoiganalysis3.pl, translateab9, and other helper scripts (Perl/R/binary).
"""

import os
import sys
import re
import math
import subprocess
import argparse
import shutil
import urllib.request
import urllib.parse
from pathlib import Path

# =====================================================================
# Helper Functions
# =====================================================================

def reverse_complement(seq):
    tab = str.maketrans("AGCTagctNn", "TCGATCGANN")
    return seq.translate(tab)[::-1].upper()

def load_codon_table(gencode_path):
    codon_table = {}
    if not os.path.exists(gencode_path):
        raise FileNotFoundError(f"Codon mapping file not found: {gencode_path}")
    with open(gencode_path, 'r', encoding='utf-8', errors='ignore') as f:
        curr_codon = None
        for line in f:
            line = line.rstrip('\r\n')
            if line.startswith('>'):
                curr_codon = line[1:].strip().upper()
            else:
                if curr_codon:
                    codon_table[curr_codon] = line.strip().upper()
    return codon_table

def translate_codon(codon, codon_table):
    if len(codon) < 3:
        return ''
    return codon_table.get(codon.upper(), 'X')

def translate_seq(dna_seq, frame, sense, codon_table):
    if sense == '-':
        working_dna = reverse_complement(dna_seq)
    else:
        working_dna = dna_seq
        
    offset = frame - 1
    aa_seq = []
    for i in range(offset, len(working_dna) - 2, 3):
        codon = working_dna[i : i + 3]
        aa_seq.append(translate_codon(codon, codon_table))
    return "".join(aa_seq)

def find_cdr3_bounds(aa_seq, libtype):
    # Find all Cysteine indices
    cys_indices = [i for i, aa in enumerate(aa_seq) if aa == 'C']
    
    # FR4 pattern: WG\wG (VH) or FG\wG (VL)
    # Match WG\wG and FG\wG in amino acid sequence
    if libtype == 0:
        motif_pat = re.compile(r'WG[A-Z]G')
    else:
        motif_pat = re.compile(r'FG[A-Z]G')
        
    motif_matches = []
    for match in motif_pat.finditer(aa_seq):
        motif_matches.append(match.start())
        
    # Search C1-C2-FR4 combination
    # C1 to C2: typically 50 to 90 amino acids (meaning 50 <= c2 - c1 - 1 <= 90)
    # C2 to motif: typically 3 to 35 amino acids (meaning 3 <= m_idx - c2 - 1 <= 35)
    for c1 in cys_indices:
        for c2 in cys_indices:
            if c2 <= c1:
                continue
            dist_c1_c2 = c2 - c1 - 1
            if 50 <= dist_c1_c2 <= 90:
                for m_idx in motif_matches:
                    if m_idx <= c2:
                        continue
                    dist_c2_m = m_idx - c2 - 1
                    if 3 <= dist_c2_m <= 35:
                        # c2+1 is 0-indexed start of CDR3; m_idx is 0-indexed start of FR4 (exclusive of CDR3)
                        # We return:
                        # - 0-indexed start (c2+1) and end (m_idx) for slicing
                        # - 1-indexed start (c2+2) and end (m_idx) for the [start-end] header formatting
                        return c2 + 1, m_idx, c2 + 2, m_idx
    return None

def translate_read(seq_id, dna_seq, codon_table, libtype):
    frames = [
        (1, '+', '1+'),
        (2, '+', '2+'),
        (3, '+', '3+'),
        (1, '-', '1-'),
        (2, '-', '2-'),
        (3, '-', '3-')
    ]
    for frame, sense, frame_str in frames:
        aa_seq = translate_seq(dna_seq, frame, sense, codon_table)
        bounds = find_cdr3_bounds(aa_seq, libtype)
        if bounds:
            c2_plus_1, m_idx, cdr_start, cdr_end = bounds
            cdr_seq = aa_seq[c2_plus_1 : m_idx]
            return {
                'frame_str': frame_str,
                'aa_seq': aa_seq,
                'cdr_seq': cdr_seq,
                'cdr_start': cdr_start,
                'cdr_end': cdr_end
            }
    return None

# =====================================================================
# Pipeline Steps
# =====================================================================

def step_join(cfg, pipelinerrorlog):
    print("Executing Step: Join")
    paired = int(cfg.get('pairedend', 0))
    if paired == 0:
        # Single-end: just run FastQC on inputs
        inputs = []
        if 'input1dir' in cfg: inputs.append((cfg['input1dir'], cfg['projectdir'] + '/' + cfg['projectname'] + '/' + ('VH' if int(cfg['libtype'])==0 else 'VL') + '/InitialRound'))
        if 'input2dir' in cfg: inputs.append((cfg['input2dir'], cfg['projectdir'] + '/' + cfg['projectname'] + '/' + ('VH' if int(cfg['libtype'])==0 else 'VL') + '/FinalRound'))
        
        for inp, out_dir in inputs:
            cmd = f"fastqc -q -o {out_dir} {inp}"
            print(f"Running: {cmd}")
            subprocess.run(cmd, shell=True, stderr=open(pipelinerrorlog, 'a'))
    else:
        # Paired-end: run FastQC on raw reads, then fastq-join
        lib_name = 'VH' if int(cfg['libtype']) == 0 else 'VL'
        round_init_dir = f"{cfg['projectdir']}/{cfg['projectname']}/{lib_name}/InitialRound"
        round_final_dir = f"{cfg['projectdir']}/{cfg['projectname']}/{lib_name}/FinalRound"
        
        # FastQC raw reads
        fastqc_files = [cfg['input1r1dir'], cfg['input1r2dir'], cfg['input2r1dir'], cfg['input2r2dir']]
        for f in fastqc_files:
            out_dir = round_init_dir if 'input1' in f else round_final_dir
            cmd = f"fastqc -q -o {out_dir} {f}"
            subprocess.run(cmd, shell=True, stderr=open(pipelinerrorlog, 'a'))
            
        # Fastq-join
        output_prefix1 = f"{round_init_dir}/{Path(cfg['input1r1dir']).stem}"
        output_prefix2 = f"{round_final_dir}/{Path(cfg['input2r1dir']).stem}"
        
        # We need to find fastq-join or fastqjoin command
        join_cmd = "fastq-join"
        # Check if fastq-join or fastqjoin exists
        if subprocess.run("command -v fastq-join", shell=True, stdout=subprocess.DEVNULL).returncode != 0:
            if subprocess.run("command -v fastqjoin", shell=True, stdout=subprocess.DEVNULL).returncode == 0:
                join_cmd = "fastqjoin"
                
        cmd1 = f"{join_cmd} {cfg['input1r1dir']} {cfg['input1r2dir']} -o {output_prefix1}.fq >> {pipelinerrorlog} 2>&1"
        cmd2 = f"{join_cmd} {cfg['input2r1dir']} {cfg['input2r2dir']} -o {output_prefix2}.fq >> {pipelinerrorlog} 2>&1"
        print(f"Running: {cmd1}")
        subprocess.run(cmd1, shell=True)
        print(f"Running: {cmd2}")
        subprocess.run(cmd2, shell=True)

def step_filter(cfg, pipelinerrorlog):
    print("Executing Step: Filter")
    lib_name = 'VH' if int(cfg['libtype']) == 0 else 'VL'
    round_init_dir = f"{cfg['projectdir']}/{cfg['projectname']}/{lib_name}/InitialRound"
    round_final_dir = f"{cfg['projectdir']}/{cfg['projectname']}/{lib_name}/FinalRound"
    paired = int(cfg.get('pairedend', 0))
    
    # Resolve join/raw output paths
    if paired == 1:
        input1 = f"{round_init_dir}/{Path(cfg['input1r1dir']).stem}.fqjoin"
        input2 = f"{round_final_dir}/{Path(cfg['input2r1dir']).stem}.fqjoin"
    else:
        input1 = cfg['input1dir']
        input2 = cfg['input2dir']
        
    # Check if we need to parser ID ( Illumina IDs with whitespace )
    def needs_parser_id(fastq_path):
        if not os.path.exists(fastq_path):
            return False
        with open(fastq_path, 'r', encoding='utf-8', errors='ignore') as f:
            for _ in range(500):
                line = f.readline()
                if not line:
                    break
                if line.startswith('@') and (' ' in line or '\t' in line):
                    return True
        return False
        
    def run_parser_id(in_path, out_path):
        print(f"Sanitizing spaces in fastq headers: {in_path} -> {out_path}")
        # Strip trailing metadata in Illumina headers
        # e.g., @instrument:run:flowcell:lane:tile:x:y space metadata -> @instrument:run:flowcell:lane:tile:x:y
        pat = re.compile(r'^@([^:]+:\d+:[^:]+:\d+:\d+:\d+:\d+)\s+(.+)')
        with open(in_path, 'r', encoding='utf-8', errors='ignore') as f_in, open(out_path, 'w', encoding='utf-8') as f_out:
            for line in f_in:
                line = line.rstrip('\r\n')
                if line.startswith('@'):
                    m = pat.match(line)
                    if m:
                        line = '@' + m.group(1)
                f_out.write(line + '\n')

    # Resolve final fastq paths before prinseq
    input1_parsed = input1
    input2_parsed = input2
    
    if needs_parser_id(input1):
        input1_parsed = f"{round_init_dir}/{Path(input1).stem}newid.fq"
        run_parser_id(input1, input1_parsed)
    if needs_parser_id(input2):
        input2_parsed = f"{round_final_dir}/{Path(input2).stem}newid.fq"
        run_parser_id(input2, input2_parsed)
        
    # Build prinseq arguments
    minqual = cfg.get('minqual', '')
    minlen = cfg.get('minlen', '')
    
    prinseq_opts = ""
    if minqual:
        prinseq_opts += f" -min_qual_mean {minqual}"
    if minlen:
        prinseq_opts += f" -min_len {minlen}"
        
    # Prinseq without filter (beforefiltering)
    before1 = f"{round_init_dir}/{Path(input1).stem}beforefiltering"
    before2 = f"{round_final_dir}/{Path(input2).stem}beforefiltering"
    
    cmd_before1 = f"prinseq-lite -fastq {input1_parsed} -out_format 1 -out_bad null -out_good {before1} >> {pipelinerrorlog} 2>&1"
    cmd_before2 = f"prinseq-lite -fastq {input2_parsed} -out_format 1 -out_bad null -out_good {before2} >> {pipelinerrorlog} 2>&1"
    print(f"Running: {cmd_before1}")
    subprocess.run(cmd_before1, shell=True)
    print(f"Running: {cmd_before2}")
    subprocess.run(cmd_before2, shell=True)
    
    # Prinseq with filter (filtered)
    filtered1 = f"{round_init_dir}/{Path(input1).stem}filtered"
    filtered2 = f"{round_final_dir}/{Path(input2).stem}filtered"
    
    cmd_filt1 = f"prinseq-lite -fastq {input1_parsed} {prinseq_opts} -out_format 5 -out_bad null -out_good {filtered1} >> {pipelinerrorlog} 2>&1"
    cmd_filt2 = f"prinseq-lite -fastq {input2_parsed} {prinseq_opts} -out_format 5 -out_bad null -out_good {filtered2} >> {pipelinerrorlog} 2>&1"
    print(f"Running: {cmd_filt1}")
    subprocess.run(cmd_filt1, shell=True)
    print(f"Running: {cmd_filt2}")
    subprocess.run(cmd_filt2, shell=True)
    
    # FastQC on filtered sequences
    subprocess.run(f"fastqc -q -o {round_init_dir} {filtered1}.fastq >> {pipelinerrorlog} 2>&1", shell=True)
    subprocess.run(f"fastqc -q -o {round_final_dir} {filtered2}.fastq >> {pipelinerrorlog} 2>&1", shell=True)

def step_translate(cfg, pipelinerrorlog):
    print("Executing Step: Translate")
    lib_name = 'VH' if int(cfg['libtype']) == 0 else 'VL'
    round_init_dir = f"{cfg['projectdir']}/{cfg['projectname']}/{lib_name}/InitialRound"
    round_final_dir = f"{cfg['projectdir']}/{cfg['projectname']}/{lib_name}/FinalRound"
    paired = int(cfg.get('pairedend', 0))
    
    # Resolve filtered paths
    if paired == 1:
        stem1 = f"{Path(cfg['input1r1dir']).stem}.fqjoin"
        stem2 = f"{Path(cfg['input2r1dir']).stem}.fqjoin"
    else:
        stem1 = Path(cfg['input1dir']).stem
        stem2 = Path(cfg['input2dir']).stem
        
    filtered1_fasta = f"{round_init_dir}/{stem1}filtered.fasta"
    filtered2_fasta = f"{round_final_dir}/{stem2}filtered.fasta"
    
    gencode_path = f"{cfg['packagedir']}/data/genetic_code.fasta"
    
    # Outputs
    out_aa1 = f"{round_init_dir}/{stem1}aa.fasta"
    out_nt1 = f"{round_init_dir}/{stem1}nt.fasta"
    out_aa2 = f"{round_final_dir}/{stem2}aa.fasta"
    out_nt2 = f"{round_final_dir}/{stem2}nt.fasta"
    
    libtype = int(cfg['libtype'])
    
    print(f"Translating {filtered1_fasta} -> {out_aa1} / {out_nt1}")
    translate_all(filtered1_fasta, gencode_path, out_aa1, out_nt1, libtype)
    print(f"Translating {filtered2_fasta} -> {out_aa2} / {out_nt2}")
    translate_all(filtered2_fasta, gencode_path, out_aa2, out_nt2, libtype)

def translate_all(input_fasta, gencode_path, out_aa_fasta, out_nt_fasta, libtype):
    codon_table = load_codon_table(gencode_path)
    
    records = []
    with open(input_fasta, 'r', encoding='utf-8', errors='ignore') as f:
        curr_id = None
        curr_seq = []
        for line in f:
            line = line.rstrip('\r\n')
            if line.startswith('>'):
                if curr_id is not None:
                    records.append((curr_id, "".join(curr_seq)))
                curr_id = line[1:]
                curr_seq = []
            else:
                curr_seq.append(line)
        if curr_id is not None:
            records.append((curr_id, "".join(curr_seq)))
            
    with open(out_aa_fasta, 'w', encoding='utf-8') as f_aa, open(out_nt_fasta, 'w', encoding='utf-8') as f_nt:
        for seq_id, dna_seq in records:
            res = translate_read(seq_id, dna_seq, codon_table, libtype)
            if res:
                frame_str = res['frame_str']
                aa_seq = res['aa_seq']
                cdr_seq = res['cdr_seq']
                cdr_start = res['cdr_start']
                cdr_end = res['cdr_end']
                
                # Write to aa.fasta
                f_aa.write(f">{seq_id}|FRAME:{frame_str}\n")
                f_aa.write(f"{aa_seq}\n")
                f_aa.write(f"#{seq_id}|FRAME:{frame_str}|[{cdr_start}-{cdr_end}]\n")
                f_aa.write(f"{cdr_seq}*\n")
                
                # Write to nt.fasta
                f_nt.write(f">{seq_id}|FRAME:{frame_str}\n")
                f_nt.write(f"{dna_seq}\n")
                
        # Append end marker '<' expected by downstream frequency counter
        f_aa.write("<\n")

def step_frequency(cfg, pipelinerrorlog):
    print("Executing Step: Frequency")
    lib_name = 'VH' if int(cfg['libtype']) == 0 else 'VL'
    round_init_dir = f"{cfg['projectdir']}/{cfg['projectname']}/{lib_name}/InitialRound"
    round_final_dir = f"{cfg['projectdir']}/{cfg['projectname']}/{lib_name}/FinalRound"
    selected_dir = f"{cfg['projectdir']}/{cfg['projectname']}/{lib_name}/SelectedSequences"
    paired = int(cfg.get('pairedend', 0))
    
    if paired == 1:
        stem1 = f"{Path(cfg['input1r1dir']).stem}.fqjoin"
        stem2 = f"{Path(cfg['input2r1dir']).stem}.fqjoin"
    else:
        stem1 = Path(cfg['input1dir']).stem
        stem2 = Path(cfg['input2dir']).stem
        
    aa1 = f"{round_init_dir}/{stem1}aa.fasta"
    aa2 = f"{round_final_dir}/{stem2}aa.fasta"
    filt1 = f"{round_init_dir}/{stem1}filtered.fasta"
    filt2 = f"{round_final_dir}/{stem2}filtered.fasta"
    freq1 = f"{round_init_dir}/{stem1}freq.txt"
    freq2 = f"{round_final_dir}/{stem2}freq.txt"
    
    # 1. Frequency counters
    print(f"Calculating frequency for {aa1}")
    frequency_counter(aa1, filt1, freq1)
    print(f"Calculating frequency for {aa2}")
    frequency_counter(aa2, filt2, freq2)
    
    # 2. Find duplicates / Enriched
    aux_name = 'vh' if int(cfg['libtype']) == 0 else 'vl'
    enriched = f"{selected_dir}/{aux_name}enriched.fasta"
    print(f"Finding duplicates (enrichment): {freq1} + {freq2} -> {enriched}")
    find_duplicates(freq1, freq2, enriched)
    
    # 3. Get N sequences
    listsize = int(cfg['listsize'])
    enriched_n = f"{selected_dir}/{aux_name}list{listsize}.fasta"
    print(f"Selecting top {listsize} sequences -> {enriched_n}")
    get_nsequences(enriched, enriched_n, listsize)

def frequency_counter(input_protein_fasta, input_filtered_fasta, output_file):
    libsize = 0
    if os.path.exists(input_filtered_fasta):
        with open(input_filtered_fasta, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                if line.startswith('>'):
                    libsize += 1
                    
    unique_cdrs = {}
    curr_id = None
    curr_seq = ""
    curr_idcdr = None
    curr_cdr = ""
    state = 0 # 0: expecting id or <, 1: expecting seq, 2: expecting idcdr, 3: expecting cdr
    
    with open(input_protein_fasta, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            line = line.rstrip('\r\n')
            if not line:
                continue
            if line.startswith('>') or line == '<':
                if curr_cdr and curr_cdr.endswith('*'):
                    cdr_clean = curr_cdr[:-1]
                    if cdr_clean not in unique_cdrs:
                        unique_cdrs[cdr_clean] = {
                            'idcdr': curr_idcdr,
                            'ids': [],
                            'seqs': []
                        }
                    unique_cdrs[cdr_clean]['ids'].append(curr_id)
                    unique_cdrs[cdr_clean]['seqs'].append(curr_seq)
                
                if line == '<':
                    break
                curr_id = line
                curr_seq = ""
                curr_idcdr = None
                curr_cdr = ""
                state = 1
            elif state == 1:
                curr_seq = line
                state = 2
            elif line.startswith('#'):
                curr_idcdr = line[1:]
                state = 3
            elif state == 3:
                curr_cdr = line
                
    sorted_cdrs = sorted(unique_cdrs.items(), key=lambda x: len(x[1]['ids']), reverse=True)
    with open(output_file, 'w', encoding='utf-8') as f_out:
        for cdr_seq, info in sorted_cdrs:
            freq_raw = len(info['ids'])
            freq = (freq_raw / libsize) * 100000000000000 if libsize > 0 else 0.0
            f_out.write(f"#{info['idcdr']}|{libsize}|{freq}\n")
            f_out.write(f"*{cdr_seq}\n")
            for seq_id, seq in zip(info['ids'], info['seqs']):
                f_out.write(f"{seq_id}\n")
                f_out.write(f"{seq}\n")

def find_duplicates(input1, input2, output_file):
    cdrs_db = {}
    libsize1 = 1
    header_pat = re.compile(r'^#(.+)\|FRAME:(\d+)[+-]?\|([^|]+)\|(\d+)\|([\d.e+-]+)')
    
    with open(input1, 'r', encoding='utf-8', errors='ignore') as f1:
        curr_freq = 0.0
        for line in f1:
            line = line.rstrip('\r\n')
            m = header_pat.match(line)
            if m:
                libsize1 = int(m.group(4))
                curr_freq = float(m.group(5))
            elif line.startswith('*'):
                cdr_seq = line[1:]
                cdrs_db[cdr_seq] = curr_freq
                
    allcdrs = {}
    libsize2 = 1
    
    with open(input2, 'r', encoding='utf-8', errors='ignore') as f2:
        curr_freq = 0.0
        curr_cdr = None
        curr_idseq = None
        state = 0
        
        for line in f2:
            line = line.rstrip('\r\n')
            m = header_pat.match(line)
            if m:
                libsize2 = int(m.group(4))
                curr_freq = float(m.group(5))
                state = 0
            elif line.startswith('*'):
                cdr_seq = line[1:]
                curr_cdr = cdr_seq
                
                if cdr_seq in cdrs_db:
                    freq1 = cdrs_db[cdr_seq]
                else:
                    freq1 = (1 / libsize1) * 100000000000000
                    
                freq2 = curr_freq
                if freq2 > freq1:
                    dif = freq2 / freq1
                    dif_str = f"{dif:.4f}"
                    allcdrs[cdr_seq] = {
                        'dif': float(dif_str),
                        'freq1': freq1,
                        'freq2': freq2,
                        'id': '',
                        'seq': ''
                    }
                state = 1
            elif state == 1:
                if line.startswith('>'):
                    curr_idseq = line[1:]
                else:
                    seq = line
                    if curr_cdr in allcdrs:
                        if len(seq) > len(allcdrs[curr_cdr]['seq']):
                            allcdrs[curr_cdr]['id'] = curr_idseq
                            allcdrs[curr_cdr]['seq'] = seq
                            
    sorted_enriched = sorted(allcdrs.items(), key=lambda x: x[1]['dif'], reverse=True)
    with open(output_file, 'w', encoding='utf-8') as f_out:
        for cdr_seq, info in sorted_enriched:
            dif_str = f"{info['dif']:.4f}"
            f_out.write(f">{info['id']}|FOLD-CHANGE:{dif_str}|P1:{info['freq1']}|P2:{info['freq2']}\n")
            f_out.write(f"{info['seq']}\n")

def get_nsequences(input_fasta, output_fasta, n):
    sequences = []
    curr_id = None
    curr_seq = []
    with open(input_fasta, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            line = line.rstrip('\r\n')
            if line.startswith('>'):
                if curr_id is not None:
                    sequences.append((curr_id, "".join(curr_seq)))
                    if len(sequences) == n:
                        break
                curr_id = line[1:]
                curr_seq = []
            else:
                curr_seq.append(line)
        else:
            if curr_id is not None:
                sequences.append((curr_id, "".join(curr_seq)))
                
    with open(output_fasta, 'w', encoding='utf-8') as f_out:
        for seq_id, seq in sequences[:n]:
            f_out.write(f">{seq_id}\n{seq}\n")

def step_number(cfg, pipelinerrorlog):
    print("Executing Step: Number")
    lib_name = 'VH' if int(cfg['libtype']) == 0 else 'VL'
    selected_dir = f"{cfg['projectdir']}/{cfg['projectname']}/{lib_name}/SelectedSequences"
    listsize = int(cfg['listsize'])
    aux_name = 'vh' if int(cfg['libtype']) == 0 else 'vl'
    
    enriched_n = f"{selected_dir}/{aux_name}list{listsize}.fasta"
    numbered = f"{selected_dir}/{aux_name}list{listsize}numbered.txt"
    numbered_f = f"{selected_dir}/{aux_name}list{listsize}numbered.fasta"
    
    print(f"Numbering residues via UCL Abnum: {enriched_n} -> {numbered}")
    numberab(enriched_n, numbered)
    print(f"Converting Abnum alignment to fasta format: {numbered} -> {numbered_f}")
    convert_to_fasta(numbered, numbered_f)

def numberab(input_fasta, output_file):
    # Reads sequences, queries UCL Abnum, and writes Kabat alignments
    seqs = []
    with open(input_fasta, 'r', encoding='utf-8', errors='ignore') as f:
        curr_id = None
        curr_seq = []
        for line in f:
            line = line.rstrip('\r\n')
            if line.startswith('>'):
                if curr_id is not None:
                    seqs.append((curr_id, "".join(curr_seq)))
                # Extract before |P1
                # Format: >idseq|FRAME:1+|FOLD-CHANGE:12.3|P1:...
                # So strip everything starting with |P1
                header = line[1:]
                if '|P1' in header:
                    header = header.split('|P1')[0]
                curr_id = header
                curr_seq = []
            else:
                curr_seq.append(line)
        if curr_id is not None:
            seqs.append((curr_id, "".join(curr_seq)))
            
    # Clear output file first
    open(output_file, 'w', encoding='utf-8').close()
    
    with open(output_file, 'a', encoding='utf-8') as f_out:
        for header, seq in seqs:
            url = f"http://www.bioinf.org.uk/abs/abnum/abnum.cgi?plain=1&aaseq={seq}&scheme=-k"
            print(f"Querying Abnum for: {header}")
            try:
                req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
                with urllib.request.urlopen(req, timeout=30) as resp:
                    content = resp.read().decode('utf-8')
                if len(content) > 1:
                    f_out.write(f">{header}\n")
                    f_out.write(content)
            except Exception as e:
                print(f"WARNING: Abnum query failed for {header}: {e}")

def convert_to_fasta(input_file, output_file):
    cdrs = []
    curr_id = None
    curr_seq = []
    pat = re.compile(r'([A-Z]?)$')
    with open(input_file, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            line = line.rstrip('\r\n')
            if line.startswith('>'):
                if curr_id is not None:
                    cdrs.append((curr_id, "".join(curr_seq)))
                curr_id = line[1:]
                curr_seq = []
            else:
                m = pat.search(line)
                if m and m.group(1):
                    curr_seq.append(m.group(1))
        if curr_id is not None:
            cdrs.append((curr_id, "".join(curr_seq)))
            
    with open(output_file, 'w', encoding='utf-8') as f_out:
        for seq_id, seq in cdrs:
            f_out.write(f">{seq_id}\n{seq}\n")

def step_nt_recovery(cfg, pipelinerrorlog):
    print("Executing Step: NT-Recovery")
    lib_name = 'VH' if int(cfg['libtype']) == 0 else 'VL'
    selected_dir = f"{cfg['projectdir']}/{cfg['projectname']}/{lib_name}/SelectedSequences"
    round_final_dir = f"{cfg['projectdir']}/{cfg['projectname']}/{lib_name}/FinalRound"
    listsize = int(cfg['listsize'])
    aux_name = 'vh' if int(cfg['libtype']) == 0 else 'vl'
    paired = int(cfg.get('pairedend', 0))
    
    if paired == 1:
        stem2 = f"{Path(cfg['input2r1dir']).stem}.fqjoin"
    else:
        stem2 = Path(cfg['input2dir']).stem
        
    numbered_f = f"{selected_dir}/{aux_name}list{listsize}numbered.fasta"
    enriched_n = f"{selected_dir}/{aux_name}list{listsize}.fasta"
    nt_final = f"{round_final_dir}/{stem2}nt.fasta"
    list_nt = f"{selected_dir}/{aux_name}list{listsize}numberednt.fasta"
    
    print(f"Recovering original nucleotides: {numbered_f} -> {list_nt}")
    get_ntsequence(numbered_f, enriched_n, nt_final, list_nt)

def get_ntsequence(numbered_fasta, enriched_fasta, nt_fasta, output_fasta):
    seqs = {}
    header_pat1 = re.compile(r'^>([^|]+)\|FRAME:(\d+)([-+])?\|FOLD-CHANGE:([\d.]*)')
    
    with open(numbered_fasta, 'r', encoding='utf-8', errors='ignore') as f1:
        curr_id = None
        curr_frame = None
        curr_sense = '+'
        curr_fold = ''
        curr_seq = []
        for line in f1:
            line = line.rstrip('\r\n')
            if line.startswith('>'):
                if curr_id is not None:
                    seqs[curr_id] = {
                        'id': f"{curr_id}|FOLD-CHANGE:{curr_fold}",
                        'subcdr': "".join(curr_seq).replace('-', ''),
                        'frame': curr_frame,
                        'sense': curr_sense,
                        'protein': '',
                        'dna': '',
                        'cds': ''
                    }
                m = header_pat1.match(line)
                if m:
                    curr_id = m.group(1)
                    curr_frame = int(m.group(2))
                    curr_sense = m.group(3) if m.group(3) else '+'
                    curr_fold = m.group(4)
                    curr_seq = []
                else:
                    curr_id = None
            else:
                if curr_id is not None:
                    curr_seq.append(line)
        if curr_id is not None:
            seqs[curr_id] = {
                'id': f"{curr_id}|FOLD-CHANGE:{curr_fold}",
                'subcdr': "".join(curr_seq).replace('-', ''),
                'frame': curr_frame,
                'sense': curr_sense,
                'protein': '',
                'dna': '',
                'cds': ''
            }
            
    header_pat2 = re.compile(r'^>([^|]+)\|FRAME')
    with open(enriched_fasta, 'r', encoding='utf-8', errors='ignore') as f2:
        curr_id = None
        curr_seq = []
        for line in f2:
            line = line.rstrip('\r\n')
            if line.startswith('>'):
                if curr_id is not None and curr_id in seqs:
                    seqs[curr_id]['protein'] = "".join(curr_seq)
                m = header_pat2.match(line)
                if m:
                    curr_id = m.group(1)
                    curr_seq = []
                else:
                    curr_id = None
            else:
                if curr_id is not None:
                    curr_seq.append(line)
        if curr_id is not None and curr_id in seqs:
            seqs[curr_id]['protein'] = "".join(curr_seq)
            
    with open(nt_fasta, 'r', encoding='utf-8', errors='ignore') as f3:
        curr_id = None
        curr_seq = []
        for line in f3:
            line = line.rstrip('\r\n')
            if line.startswith('>'):
                if curr_id is not None and curr_id in seqs:
                    seqs[curr_id]['dna'] = "".join(curr_seq)
                m = header_pat2.match(line)
                if m:
                    curr_id = m.group(1)
                    curr_seq = []
                else:
                    curr_id = None
            else:
                if curr_id is not None:
                    curr_seq.append(line)
        if curr_id is not None and curr_id in seqs:
            seqs[curr_id]['dna'] = "".join(curr_seq)
            
    # Process slices
    for seq_id, info in seqs.items():
        if not info['protein'] or not info['dna']:
            continue
        subcdr = info['subcdr']
        protein = info['protein']
        start_idx = protein.find(subcdr)
        if start_idx != -1:
            nt_start = (start_idx * 3) + (info['frame'] - 1)
            nt_len = len(subcdr) * 3
            if info['sense'] == '+':
                info['cds'] = info['dna'][nt_start : nt_start + nt_len]
            else:
                rc_dna = reverse_complement(info['dna'])
                rc_cds = rc_dna[nt_start : nt_start + nt_len]
                info['cds'] = reverse_complement(rc_cds)
                
    with open(output_fasta, 'w', encoding='utf-8') as f_out:
        for seq_id, info in seqs.items():
            f_out.write(f">{info['id']}\n")
            cds = info['cds']
            if len(cds) > 1:
                for j in range(0, len(cds), 70):
                    f_out.write(cds[j : j + 70] + '\n')
            else:
                f_out.write("Could not retrieve nucleotide sequence.")
                if 'X' in info['protein']:
                    f_out.write(" Unknown base (N).\n")
                else:
                    f_out.write("\n")

def step_stats(cfg, pipelinerrorlog):
    print("Executing Step: Stats")
    lib_name = 'VH' if int(cfg['libtype']) == 0 else 'VL'
    round_init_dir = f"{cfg['projectdir']}/{cfg['projectname']}/{lib_name}/InitialRound"
    round_final_dir = f"{cfg['projectdir']}/{cfg['projectname']}/{lib_name}/FinalRound"
    selected_dir = f"{cfg['projectdir']}/{cfg['projectname']}/{lib_name}/SelectedSequences"
    report_dir = f"{cfg['projectdir']}/{cfg['projectname']}/Report"
    listsize = int(cfg['listsize'])
    aux_name = 'vh' if int(cfg['libtype']) == 0 else 'vl'
    paired = int(cfg.get('pairedend', 0))
    
    if paired == 1:
        stem1 = f"{Path(cfg['input1r1dir']).stem}.fqjoin"
        stem2 = f"{Path(cfg['input2r1dir']).stem}.fqjoin"
    else:
        stem1 = Path(cfg['input1dir']).stem
        stem2 = Path(cfg['input2dir']).stem
        
    filtered1 = f"{round_init_dir}/{stem1}filtered.fasta"
    filtered2 = f"{round_final_dir}/{stem2}filtered.fasta"
    
    enriched_n = f"{selected_dir}/{aux_name}list{listsize}.fasta"
    numbered_f = f"{selected_dir}/{aux_name}list{listsize}numbered.fasta"
    output_txt = f"{report_dir}/{aux_name}outputRstats.txt"
    
    # 1. Count filtered reads
    sizer0 = count_fasta_records(filtered1)
    sizern = count_fasta_records(filtered2)
    ncandidates = count_fasta_records(numbered_f)
    
    print(f"Calculating proportion z-test: sizer0={sizer0}, sizern={sizern}, ncandidates={ncandidates}")
    calculate_z_test(enriched_n, numbered_f, sizer0, sizern, ncandidates, output_txt)

def count_fasta_records(fasta_path):
    if not os.path.exists(fasta_path):
        return 0
    cnt = 0
    with open(fasta_path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if line.startswith('>'):
                cnt += 1
    return cnt

def calculate_z_test(enriched_fasta, abnum_fasta, sizer0, sizern, ncandidates, output_txt):
    seqs = {}
    count = 0
    header_pat = re.compile(r'^>([^|]+)\|FRAME')
    with open(abnum_fasta, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            line = line.rstrip('\r\n')
            m = header_pat.match(line)
            if m:
                seqs[m.group(1)] = {
                    'freq1': 0.0,
                    'freq2': 0.0,
                    'count': count
                }
                count += 1
                
    enriched_pat = re.compile(r'^>([^|]+)\|FRAME:\d[+-]?\|FOLD-CHANGE:[\d.]+\|P1:([\d.]+)\|P2:([\d.]+)')
    with open(enriched_fasta, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            line = line.rstrip('\r\n')
            m = enriched_pat.match(line)
            if m:
                seq_id = m.group(1)
                if seq_id in seqs:
                    seqs[seq_id]['freq1'] = float(m.group(2)) / 100000000000000.0
                    seqs[seq_id]['freq2'] = float(m.group(3)) / 100000000000000.0
                    
    def pnorm(z):
        return 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))
        
    with open(output_txt, 'w', encoding='utf-8') as f_out:
        if ncandidates > 0:
            alfa_bonf = 0.05 / ncandidates
            f_out.write(f"alfaBonf:\t{alfa_bonf}\n")
            
            sorted_seqs = sorted(seqs.items(), key=lambda x: x[1]['count'])
            for seq_id, info in sorted_seqs:
                f_out.write(f">{seq_id}\n")
                pantes = info['freq1']
                pdepois = info['freq2']
                pantes_menos_pdepois = pdepois - pantes
                n1 = sizer0
                n2 = sizern
                
                pantesdepois = (n1 * pantes + n2 * pdepois) / (n1 + n2)
                if pantesdepois <= 0.0 or pantesdepois >= 1.0:
                    erropadiff = 0.0
                else:
                    erropadiff = math.sqrt(pantesdepois * (1.0 - pantesdepois) * (1.0/n1 + 1.0/n2))
                    
                if erropadiff > 0:
                    z_obs = pantes_menos_pdepois / erropadiff
                    p_value = 1.0 - pnorm(z_obs)
                else:
                    p_value = 1.0
                    
                ep_diff_ic = math.sqrt((pantes * (1.0 - pantes) / n1) + (pdepois * (1.0 - pdepois) / n2))
                inf_ic = pantes - pdepois - 1.96 * ep_diff_ic
                sup_ic = pantes - pdepois + 1.96 * ep_diff_ic
                
                f_out.write(f"p-value:\t{p_value}\n")
                f_out.write(f"infIC:\t{inf_ic}\n")
                f_out.write(f"supIC:\t{sup_ic}\n")
        else:
            f_out.write("No candidates\n")

def step_germline(cfg, pipelinerrorlog):
    print("Executing Step: Germline")
    libtype = int(cfg['libtype'])
    lib_name = 'VH' if libtype == 0 else 'VL'
    selected_dir = f"{cfg['projectdir']}/{cfg['projectname']}/{lib_name}/SelectedSequences"
    listsize = int(cfg['listsize'])
    aux_name = 'vh' if libtype == 0 else 'vl'
    
    numbered_f = f"{selected_dir}/{aux_name}list{listsize}numbered.fasta"
    germclass = f"{selected_dir}/{aux_name}list{listsize}numberedgermlineclassification.txt"
    
    igblastdir = cfg['igblastdir']
    germlinedb = f"{cfg['packagedir']}/data/germline_human/{lib_name}_germline_human"

    igblastp_bin = shutil.which("igblastp") or os.path.join(igblastdir, "bin", "igblastp")
    igblast_install_dir = str(Path(igblastp_bin).parent.parent)
    cmd = f"{igblastp_bin} -germline_db_V {germlinedb} -query {numbered_f} -organism human -domain_system kabat -outfmt 7 > {germclass}"
    print(f"Running: {cmd}")
    env = os.environ.copy()
    env['IGDATA'] = igblast_install_dir
    subprocess.run(cmd, shell=True, check=True, env=env)

def step_report(cfg, pipelinerrorlog):
    print("Executing Step: Report")
    lib_name = 'VH' if int(cfg['libtype']) == 0 else 'VL'
    round_init_dir = f"{cfg['projectdir']}/{cfg['projectname']}/{lib_name}/InitialRound"
    round_final_dir = f"{cfg['projectdir']}/{cfg['projectname']}/{lib_name}/FinalRound"
    selected_dir = f"{cfg['projectdir']}/{cfg['projectname']}/{lib_name}/SelectedSequences"
    report_dir = f"{cfg['projectdir']}/{cfg['projectname']}/Report"
    paired = int(cfg.get('pairedend', 0))
    
    if paired == 1:
        stem1 = f"{Path(cfg['input1r1dir']).stem}.fqjoin"
        stem2 = f"{Path(cfg['input2r1dir']).stem}.fqjoin"
    else:
        stem1 = Path(cfg['input1dir']).stem
        stem2 = Path(cfg['input2dir']).stem
        
    before1 = f"{round_init_dir}/{stem1}beforefiltering.fasta"
    before2 = f"{round_final_dir}/{stem2}beforefiltering.fasta"
    csv_path = f"{cfg['projectdir']}/{cfg['projectname']}/{lib_name}/{lib_name.lower()}SequenceCounting.csv"
    
    plot1_path = f"{report_dir}/plot1{lib_name.lower()}.png"
    plot2_path = f"{report_dir}/plot2{lib_name.lower()}.png"
    readlen = int(cfg.get('minlen', 300))
    
    # 1. Plotting statistics (pure Python + fallback to R)
    print("Plotting statistics...")
    plot_statistics(before1, before2, csv_path, plot1_path, plot2_path, readlen, lib_name.lower(), report_dir)
    
    # Clean up beforefiltering fasta files to save space
    for f in [before1, before2]:
        if os.path.exists(f):
            os.remove(f)

def plot_statistics(initial_fasta, final_fasta, csv_path, plot1_path, plot2_path, readlen, libtype, report_dir):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import numpy as np
        
        # Calculate good/bad reads for initial/final fasta
        def get_good_bad(fasta_path):
            total = 0
            good = 0
            if not os.path.exists(fasta_path):
                return 0, 0, 0
            with open(fasta_path, 'r', encoding='utf-8', errors='ignore') as f:
                curr_len = 0
                for line in f:
                    line = line.rstrip('\r\n')
                    if line.startswith('>'):
                        if total > 0:
                            if curr_len >= readlen:
                                good += 1
                        total += 1
                        curr_len = 0
                    else:
                        curr_len += len(line)
                if total > 0:
                    if curr_len >= readlen:
                        good += 1
            bad = total - good
            return total, good, bad
            
        t_init, g_init, b_init = get_good_bad(initial_fasta)
        t_final, g_final, b_final = get_good_bad(final_fasta)
        
        # Plot 1: Stacked Bar Chart
        libraries = ['Initial', 'Final']
        adequate = [g_init / t_init if t_init > 0 else 0.0, g_final / t_final if t_final > 0 else 0.0]
        inadequate = [b_init / t_init if t_init > 0 else 0.0, b_final / t_final if t_final > 0 else 0.0]
        
        fig, ax = plt.subplots(figsize=(9, 7))
        width = 0.35
        ax.bar(libraries, inadequate, width, label='Inadequate', color='#ff7f0e')
        ax.bar(libraries, adequate, width, bottom=inadequate, label='Adequate', color='#1f77b4')
        ax.set_ylabel('Reads')
        ax.set_xlabel('Library')
        ax.set_title('Proportion of Reads with Adequate Length')
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: '{:.0%}'.format(y)))
        ax.legend(title='Read Length')
        plt.tight_layout()
        plt.savefig(plot1_path, dpi=300)
        plt.close()
        
        # Plot 2: Total Reads by Task
        counts = {}
        with open(csv_path, 'r', encoding='utf-8') as f:
            next(f) # skip header
            for line in f:
                line = line.rstrip('\r\n')
                if not line:
                    continue
                lib, val, step = line.split(',')
                val = int(val)
                if step not in counts:
                    counts[step] = {}
                counts[step][lib] = val
                
        step_mapping = [
            ('raw', 'None'),
            ('joining', 'Joining'),
            ('filtering', 'Filtering'),
            ('translation', 'Translation'),
            ('frequency', 'Frequency_Calculation'),
            ('enrichment', 'Enrichment')
        ]
        
        plot_steps = []
        r0_vals = []
        rn_vals = []
        
        for step_key, step_name in step_mapping:
            if step_key in counts or step_key == 'enrichment':
                plot_steps.append(step_name)
                r0_vals.append(counts.get(step_key, {}).get('R0', 0))
                if step_key == 'enrichment':
                    rn_vals.append(counts.get(step_key, {}).get('Selected', 0))
                else:
                    rn_vals.append(counts.get(step_key, {}).get('RN', 0))
                    
        x = np.arange(len(plot_steps))
        width = 0.35
        
        fig, ax = plt.subplots(figsize=(9, 7))
        ax.bar(x - width/2, r0_vals, width, label='Initial', color='#1f77b4')
        ax.bar(x + width/2, rn_vals, width, label='Final', color='#ff7f0e')
        ax.set_xlabel('Task')
        ax.set_ylabel('Reads')
        ax.set_title('Total Number of Reads by Task')
        ax.set_xticks(x)
        ax.set_xticklabels(plot_steps, rotation=15)
        ax.legend(title='Library')
        ax.ticklabel_format(style='plain', axis='y')
        
        plt.tight_layout()
        plt.savefig(plot2_path, dpi=300)
        plt.close()
        print("Plots generated via Matplotlib successfully.")
        
    except Exception as e:
        print(f"Matplotlib plotting skipped/failed: {e}. Generating R scripts as fallback...")
        # Fallback to creating and running R script
        scriptlen = f"{report_dir}/plot1{libtype}.r"
        scriptask = f"{report_dir}/plot2{libtype}.r"
        
        # Write R plot 1
        with open(scriptlen, 'w', encoding='utf-8') as f_r1:
            f_r1.write("#! /usr/bin/env Rscript\n")
            f_r1.write("library(\"ggplot2\")\n")
            f_r1.write("library(\"scales\")\n")
            f_r1.write("input <- read.table(header= TRUE, stringsAsFactors = FALSE, text = '\n")
            f_r1.write("Library\t\t\t\t\tLength\t\t\t\t\t\treads\n")
            
            def get_good_bad_raw(fasta_path):
                total = 0
                good = 0
                if not os.path.exists(fasta_path):
                    return 0, 0
                with open(fasta_path, 'r', encoding='utf-8', errors='ignore') as f:
                    curr_len = 0
                    for line in f:
                        line = line.rstrip('\r\n')
                        if line.startswith('>'):
                            if total > 0:
                                if curr_len >= readlen:
                                    good += 1
                            total += 1
                            curr_len = 0
                        else:
                            curr_len += len(line)
                    if total > 0:
                        if curr_len >= readlen:
                            good += 1
                return total, good
                
            t_init, g_init = get_good_bad_raw(initial_fasta)
            t_final, g_final = get_good_bad_raw(final_fasta)
            
            pg_init = g_init / t_init if t_init > 0 else 0
            pb_init = 1 - pg_init
            pg_final = g_final / t_final if t_final > 0 else 0
            pb_final = 1 - pg_final
            
            f_r1.write(f"Initial\t\t\t\t\tInadequate\t\t\t\t\t{pb_init:.4f}\n")
            f_r1.write(f"Initial\t\t\t\t\tAdequate\t\t\t\t\t{pg_init:.4f}\n")
            f_r1.write(f"Final\t\t\t\t\tInadequate\t\t\t\t\t{pb_final:.4f}\n")
            f_r1.write(f"Final\t\t\t\t\tAdequate\t\t\t\t\t{pg_final:.4f}\n")
            f_r1.write("')\n")
            f_r1.write(f"png(file=\"{plot1_path}\", res=300, width=7, height=7, units = 'in')\n")
            f_r1.write("m <- ggplot(data = input,aes(x= factor(Library, levels = unique(Library)),y=reads,fill=Length)) + geom_bar(width=0.3,stat=\"identity\") + scale_y_continuous(labels=percent)\n")
            f_r1.write("m <- m +  scale_fill_discrete(name=\"Read Length\")\n")
            f_r1.write("m <- m + xlab(\"Library\") + ylab(\"Reads\")\n")
            f_r1.write("m <- m + ggtitle(\"Proportion of Reads with Adequate Length\")\n")
            f_r1.write("print(m)\n")
            f_r1.write("garbage <- dev.off()\n")
            
        # Write R plot 2
        with open(scriptask, 'w', encoding='utf-8') as f_r2:
            f_r2.write("#! /usr/bin/env Rscript\n")
            f_r2.write("library(\"ggplot2\")\n")
            f_r2.write("library(\"scales\")\n")
            f_r2.write("input <- read.table(header= TRUE, stringsAsFactors = FALSE, text = '\n")
            f_r2.write("Library\t\t\t\t\tProcess\t\t\t\t\t\treads\n")
            
            with open(csv_path, 'r', encoding='utf-8') as f_csv:
                next(f_csv) # skip header
                for line in f_csv:
                    line = line.rstrip('\r\n')
                    if not line:
                        continue
                    lib, val, step = line.split(',')
                    if step == 'raw':
                        if lib == 'R0': f_r2.write(f"Initial\t\t\t\t\tNone\t\t\t\t\t{val}\n")
                        elif lib == 'RN': f_r2.write(f"Final\t\t\t\t\tNone\t\t\t\t\t{val}\n")
                    elif step == 'joining':
                        if lib == 'R0': f_r2.write(f"Initial\t\t\t\t\tJoining\t\t\t\t\t{val}\n")
                        elif lib == 'RN': f_r2.write(f"Final\t\t\t\t\tJoining\t\t\t\t\t{val}\n")
                    elif step == 'filtering':
                        if lib == 'R0': f_r2.write(f"Initial\t\t\t\t\tFiltering\t\t\t\t\t{val}\n")
                        elif lib == 'RN': f_r2.write(f"Final\t\t\t\t\tFiltering\t\t\t\t\t{val}\n")
                    elif step == 'translation':
                        if lib == 'R0': f_r2.write(f"Initial\t\t\t\t\tTranslation\t\t\t\t\t{val}\n")
                        elif lib == 'RN': f_r2.write(f"Final\t\t\t\t\tTranslation\t\t\t\t\t{val}\n")
                    elif step == 'frequency':
                        if lib == 'R0': f_r2.write(f"Initial\t\t\t\t\tFrequency_Calculation\t\t\t\t\t{val}\n")
                        elif lib == 'RN': f_r2.write(f"Final\t\t\t\t\tFrequency_Calculation\t\t\t\t\t{val}\n")
                    elif step == 'enrichment':
                        if lib == 'Selected': f_r2.write(f"Final\t\t\t\t\tEnrichment\t\t\t\t\t{val}\n")
            f_r2.write("')\n")
            f_r2.write(f"png(file=\"{plot2_path}\", res=300, width=9, height=7, units = 'in')\n")
            f_r2.write("m <- ggplot(data = input,aes(x= factor(Process, levels = unique(Process)),y=reads,fill=Library)) + geom_bar(stat=\"identity\",position=position_dodge())\n")
            f_r2.write("m <- m + xlab(\"Task\") + ylab(\"Reads\")\n")
            f_r2.write("m <- m + ggtitle(\"Total Number of Reads by Task\")\n")
            f_r2.write("options(scipen=999)\n")
            f_r2.write("print(m)\n")
            f_r2.write("garbage <- dev.off()\n")
            
        # Run scripts
        subprocess.run(f"Rscript {scriptlen}", shell=True)
        subprocess.run(f"Rscript {scriptask}", shell=True)

# =====================================================================
# Database / Counting logic
# =====================================================================

def count_sequences_step(input_file, is_fastq=True):
    if not os.path.exists(input_file):
        return 0
    cnt = 0
    with open(input_file, 'r', encoding='utf-8', errors='ignore') as f:
        if is_fastq:
            for _ in f:
                cnt += 1
            return cnt // 4
        else:
            for line in f:
                if line.startswith('>'):
                    cnt += 1
            return cnt

def count_cdr_headers(input_file):
    if not os.path.exists(input_file):
        return 0
    cnt = 0
    with open(input_file, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if line.startswith('#'):
                cnt += 1
    return cnt

def output_csv_counting(cfg):
    libtype = int(cfg['libtype'])
    lib_name = 'VH' if libtype == 0 else 'VL'
    round_init_dir = f"{cfg['projectdir']}/{cfg['projectname']}/{lib_name}/InitialRound"
    round_final_dir = f"{cfg['projectdir']}/{cfg['projectname']}/{lib_name}/FinalRound"
    selected_dir = f"{cfg['projectdir']}/{cfg['projectname']}/{lib_name}/SelectedSequences"
    paired = int(cfg.get('pairedend', 0))
    listsize = int(cfg['listsize'])
    aux_name = 'vh' if libtype == 0 else 'vl'
    
    if paired == 1:
        stem1 = f"{Path(cfg['input1r1dir']).stem}.fqjoin"
        stem2 = f"{Path(cfg['input2r1dir']).stem}.fqjoin"
    else:
        stem1 = Path(cfg['input1dir']).stem
        stem2 = Path(cfg['input2dir']).stem
        
    raw_init = cfg['input1r1dir'] if paired == 1 else cfg['input1dir']
    raw_final = cfg['input2r1dir'] if paired == 1 else cfg['input2dir']
    
    counts = {}
    # step raw
    counts['raw'] = {
        'R0': count_sequences_step(raw_init, is_fastq=True),
        'RN': count_sequences_step(raw_final, is_fastq=True)
    }
    
    # step joining
    if paired == 1:
        counts['joining'] = {
            'R0': count_sequences_step(f"{round_init_dir}/{stem1}.fq", is_fastq=True),
            'RN': count_sequences_step(f"{round_final_dir}/{stem2}.fq", is_fastq=True)
        }
        
    # step filtering
    counts['filtering'] = {
        'R0': count_sequences_step(f"{round_init_dir}/{stem1}filtered.fasta", is_fastq=False),
        'RN': count_sequences_step(f"{round_final_dir}/{stem2}filtered.fasta", is_fastq=False)
    }
    
    # step translation
    counts['translation'] = {
        'R0': count_sequences_step(f"{round_init_dir}/{stem1}aa.fasta", is_fastq=False),
        'RN': count_sequences_step(f"{round_final_dir}/{stem2}aa.fasta", is_fastq=False)
    }
    
    # step frequency
    counts['frequency'] = {
        'R0': count_cdr_headers(f"{round_init_dir}/{stem1}freq.txt"),
        'RN': count_cdr_headers(f"{round_final_dir}/{stem2}freq.txt")
    }
    
    # step enrichment
    counts['enrichment'] = {
        'Selected': count_sequences_step(f"{selected_dir}/{aux_name}enriched.fasta", is_fastq=False)
    }
    
    # step numeration
    counts['numeration'] = {
        'Selected': count_sequences_step(f"{selected_dir}/{aux_name}list{listsize}numbered.fasta", is_fastq=False)
    }
    
    csv_file = f"{cfg['projectdir']}/{cfg['projectname']}/{lib_name}/{lib_name.lower()}SequenceCounting.csv"
    with open(csv_file, 'w', encoding='utf-8') as f_out:
        f_out.write("library,reads,step\n")
        
        # Sort keys based on pipeline order index
        order = ['raw', 'joining', 'filtering', 'translation', 'frequency', 'enrichment', 'numeration']
        for step in order:
            if step not in counts:
                continue
            step_counts = counts[step]
            if step in ['enrichment', 'numeration']:
                f_out.write(f"Selected,{step_counts.get('Selected', 0)},{step}\n")
            else:
                f_out.write(f"R0,{step_counts.get('R0', 0)},{step}\n")
                f_out.write(f"RN,{step_counts.get('RN', 0)},{step}\n")
                
    print(f"Created sequence counts CSV: {csv_file}")

# =====================================================================
# Main Orchestrator Logic
# =====================================================================

def parse_cfg(config_file):
    cfg = {}
    if not os.path.exists(config_file):
        raise FileNotFoundError(f"Configuration file not found: {config_file}")
    with open(config_file, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.rstrip('\r\n')
            if not line or line.startswith('#') or line.startswith('-'):
                continue
            if ':' in line:
                key, val = line.split(':', 1)
                cfg[key.strip()] = val.strip()
    return cfg

def main():
    parser = argparse.ArgumentParser(description="ATTILA Immunoglobulin Sequence Analysis Pipeline")
    parser.add_argument('config', help="Path to the configuration file (.cfg)")
    parser.add_argument('log', help="Path to the error/log output file")
    
    args = parser.parse_args()
    
    # Load config
    cfg = parse_cfg(args.config)
    
    # Create required project directories
    libtype = int(cfg['libtype'])
    lib_name = 'VH' if libtype == 0 else 'VL'
    proj_dir = cfg['projectdir']
    proj_name = cfg['projectname']
    
    vh_vl_dir = f"{proj_dir}/{proj_name}/{lib_name}"
    round_init_dir = f"{vh_vl_dir}/InitialRound"
    round_final_dir = f"{vh_vl_dir}/FinalRound"
    selected_dir = f"{vh_vl_dir}/SelectedSequences"
    report_dir = f"{proj_dir}/{proj_name}/Report"
    
    for d in [round_init_dir, round_final_dir, selected_dir, report_dir]:
        os.makedirs(d, exist_ok=True)
        
    pipelinerrorlog = f"{proj_dir}/{proj_name}/{lib_name.lower()}PipelineError.log"
    # Touch error log
    open(pipelinerrorlog, 'a').close()
    
    # Redirect stderr to project log
    sys.stderr = open(args.log, 'a')
    
    # Identify steps to run
    steps_env = os.environ.get('ATTILA_STEPS', 'all').lower()
    run_all = (steps_env == 'all')
    enabled_steps = [s.strip() for s in steps_env.split(',') if s.strip()]
    
    def should_run(step_name):
        return run_all or (step_name in enabled_steps)
        
    # Execute Steps
    if should_run('join'):
        step_join(cfg, pipelinerrorlog)
        
    if should_run('filter'):
        step_filter(cfg, pipelinerrorlog)
        
    if should_run('translate'):
        step_translate(cfg, pipelinerrorlog)
        
    if should_run('frequency'):
        step_frequency(cfg, pipelinerrorlog)
        
    if should_run('number'):
        step_number(cfg, pipelinerrorlog)
        
    if should_run('nt-recovery'):
        step_nt_recovery(cfg, pipelinerrorlog)
        
    if should_run('stats'):
        step_stats(cfg, pipelinerrorlog)
        
    if should_run('germline'):
        step_germline(cfg, pipelinerrorlog)
        
    # Output counts CSV (happens always after frequency/number stages to ensure files exist)
    output_csv_counting(cfg)
    
    if should_run('report'):
        step_report(cfg, pipelinerrorlog)

if __name__ == '__main__':
    main()
