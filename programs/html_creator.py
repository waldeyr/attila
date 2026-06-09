#!/usr/bin/env python3
"""
html_creator.py - Translates html_creator.pl into native Python 3.
Generates the final Report.html summarizing candidate clones, proportion of reads, and germline classification.
"""

import sys
import os
import re
from pathlib import Path

def parse_fasta(input_fasta, table1, table2, counter):
    # Match ID, Frame, and Fold Change
    # Regex matching: ^>(.+)\|FRAME:(\d)[+-]?\|FOLD-CHANGE:(\d+\.\d+)$
    pat = re.compile(r'^>([^|]+)\|FRAME:(\d+)[+-]?\|FOLD-CHANGE:([\d.]+)')
    if not os.path.exists(input_fasta):
        return counter
        
    with open(input_fasta, 'r', encoding='utf-8', errors='ignore') as f:
        curr_id = None
        curr_seq = []
        for line in f:
            line = line.rstrip('\r\n')
            if line.startswith('>'):
                if curr_id is not None:
                    table2[curr_id]['seq'] = "".join(curr_seq).replace('-', '')
                m = pat.match(line)
                if m:
                    curr_id = m.group(1)
                    table1[curr_id] = {
                        'fc': f"{float(m.group(3)):.2f}",
                        'count': counter,
                        'germline': '',
                        'identity': 0.0
                    }
                    table2[curr_id] = {
                        'count': counter,
                        'seq': '',
                        'fr1': '',
                        'cdr1': '',
                        'fr2': '',
                        'cdr2': '',
                        'fr3': '',
                        'cdr3': '',
                        'fr4': ''
                    }
                    counter += 1
                    curr_seq = []
                else:
                    curr_id = None
            else:
                if curr_id is not None:
                    curr_seq.append(line)
        if curr_id is not None:
            table2[curr_id]['seq'] = "".join(curr_seq).replace('-', '')
            
    return counter

def parse_igblast(input_igblast, table1, table2):
    if not os.path.exists(input_igblast):
        return
        
    # Match query ID in IgBlast format 7
    # e.g., # Query: query_id|FRAME:2+|FOLD-CHANGE:12.34
    query_pat = re.compile(r'^#\s+Query:\s+([^|]+)\|FRAME:(\d+)[+-]?\|FOLD-CHANGE:([\d.]+)')
    
    # Match region lines: FR1/CDR1/FR2/CDR2/FR3
    # e.g., FR1  1  25  ...
    region_pat = re.compile(r'^(FR1|CDR1|FR2|CDR2|FR3)\s+(\d+)\s+(\d+)')
    
    # Match V gene classification (germline and identity)
    # e.g., V   IGHV3-23*01   98.2   ...
    v_pat = re.compile(r'^V.*?FOLD-CHANGE:[\d.]+\s+(\S+)?\s+([\d.]+)?\s+(\S+)')
    
    with open(input_igblast, 'r', encoding='utf-8', errors='ignore') as f:
        curr_id = None
        flag = -1
        for line in f:
            line = line.rstrip('\r\n')
            m_q = query_pat.match(line)
            if m_q:
                curr_id = m_q.group(1)
                flag = 0
                continue
                
            if not curr_id or curr_id not in table2:
                continue
                
            m_r = region_pat.match(line)
            if m_r:
                region_name = m_r.group(1).lower()
                start = int(m_r.group(2)) - 1
                end = int(m_r.group(3))
                seq = table2[curr_id]['seq']
                table2[curr_id][region_name] = seq[start : end]
                
                # If FR3, compute CDR3 and FR4
                if region_name == 'fr3':
                    # Search for WG\wG or FG\wG motif
                    matches = list(re.finditer(r'WG\wG|FG\wG', seq))
                    if matches:
                        cdr3end = matches[-1].start()
                        table2[curr_id]['cdr3'] = seq[end : cdr3end]
                        table2[curr_id]['fr4'] = seq[cdr3end :]
            else:
                if flag == 0:
                    m_v = v_pat.match(line)
                    if m_v:
                        table1[curr_id]['germline'] = m_v.group(1) if m_v.group(1) else ''
                        table1[curr_id]['identity'] = float(m_v.group(2)) if m_v.group(2) else 0.0
                        flag = 1

def parse_csv_loss(csv_file, aux_name, loss_dict):
    if not os.path.exists(csv_file):
        return
    loss_dict[aux_name] = {'r0': {'raw': 0, 'filtering': 0}, 'rn': {'raw': 0, 'filtering': 0}}
    with open(csv_file, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.rstrip('\r\n')
            if line.startswith('R0,') and line.endswith(',raw'):
                loss_dict[aux_name]['r0']['raw'] = int(line.split(',')[1])
            elif line.startswith('R0,') and line.endswith(',filtering'):
                loss_dict[aux_name]['r0']['filtering'] = int(line.split(',')[1])
            elif line.startswith('RN,') and line.endswith(',raw'):
                loss_dict[aux_name]['rn']['raw'] = int(line.split(',')[1])
            elif line.startswith('RN,') and line.endswith(',filtering'):
                loss_dict[aux_name]['rn']['filtering'] = int(line.split(',')[1])

def parse_stats(stat_file, stat_dict):
    if not os.path.exists(stat_file):
        return '0.05'
    alfa_bonf = '0.05'
    curr_id = None
    with open(stat_file, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            line = line.rstrip('\r\n')
            if line.startswith('alfaBonf:'):
                alfa_bonf = line.split('\t')[1].strip()
            elif line.startswith('>'):
                curr_id = line[1:].strip()
                stat_dict[curr_id] = {'pvalue': '', 'infIC': '', 'supIC': ''}
            elif curr_id:
                if line.startswith('p-value:'):
                    stat_dict[curr_id]['pvalue'] = line.split('\t')[1].strip()
                elif line.startswith('infIC:'):
                    stat_dict[curr_id]['infIC'] = line.split('\t')[1].strip()
                elif line.startswith('supIC:'):
                    stat_dict[curr_id]['supIC'] = line.split('\t')[1].strip()
    return alfa_bonf

def main():
    if len(sys.argv) < 14:
        print("Usage: python3 html_creator.py <in1> <in2> <in3> <in4> <out> <dir1> <dir2> <dir3> <dir4> <csvfilevh> <csvfilevl> <statfilevh> <statfilevl>")
        sys.exit(1)
        
    in1, in2, in3, in4, out, dir1, dir2, dir3, dir4, csvfilevh, csvfilevl, statfilevh, statfilevl = sys.argv[1:14]

    # Detect log files relative to project root (= dirname of index.html)
    project_root = os.path.dirname(out)
    known_logs = [
        ("vherror.log",         "VH Error Log"),
        ("vlerror.log",         "VL Error Log"),
        ("vhPipelineError.log", "VH Pipeline Error Log"),
        ("vlPipelineError.log", "VL Pipeline Error Log"),
        ("Report/webpage.log",  "Report Generation Log"),
    ]
    log_entries = [
        (rel, label) for rel, label in known_logs
        if os.path.exists(os.path.join(project_root, rel))
    ]

    # 1. Parse Fastas
    table1 = {}
    table2 = {}
    
    countvh = parse_fasta(in1, table1, table2, 0)
    countvl = parse_fasta(in2, table1, table2, countvh)
    
    # 2. Parse IgBlast
    parse_igblast(in3, table1, table2)
    parse_igblast(in4, table1, table2)
    
    # 3. Parse CSVs for loss of reads
    loss_dict = {}
    parse_csv_loss(csvfilevh, 'vh', loss_dict)
    parse_csv_loss(csvfilevl, 'vl', loss_dict)
    
    # Loss percentage calculations
    def calc_loss_pct(raw, filt):
        if raw <= 0: return "0.00"
        return f"{((raw - filt) / raw) * 100:.2f}"
        
    perdavhr0 = calc_loss_pct(loss_dict.get('vh', {}).get('r0', {}).get('raw', 0), loss_dict.get('vh', {}).get('r0', {}).get('filtering', 0))
    perdavhrn = calc_loss_pct(loss_dict.get('vh', {}).get('rn', {}).get('raw', 0), loss_dict.get('vh', {}).get('rn', {}).get('filtering', 0))
    perdavlr0 = calc_loss_pct(loss_dict.get('vl', {}).get('r0', {}).get('raw', 0), loss_dict.get('vl', {}).get('r0', {}).get('filtering', 0))
    perdavlrn = calc_loss_pct(loss_dict.get('vl', {}).get('rn', {}).get('raw', 0), loss_dict.get('vl', {}).get('rn', {}).get('filtering', 0))
    
    # 4. Parse Stats
    stat_dict = {}
    alfavh = parse_stats(statfilevh, stat_dict)
    alfavl = parse_stats(statfilevl, stat_dict)
    
    # 5. Write index.html at project root
    os.makedirs(project_root, exist_ok=True)
    with open(out, 'w', encoding='utf-8') as f_out:
        f_out.write("<!DOCTYPE html>\n")
        f_out.write("<html>\n")
        f_out.write("<head>\n")
        f_out.write("<meta charset=\"utf-8\">\n")
        f_out.write("<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n")
        f_out.write("<link href=\"https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css\" rel=\"stylesheet\">\n")
        f_out.write("<link href=\"https://cdn.jsdelivr.net/npm/bootstrap-icons@1.10.5/font/bootstrap-icons.css\" rel=\"stylesheet\">\n")
        f_out.write("<link href=\"https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&display=swap\" rel=\"stylesheet\">\n")
        f_out.write("<script src=\"https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/js/bootstrap.bundle.min.js\"></script>\n")
        f_out.write("<title>Selected sequences - ATTILA Report</title>\n")
        f_out.write("<style>\n")
        f_out.write(":root {\n")
        f_out.write("    --primary-color: #0f172a;\n")
        f_out.write("    --primary-light: #1e293b;\n")
        f_out.write("    --accent-color: #2563eb;\n")
        f_out.write("    --bg-light: #f8fafc;\n")
        f_out.write("    --border-color: #e2e8f0;\n")
        f_out.write("    /* Sequence Region Colors */\n")
        f_out.write("    --color-fr: #334155;\n")
        f_out.write("    --bg-fr: #f1f5f9;\n")
        f_out.write("    --color-cdr1: #0369a1;\n")
        f_out.write("    --bg-cdr1: #e0f2fe;\n")
        f_out.write("    --color-cdr2: #15803d;\n")
        f_out.write("    --bg-cdr2: #dcfce7;\n")
        f_out.write("    --color-cdr3: #7e22ce;\n")
        f_out.write("    --bg-cdr3: #f3e8ff;\n")
        f_out.write("}\n")
        f_out.write("body {\n")
        f_out.write("    font-family: 'Inter', -apple-system, BlinkMacSystemFont, \"Segoe UI\", Roboto, sans-serif;\n")
        f_out.write("    background-color: var(--bg-light);\n")
        f_out.write("    color: var(--primary-color);\n")
        f_out.write("    padding-top: 1.5rem;\n")
        f_out.write("    padding-bottom: 2rem;\n")
        f_out.write("}\n")
        f_out.write(".container-fluid {\n")
        f_out.write("    max-width: 1400px;\n")
        f_out.write("    padding: 0 1.5rem;\n")
        f_out.write("}\n")
        f_out.write("header {\n")
        f_out.write("    background: linear-gradient(135deg, #0f172a 0%, #1e293b 100%);\n")
        f_out.write("    border-radius: 12px;\n")
        f_out.write("    padding: 2rem;\n")
        f_out.write("    box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.1);\n")
        f_out.write("    margin-bottom: 2rem;\n")
        f_out.write("}\n")
        f_out.write("#application {\n")
        f_out.write("    display: flex;\n")
        f_out.write("    justify-content: center;\n")
        f_out.write("    align-items: center;\n")
        f_out.write("}\n")
        f_out.write("#application img {\n")
        f_out.write("    max-height: 80px;\n")
        f_out.write("    object-fit: contain;\n")
        f_out.write("    filter: drop-shadow(0 4px 6px rgba(0,0,0,0.2));\n")
        f_out.write("}\n")
        f_out.write(".nav-tabs {\n")
        f_out.write("    border-bottom: 2px solid var(--border-color);\n")
        f_out.write("    margin-bottom: 2rem;\n")
        f_out.write("    justify-content: center;\n")
        f_out.write("}\n")
        f_out.write(".nav-tabs .nav-link {\n")
        f_out.write("    color: #64748b;\n")
        f_out.write("    font-weight: 500;\n")
        f_out.write("    border: none !important;\n")
        f_out.write("    padding: 0.75rem 1.5rem;\n")
        f_out.write("    position: relative;\n")
        f_out.write("    transition: color 0.2s ease;\n")
        f_out.write("}\n")
        f_out.write(".nav-tabs .nav-link:hover {\n")
        f_out.write("    color: var(--accent-color);\n")
        f_out.write("}\n")
        f_out.write(".nav-tabs .nav-link.active {\n")
        f_out.write("    color: var(--accent-color) !important;\n")
        f_out.write("    background-color: transparent !important;\n")
        f_out.write("}\n")
        f_out.write(".nav-tabs .nav-link.active::after {\n")
        f_out.write("    content: '';\n")
        f_out.write("    position: absolute;\n")
        f_out.write("    bottom: -2px;\n")
        f_out.write("    left: 0;\n")
        f_out.write("    width: 100%;\n")
        f_out.write("    height: 3px;\n")
        f_out.write("    background-color: var(--accent-color);\n")
        f_out.write("    border-radius: 3px;\n")
        f_out.write("}\n")
        f_out.write(".card-section {\n")
        f_out.write("    background-color: #ffffff;\n")
        f_out.write("    border-radius: 12px;\n")
        f_out.write("    border: 1px solid var(--border-color);\n")
        f_out.write("    padding: 2rem;\n")
        f_out.write("    margin-bottom: 2rem;\n")
        f_out.write("    box-shadow: 0 1px 3px 0 rgba(0, 0, 0, 0.05);\n")
        f_out.write("    text-align: left;\n")
        f_out.write("}\n")
        f_out.write(".section-header {\n")
        f_out.write("    font-size: 1.25rem;\n")
        f_out.write("    font-weight: 600;\n")
        f_out.write("    color: var(--primary-color);\n")
        f_out.write("    margin-bottom: 1.5rem;\n")
        f_out.write("    border-bottom: 2px solid var(--border-color);\n")
        f_out.write("    padding-bottom: 0.5rem;\n")
        f_out.write("    display: flex;\n")
        f_out.write("    justify-content: space-between;\n")
        f_out.write("    align-items: center;\n")
        f_out.write("    cursor: pointer;\n")
        f_out.write("    user-select: none;\n")
        f_out.write("}\n")
        f_out.write(".section-header:hover {\n")
        f_out.write("    color: var(--accent-color);\n")
        f_out.write("}\n")
        f_out.write(".section-header .chevron::after {\n")
        f_out.write("    content: '\\F282';\n")
        f_out.write("    font-family: 'bootstrap-icons';\n")
        f_out.write("    transition: transform 0.2s ease;\n")
        f_out.write("    display: inline-block;\n")
        f_out.write("}\n")
        f_out.write(".section-header.collapsed .chevron::after {\n")
        f_out.write("    transform: rotate(-90deg);\n")
        f_out.write("}\n")
        f_out.write(".table {\n")
        f_out.write("    margin-bottom: 0;\n")
        f_out.write("    border-collapse: separate;\n")
        f_out.write("    border-spacing: 0;\n")
        f_out.write("}\n")
        f_out.write(".table th {\n")
        f_out.write("    background-color: #f8fafc;\n")
        f_out.write("    color: #475569;\n")
        f_out.write("    font-weight: 600;\n")
        f_out.write("    border-bottom: 2px solid var(--border-color);\n")
        f_out.write("    text-align: center;\n")
        f_out.write("    vertical-align: middle;\n")
        f_out.write("    padding: 0.75rem;\n")
        f_out.write("}\n")
        f_out.write(".table td {\n")
        f_out.write("    padding: 0.75rem;\n")
        f_out.write("    vertical-align: middle;\n")
        f_out.write("    text-align: center;\n")
        f_out.write("}\n")
        f_out.write(".table-hover tbody tr:hover {\n")
        f_out.write("    background-color: #f1f5f9;\n")
        f_out.write("    transition: background-color 0.15s ease;\n")
        f_out.write("}\n")
        f_out.write(".table-bordered {\n")
        f_out.write("    border: 1px solid var(--border-color);\n")
        f_out.write("    border-radius: 8px;\n")
        f_out.write("    overflow: hidden;\n")
        f_out.write("}\n")
        f_out.write(".bg-fr {\n")
        f_out.write("    color: var(--color-fr) !important;\n")
        f_out.write("    background-color: var(--bg-fr) !important;\n")
        f_out.write("    box-shadow: none !important;\n")
        f_out.write("    font-family: SFMono-Regular, Menlo, Monaco, Consolas, \"Liberation Mono\", \"Courier New\", monospace;\n")
        f_out.write("    font-size: 0.85rem;\n")
        f_out.write("    word-break: break-all;\n")
        f_out.write("}\n")
        f_out.write(".bg-cdr1 {\n")
        f_out.write("    color: var(--color-cdr1) !important;\n")
        f_out.write("    background-color: var(--bg-cdr1) !important;\n")
        f_out.write("    box-shadow: none !important;\n")
        f_out.write("    font-family: SFMono-Regular, Menlo, Monaco, Consolas, \"Liberation Mono\", \"Courier New\", monospace;\n")
        f_out.write("    font-size: 0.85rem;\n")
        f_out.write("    font-weight: 600;\n")
        f_out.write("    word-break: break-all;\n")
        f_out.write("}\n")
        f_out.write(".bg-cdr2 {\n")
        f_out.write("    color: var(--color-cdr2) !important;\n")
        f_out.write("    background-color: var(--bg-cdr2) !important;\n")
        f_out.write("    box-shadow: none !important;\n")
        f_out.write("    font-family: SFMono-Regular, Menlo, Monaco, Consolas, \"Liberation Mono\", \"Courier New\", monospace;\n")
        f_out.write("    font-size: 0.85rem;\n")
        f_out.write("    font-weight: 600;\n")
        f_out.write("    word-break: break-all;\n")
        f_out.write("}\n")
        f_out.write(".bg-cdr3 {\n")
        f_out.write("    color: var(--color-cdr3) !important;\n")
        f_out.write("    background-color: var(--bg-cdr3) !important;\n")
        f_out.write("    box-shadow: none !important;\n")
        f_out.write("    font-family: SFMono-Regular, Menlo, Monaco, Consolas, \"Liberation Mono\", \"Courier New\", monospace;\n")
        f_out.write("    font-size: 0.85rem;\n")
        f_out.write("    font-weight: 600;\n")
        f_out.write("    word-break: break-all;\n")
        f_out.write("}\n")
        f_out.write("#BFR1, #BFR2, #BFR3, #BFR4 {\n")
        f_out.write("    background-color: #64748b !important;\n")
        f_out.write("    color: #ffffff !important;\n")
        f_out.write("}\n")
        f_out.write("#BCDR1 {\n")
        f_out.write("    background-color: #0284c7 !important;\n")
        f_out.write("    color: #ffffff !important;\n")
        f_out.write("}\n")
        f_out.write("#BCDR2 {\n")
        f_out.write("    background-color: #16a34a !important;\n")
        f_out.write("    color: #ffffff !important;\n")
        f_out.write("}\n")
        f_out.write("#BCDR3 {\n")
        f_out.write("    background-color: #9333ea !important;\n")
        f_out.write("    color: #ffffff !important;\n")
        f_out.write("}\n")
        f_out.write("#idread, #idcolumn {\n")
        f_out.write("    background-color: #475569 !important;\n")
        f_out.write("    color: #ffffff !important;\n")
        f_out.write("}\n")
        f_out.write(".germlinerow {\n")
        f_out.write("    background-color: #475569 !important;\n")
        f_out.write("    color: #ffffff !important;\n")
        f_out.write("}\n")
        f_out.write(".popover {\n")
        f_out.write("    font-family: 'Inter', sans-serif;\n")
        f_out.write("    box-shadow: 0 10px 15px -3px rgba(0, 0, 0, 0.1);\n")
        f_out.write("    border: 1px solid var(--border-color);\n")
        f_out.write("}\n")
        f_out.write(".plot-card {\n")
        f_out.write("    background-color: #ffffff;\n")
        f_out.write("    border: 1px solid var(--border-color);\n")
        f_out.write("    border-radius: 8px;\n")
        f_out.write("    padding: 1rem;\n")
        f_out.write("    box-shadow: 0 1px 2px 0 rgba(0, 0, 0, 0.05);\n")
        f_out.write("    text-align: center;\n")
        f_out.write("}\n")
        f_out.write(".plot-card img {\n")
        f_out.write("    max-width: 100%;\n")
        f_out.write("    height: auto;\n")
        f_out.write("    border-radius: 4px;\n")
        f_out.write("}\n")
        f_out.write("footer {\n")
        f_out.write("    margin-top: 4rem;\n")
        f_out.write("    border-top: 1px solid var(--border-color);\n")
        f_out.write("    padding-top: 2rem;\n")
        f_out.write("    padding-bottom: 2rem;\n")
        f_out.write("}\n")
        f_out.write("#tablefooter {\n")
        f_out.write("    width: 100%;\n")
        f_out.write("}\n")
        f_out.write("#bioinfo img {\n")
        f_out.write("    max-height: 50px;\n")
        f_out.write("    object-fit: contain;\n")
        f_out.write("}\n")
        f_out.write("#unb img {\n")
        f_out.write("    max-height: 50px;\n")
        f_out.write("    object-fit: contain;\n")
        f_out.write("    float: right;\n")
        f_out.write("}\n")
        f_out.write("#address {\n")
        f_out.write("    text-align: left;\n")
        f_out.write("    font-size: 0.85rem;\n")
        f_out.write("    color: #64748b;\n")
        f_out.write("    line-height: 1.5;\n")
        f_out.write("}\n")
        f_out.write(".no-clones-msg {\n")
        f_out.write("    color: #ef4444;\n")
        f_out.write("    font-weight: 500;\n")
        f_out.write("    text-align: center;\n")
        f_out.write("    margin: 2rem 0;\n")
        f_out.write("}\n")
        f_out.write("</style>\n")
        f_out.write("</head>\n")
        f_out.write("<body>\n")
        f_out.write("<div class=\"container-fluid\">\n")
        f_out.write("<header>\n")
        f_out.write("<div id=\"application\">\n")
        f_out.write("<img src=\"Report/app5.png\" alt=\"ATTILA Logo\">\n")
        f_out.write("</div>\n")
        f_out.write("</header>\n")
        f_out.write("<ul class=\"nav nav-tabs\" id=\"libraryTabs\" role=\"tablist\">\n")
        f_out.write("<li class=\"nav-item\" role=\"presentation\">\n")
        f_out.write("<button class=\"nav-link active\" id=\"vh-tab\" data-bs-toggle=\"tab\" data-bs-target=\"#vh\" type=\"button\" role=\"tab\" aria-controls=\"vh\" aria-selected=\"true\">VH Library</button>\n")
        f_out.write("</li>\n")
        f_out.write("<li class=\"nav-item\" role=\"presentation\">\n")
        f_out.write("<button class=\"nav-link\" id=\"vl-tab\" data-bs-toggle=\"tab\" data-bs-target=\"#vl\" type=\"button\" role=\"tab\" aria-controls=\"vl\" aria-selected=\"false\">VL Library</button>\n")
        f_out.write("</li>\n")
        f_out.write("<li class=\"nav-item\" role=\"presentation\">\n")
        f_out.write("<button class=\"nav-link\" id=\"logs-tab\" data-bs-toggle=\"tab\" data-bs-target=\"#logs\" type=\"button\" role=\"tab\" aria-controls=\"logs\" aria-selected=\"false\">Log Files</button>\n")
        f_out.write("</li>\n")
        f_out.write("</ul>\n")
        f_out.write("<div class=\"tab-content\" id=\"libraryTabsContent\">\n")
        
        for lib in [0, 1]:
            if lib == 0:
                f_out.write("<div id=\"vh\" class=\"tab-pane fade show active\" role=\"tabpanel\" aria-labelledby=\"vh-tab\">\n")
                start = 0
                end = countvh - 1
                path1 = dir1
                path2 = dir2
                input_f = in1
                aux_n = 'vh'
            else:
                f_out.write("<div id=\"vl\" class=\"tab-pane fade\" role=\"tabpanel\" aria-labelledby=\"vl-tab\">\n")
                start = countvh
                end = countvl - 1
                path1 = dir3
                path2 = dir4
                input_f = in2
                aux_n = 'vl'
                
            has_candidates = os.path.exists(input_f) and count_fasta_records(input_f) > 0
            
            # Reads Information section
            loss_id = "lossreads" if lib == 0 else "lossreads1"
            f_out.write("<div class=\"card-section\">\n")
            f_out.write(f"<div class=\"section-header\" data-bs-toggle=\"collapse\" data-bs-target=\"#{loss_id}\" aria-expanded=\"true\" aria-controls=\"{loss_id}\">\n")
            f_out.write("<span>Reads Information</span>\n")
            f_out.write("<span class=\"chevron\"></span>\n")
            f_out.write("</div>\n")
            f_out.write(f"<div id=\"{loss_id}\" class=\"collapse show\">\n")
            f_out.write("<div class=\"table-responsive mb-4\">\n")
            f_out.write("<table class=\"table table-bordered table-striped table-hover align-middle\">\n")
            f_out.write("<thead>\n")
            f_out.write("<tr class=\"germlinerow text-white\">\n")
            f_out.write("<th>Dataset</th>\n")
            f_out.write("<th>Loss of Reads (Initial Library)</th>\n")
            f_out.write("<th>Loss of Reads (Final Library)</th>\n")
            f_out.write("</tr>\n")
            f_out.write("</thead>\n")
            f_out.write("<tbody>\n")
            f_out.write("<tr>\n")
            if lib == 0:
                f_out.write("<td>VH</td>\n")
                f_out.write(f"<td>{perdavhr0}%</td>\n")
                f_out.write(f"<td>{perdavhrn}%</td>\n")
            else:
                f_out.write("<td>VL</td>\n")
                f_out.write(f"<td>{perdavlr0}%</td>\n")
                f_out.write(f"<td>{perdavlrn}%</td>\n")
            f_out.write("</tr>\n")
            f_out.write("</tbody>\n")
            f_out.write("</table>\n")
            f_out.write("</div>\n")
            
            f_out.write("<div class=\"row g-4\">\n")
            f_out.write("<div class=\"col-md-6\">\n")
            f_out.write("<div class=\"plot-card\">\n")
            f_out.write("<h6 class=\"text-muted mb-3\">Reads Loss & Quality Statistics</h6>\n")
            f_out.write(f"<img src=\"{path1}\" alt=\"Plot 1\">\n")
            f_out.write("</div>\n")
            f_out.write("</div>\n")
            f_out.write("<div class=\"col-md-6\">\n")
            f_out.write("<div class=\"plot-card\">\n")
            f_out.write("<h6 class=\"text-muted mb-3\">Sequence Duplication Level</h6>\n")
            f_out.write(f"<img src=\"{path2}\" alt=\"Plot 2\">\n")
            f_out.write("</div>\n")
            f_out.write("</div>\n")
            f_out.write("</div>\n")
            f_out.write("</div>\n")
            f_out.write("</div>\n")
            
            if has_candidates:
                candidates_id = "candidates" if lib == 0 else "candidates1"
                alfa_val = alfavh if lib == 0 else alfavl
                f_out.write("<div class=\"card-section\">\n")
                f_out.write(f"<div class=\"section-header\" data-bs-toggle=\"collapse\" data-bs-target=\"#{candidates_id}\" aria-expanded=\"true\" aria-controls=\"{candidates_id}\">\n")
                f_out.write("<span>Candidate Clones</span>\n")
                f_out.write("<span class=\"chevron\"></span>\n")
                f_out.write("</div>\n")
                f_out.write(f"<div id=\"{candidates_id}\" class=\"collapse show\">\n")
                f_out.write("<div class=\"table-responsive\">\n")
                f_out.write("<table class=\"table table-bordered table-hover align-middle\">\n")
                f_out.write("<thead>\n")
                f_out.write("<tr class=\"germlinerow text-white\">\n")
                f_out.write("<th style=\"width: 12%\">Sequence ID</th>\n")
                f_out.write("<th style=\"width: 13%\">Fold Change</th>\n")
                f_out.write(f"<th style=\"width: 15%\">p-value (&alpha;<sub>Bonf</sub>={alfa_val})</th>\n")
                f_out.write("<th style=\"width: 15%\">InfIC (IC=95%)</th>\n")
                f_out.write("<th style=\"width: 15%\">SupIC (IC=95%)</th>\n")
                f_out.write("<th style=\"width: 18%\">Germline</th>\n")
                f_out.write("<th style=\"width: 12%\">Identity (%)</th>\n")
                f_out.write("</tr>\n")
                f_out.write("</thead>\n")
                f_out.write("<tbody>\n")
                
                counter = 0
                sorted_seqs = sorted(table1.items(), key=lambda x: x[1]['count'])
                for seq_id, info in sorted_seqs:
                    if start <= info['count'] <= end:
                        counter += 1
                        position = f"{counter}\xB0"
                        stats = stat_dict.get(seq_id, {'pvalue': '', 'infIC': '', 'supIC': ''})
                        f_out.write("<tr>\n")
                        f_out.write(f"<td><a class=\"btn btn-sm btn-outline-secondary py-0 px-2\" href=\"#\" data-bs-toggle=\"popover\" data-bs-trigger=\"hover\" data-bs-placement=\"top\" data-bs-content=\"{seq_id}\">{position}</a></td>\n")
                        f_out.write(f"<td>{info['fc']}</td>\n")
                        f_out.write(f"<td>{stats['pvalue']}</td>\n")
                        f_out.write(f"<td>{stats['infIC']}</td>\n")
                        f_out.write(f"<td>{stats['supIC']}</td>\n")
                        f_out.write(f"<td>{info['germline']}</td>\n")
                        f_out.write(f"<td>{info['identity']}</td>\n")
                        f_out.write("</tr>\n")
                f_out.write("</tbody>\n")
                f_out.write("</table>\n")
                f_out.write("</div>\n")
                f_out.write("</div>\n")
                f_out.write("</div>\n")
                
                # Regions of V Domain
                dom_id = "variabledom" if lib == 0 else "variabledom1"
                f_out.write("<div class=\"card-section\">\n")
                f_out.write(f"<div class=\"section-header\" data-bs-toggle=\"collapse\" data-bs-target=\"#{dom_id}\" aria-expanded=\"true\" aria-controls=\"{dom_id}\">\n")
                f_out.write("<span>Regions of Variable Domain of Candidate Clones</span>\n")
                f_out.write("<span class=\"chevron\"></span>\n")
                f_out.write("</div>\n")
                f_out.write(f"<div id=\"{dom_id}\" class=\"collapse show\">\n")
                f_out.write("<div class=\"table-responsive\">\n")
                f_out.write("<table class=\"table table-bordered align-middle\">\n")
                f_out.write("<thead>\n")
                f_out.write("<tr>\n")
                f_out.write("<th id=\"idread\" style=\"width: 8%\">Sequence ID</th>\n")
                f_out.write("<th id=\"BFR1\" style=\"width: 16%\">FR1</th>\n")
                f_out.write("<th id=\"BCDR1\" style=\"width: 10%\">CDR1</th>\n")
                f_out.write("<th id=\"BFR2\" style=\"width: 10%\">FR2</th>\n")
                f_out.write("<th id=\"BCDR2\" style=\"width: 10%\">CDR2</th>\n")
                f_out.write("<th id=\"BFR3\" style=\"width: 20%\">FR3</th>\n")
                f_out.write("<th id=\"BCDR3\" style=\"width: 12%\">CDR3</th>\n")
                f_out.write("<th id=\"BFR4\" style=\"width: 14%\">FR4</th>\n")
                f_out.write("</tr>\n")
                f_out.write("</thead>\n")
                f_out.write("<tbody>\n")
                
                counter = 0
                sorted_seqs2 = sorted(table2.items(), key=lambda x: x[1]['count'])
                for seq_id, info in sorted_seqs2:
                    if start <= info['count'] <= end:
                        counter += 1
                        position = f"{counter}\xB0"
                        f_out.write("<tr>\n")
                        f_out.write(f"<td id=\"idcolumn\"><a class=\"btn btn-sm btn-light py-0 px-2\" href=\"#\" data-bs-toggle=\"popover\" data-bs-trigger=\"hover\" data-bs-placement=\"top\" data-bs-content=\"{seq_id}\">{position}</a></td>\n")
                        f_out.write(f"<td class=\"bg-fr\">{info['fr1']}</td>\n")
                        f_out.write(f"<td class=\"bg-cdr1\">{info['cdr1']}</td>\n")
                        f_out.write(f"<td class=\"bg-fr\">{info['fr2']}</td>\n")
                        f_out.write(f"<td class=\"bg-cdr2\">{info['cdr2']}</td>\n")
                        f_out.write(f"<td class=\"bg-fr\">{info['fr3']}</td>\n")
                        f_out.write(f"<td class=\"bg-cdr3\">{info['cdr3']}</td>\n")
                        f_out.write(f"<td class=\"bg-fr\">{info['fr4']}</td>\n")
                        f_out.write("</tr>\n")
                f_out.write("</tbody>\n")
                f_out.write("</table>\n")
                f_out.write("</div>\n")
                f_out.write("</div>\n")
                f_out.write("</div>\n")
            else:
                f_out.write("<div class=\"card-section text-center\">\n")
                f_out.write("<p class=\"no-clones-msg\">No candidate clones found for this library.</p>\n")
                f_out.write("</div>\n")
                
            f_out.write("</div>\n")

        # Logs tab pane
        f_out.write("<div id=\"logs\" class=\"tab-pane fade\" role=\"tabpanel\" aria-labelledby=\"logs-tab\">\n")
        f_out.write("<div class=\"card-section\">\n")
        if log_entries:
            for rel, label in log_entries:
                f_out.write(f"<a href=\"{rel}\" target=\"_blank\" rel=\"noopener\" class=\"d-flex align-items-center gap-3 py-2 border-bottom text-decoration-none text-dark\">\n")
                f_out.write(f"<i class=\"bi bi-file-earmark-text text-secondary fs-5\"></i>\n")
                f_out.write(f"<span>{label}</span>\n")
                f_out.write(f"<i class=\"bi bi-box-arrow-up-right ms-auto text-muted\" style=\"font-size:0.85rem\"></i>\n")
                f_out.write(f"</a>\n")
        else:
            f_out.write("<p class=\"text-muted\">No log files found.</p>\n")
        f_out.write("</div>\n")
        f_out.write("</div>\n")

        f_out.write("</div>\n")
        f_out.write("</div>\n")

        # Footer
        f_out.write("<footer>\n")
        f_out.write("<div class=\"container-fluid\">\n")
        f_out.write("<table id=\"tablefooter\">\n")
        f_out.write("<tr>\n")
        f_out.write("<td id=\"bioinfo\" style=\"width: 10%\">\n")
        f_out.write("<img src=\"Report/bioinfo.png\" alt=\"Bioinfo\">\n")
        f_out.write("</td>\n")
        f_out.write("<td id=\"address\" style=\"width: 60%\">\n")
        f_out.write("<strong>Universidade de Bras&iacute;lia</strong><br/>\n")
        f_out.write("Instituto de Ci&ecirc;ncias Biol&oacute;gicas<br/>\n")
        f_out.write("Laborat&oacute;rio de Bioinform&aacute;tica, Bloco K<br/>\n")
        f_out.write("CEP 70.910-900 &middot; Bras&iacute;lia, DF, Brasil\n")
        f_out.write("</td>\n")
        f_out.write("<td id=\"unb\" style=\"width: 30%\">\n")
        f_out.write("<img src=\"Report/unb_logo1.png\" alt=\"UnB Logo\">\n")
        f_out.write("</td>\n")
        f_out.write("</tr>\n")
        f_out.write("</table>\n")
        f_out.write("</div>\n")
        f_out.write("</footer>\n")
        
        # Script initialization for Popovers and smooth toggles
        f_out.write("<script>\n")
        f_out.write("document.addEventListener('DOMContentLoaded', function () {\n")
        f_out.write("    // Initialize all Bootstrap popovers\n")
        f_out.write("    var popoverTriggerList = [].slice.call(document.querySelectorAll('[data-bs-toggle=\"popover\"]'))\n")
        f_out.write("    var popoverList = popoverTriggerList.map(function (popoverTriggerEl) {\n")
        f_out.write("        return new bootstrap.Popover(popoverTriggerEl)\n")
        f_out.write("    });\n")
        f_out.write("\n")
        f_out.write("    // Add toggle support to rotate collapse chevrons\n")
        f_out.write("    var collapseHeaders = document.querySelectorAll('.section-header');\n")
        f_out.write("    collapseHeaders.forEach(function (header) {\n")
        f_out.write("        var targetId = header.getAttribute('data-bs-target');\n")
        f_out.write("        var targetEl = document.querySelector(targetId);\n")
        f_out.write("        if (targetEl) {\n")
        f_out.write("            targetEl.addEventListener('show.bs.collapse', function () {\n")
        f_out.write("                header.classList.remove('collapsed');\n")
        f_out.write("            });\n")
        f_out.write("            targetEl.addEventListener('hide.bs.collapse', function () {\n")
        f_out.write("                header.classList.add('collapsed');\n")
        f_out.write("            });\n")
        f_out.write("        }\n")
        f_out.write("    });\n")
        f_out.write("});\n")
        f_out.write("</script>\n")
        f_out.write("</body>\n")
        f_out.write("</html>\n")
        
    print(f"index.html generated successfully at: {out}")

def count_fasta_records(fasta_path):
    if not os.path.exists(fasta_path):
        return 0
    cnt = 0
    with open(fasta_path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if line.startswith('>'):
                cnt += 1
    return cnt

if __name__ == '__main__':
    main()
