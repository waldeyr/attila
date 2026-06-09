#!/usr/bin/env python3
"""
generate_synthetic_data.py - Generates realistic synthetic heavy (VH) and light (VL)
chain FASTQ and FASTA reads using human antibody framework regions to ensure
successful UCL Abnum alignment and ranking.
"""

import os

def create_synthetic_data(out_dir):
    os.makedirs(out_dir, exist_ok=True)
    
    # 1. Realistic human heavy chain framework (based on VH3-23 germline)
    # EVQLLESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAK (98 residues)
    fw_vh = (
        "GAGGTGCAGCTGTTGGAGTCTGGGGGAGGCTTGGTACAGCCTGGGGGGTCCCTGAGACTCTCCTGTG"
        "CAGCCTCTGGATTCACCTTTAGCAGCTATGCCATGAGCTGGGTCCGCCAGGCTCCAGGGAAGGGGCT"
        "GGAGTGGGTCAGCGCTATTAGTGGTAGTGGTGGTAGCACATACTACGCAGACTCCGTGAAGGGCCGG"
        "TTCACCATCTCCAGAGACAATTCCAAGAACACGCTGTATCTGCAAATGAACAGCCTGAGAGCCGAGG"
        "ACACGGCCGTGTATTACTGTGCGAAA"
    )
    # Standard human VH FR4: WGQGTLVTVSS (11 residues)
    fr4_vh = "TGGGGTCAAGGCACCCTGGTCACCGTCTCCTCA"
    
    # 2. Realistic human light chain framework (based on V1-13 lambda germline)
    # QSVLTQPPSVSGAPGQRVTISCTGSSSNIGNNYDVHWYQQLPGTAPKLLIYGNSNRPSGVPDRFSGSKSGTSASLAITGLQAEDEADYYC (90 residues)
    fw_vl = (
        "CAGTCTGTGCTGACTCAGCCACCCTCAGTGTCTGGGGCCCCAGGGCAGAGGGTCACCATCTCCTGCA"
        "CTGGCTCCAGCTCGAACATCGGGGCAGGTTATGATGTACACTGGTACCAGCAGCTTCCAGGAACAGC"
        "CCCCAAACTCCTCATCTATGGTAACAGCAATCGGCCCAGCGGGGTCCCTGACCGATTCTCTGGCTCC"
        "AAGTCTGGCACCTCAGCCTCCCTGGCCATCACTGGGCTCCAGGCTGAGGATGAGGCTGATTATTACT"
        "GC"
    )
    # Standard human VL FR4: FGQGTKLTVLG (11 residues)
    fr4_vl = "TTCGGCGGAGGGACCAAGCTGACCGTCCTAGGT"
    
    # CDR3 variants (DNA)
    # Variant A (Enriched): GGGGGGG (21 bp)
    cdr3_a = "GGTGGTGGTGGTGGTGGTGGT"
    # Variant B (Decreased): SSSSSSS (21 bp)
    cdr3_b = "TCTTCTTCTTCTTCTTCTTCT"
    # Variant C (Static): AAAAAAA (21 bp)
    cdr3_c = "GCTGCTGCTGCTGCTGCTGCT"
    # Variant D (Enriched): QQQQQQQ (21 bp)
    cdr3_d = "CAGCAGCAGCAGCAGCAGCAG"
    # Variant E (Enriched): YYYYYYY (21 bp)
    cdr3_e = "TATTATTATTATTATTATTAT"
    
    # Assembly functions
    def make_vh_seq(cdr3_dna):
        return fw_vh + cdr3_dna + fr4_vh
        
    def make_vl_seq(cdr3_dna):
        return fw_vl + cdr3_dna + fr4_vl

    # Quality strings (Q40, length matching reads)
    # VH length: 294 + 21 + 33 = 348 bp
    # VL length: 270 + 21 + 33 = 324 bp
    qual_vh = "I" * 348
    qual_vl = "I" * 324
    
    # Define datasets
    # Initial round: A=1, B=5, C=3, D=1, E=2  (A, D, E are rare — will be selected)
    # Final round:   A=10, B=1, C=3, D=8, E=6 (3 enriched clones: A, D, E)
    vh_r0 = (
        [("read_VH_A_1", cdr3_a)] * 1 +
        [("read_VH_B_" + str(i), cdr3_b) for i in range(1, 6)] +
        [("read_VH_C_" + str(i), cdr3_c) for i in range(1, 4)] +
        [("read_VH_D_1", cdr3_d)] * 1 +
        [("read_VH_E_" + str(i), cdr3_e) for i in range(1, 3)]
    )
    vh_rn = (
        [("read_VH_A_" + str(i), cdr3_a) for i in range(1, 11)] +
        [("read_VH_B_1", cdr3_b)] * 1 +
        [("read_VH_C_" + str(i), cdr3_c) for i in range(1, 4)] +
        [("read_VH_D_" + str(i), cdr3_d) for i in range(1, 9)] +
        [("read_VH_E_" + str(i), cdr3_e) for i in range(1, 7)]
    )

    vl_r0 = (
        [("read_VL_A_1", cdr3_a)] * 1 +
        [("read_VL_B_" + str(i), cdr3_b) for i in range(1, 6)] +
        [("read_VL_C_" + str(i), cdr3_c) for i in range(1, 4)] +
        [("read_VL_D_1", cdr3_d)] * 1 +
        [("read_VL_E_" + str(i), cdr3_e) for i in range(1, 3)]
    )
    vl_rn = (
        [("read_VL_A_" + str(i), cdr3_a) for i in range(1, 11)] +
        [("read_VL_B_1", cdr3_b)] * 1 +
        [("read_VL_C_" + str(i), cdr3_c) for i in range(1, 4)] +
        [("read_VL_D_" + str(i), cdr3_d) for i in range(1, 9)] +
        [("read_VL_E_" + str(i), cdr3_e) for i in range(1, 7)]
    )
    
    tasks = [
        ("vh_initial", vh_r0, make_vh_seq, qual_vh),
        ("vh_final", vh_rn, make_vh_seq, qual_vh),
        ("vl_initial", vl_r0, make_vl_seq, qual_vl),
        ("vl_final", vl_rn, make_vl_seq, qual_vl)
    ]
    
    for prefix, reads, make_seq_fn, qual in tasks:
        # A. Write FASTQ file
        fastq_path = os.path.join(out_dir, f"{prefix}.fastq")
        with open(fastq_path, "w", encoding="utf-8") as f_fq:
            for r_id, cdr3_dna in reads:
                seq = make_seq_fn(cdr3_dna)
                f_fq.write(f"@{r_id}\n{seq}\n+\n{qual}\n")
        print(f"Generated: {fastq_path}")
        
        # B. Write pre-filtered FASTA file
        fasta_path = os.path.join(out_dir, f"{prefix}filtered.fasta")
        with open(fasta_path, "w", encoding="utf-8") as f_fa:
            for r_id, cdr3_dna in reads:
                seq = make_seq_fn(cdr3_dna)
                f_fa.write(f">{r_id}\n{seq}\n")
        print(f"Generated: {fasta_path}")

if __name__ == "__main__":
    import sys
    out = os.path.dirname(os.path.abspath(__file__))
    if len(sys.argv) > 1:
        out = sys.argv[1]
    create_synthetic_data(out)
