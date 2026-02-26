import random
import sys
import subprocess
import os
from pathlib import Path
from itertools import islice

if len(sys.argv) < 4:
    print("Usage: python simulate.py <normal.fa> <edited.fa> <out_dir>")
    sys.exit(1)

normal_fasta = sys.argv[1]
edited_fasta = sys.argv[2]
out_dir = Path(sys.argv[3])
out_dir.mkdir(parents=True, exist_ok=True)

READ_LEN = 100
FOLD = 10000
N_TARGET = 100000 

def run_art():
    subprocess.run(['art_illumina', '-p', '-i', normal_fasta, '-f', str(FOLD), 
                    '-l', str(READ_LEN), '-m', '200', '-s', '10', '-o', str(out_dir / "POOL_NORM_")])
    subprocess.run(['art_illumina', '-p', '-i', edited_fasta, '-f', str(FOLD), 
                    '-l', str(READ_LEN), '-m', '200', '-s', '10', '-o', str(out_dir / "POOL_EDITED_")])

run_art()

def generate_sample(n_target, out_prefix, ratio=0.0):
    """
    ratio=0.0 means 100% normal.
    ratio=0.1 means 10% edited, 90% normal.
    """
    reservoir = [] 
    
    files = {
        'n1': open(out_dir / "POOL_NORM_1.fq"),
        'n2': open(out_dir / "POOL_NORM_2.fq"),
        'e1': open(out_dir / "POOL_EDITED_1.fq"),
        'e2': open(out_dir / "POOL_EDITED_2.fq")
    }

    try:
        count = 0
        while True:
            # Decide pool
            is_edited = random.random() < ratio
            r1_fh = files['e1'] if is_edited else files['n1']
            r2_fh = files['e2'] if is_edited else files['n2']

            block1 = list(islice(r1_fh, 4))
            block2 = list(islice(r2_fh, 4))

            if not block1 or not block2:
                break
            
            count += 1
            if len(reservoir) < n_target:
                reservoir.append((block1, block2))
            else:
                m = random.randint(0, count - 1)
                if m < n_target:
                    reservoir[m] = (block1, block2)
        
        with open(out_dir / f"{out_prefix}_R1.fq", 'w') as f1, \
             open(out_dir / f"{out_prefix}_R2.fq", 'w') as f2:
            for r1, r2 in reservoir:
                f1.writelines(r1)
                f2.writelines(r2)
                
    finally:
        for fh in files.values():
            fh.close()

for i in range(1, 10):
    sample_name = f"ctrl_DS{i}"
    current_ratio = 0.0 if i <= 6 else 0.3
    
    print(f"Generating {sample_name} (Ratio: {current_ratio})...")
    generate_sample(N_TARGET, sample_name, ratio=current_ratio)
