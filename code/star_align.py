import subprocess
from pathlib import Path

genome_index = Path("/home/emadeldin/HTRIBE_ncRNA/m38_20Feb2026/m38_index")
input_path   = Path("/home/emadeldin/HTRIBE_ncRNA/raw/mRNA_seq/processed/extract.trim")
output_path  = Path("/home/emadeldin/HTRIBE_ncRNA/m38_20Feb2026/aligned_reads")
output_path.mkdir(parents=True, exist_ok=True)

all_files = sorted(f.name for f in input_path.glob("*.fastq"))

read_pairs = {}
for fname in all_files:
    if fname.endswith("_R1.fastq"):
        r2 = fname.replace("_R1.fastq", "_R2.fastq")
        if r2 in all_files:
            read_pairs[fname] = r2

for forward_read, reverse_read in read_pairs.items():
    sample_id = forward_read.replace("_R1.fastq", "")
    out_prefix = output_path / f"{sample_id}_"

    command = [
        "STAR",
        "--genomeDir", str(genome_index),
	"--outSAMtype", "BAM", "SortedByCoordinate",
        "--runThreadN", "20", # this will error if you change it and chatgpt won't help you out
        "--readFilesIn",
        str(input_path / forward_read),
        str(input_path / reverse_read),
        "--outFileNamePrefix", str(out_prefix)
    ]

    print("Running:", " ".join(command))
    subprocess.run(command, check=True)
