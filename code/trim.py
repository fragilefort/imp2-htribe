from pathlib import Path
import subprocess

input_path = Path("/home/emadeldin/HTRIBE_ncRNA/raw/mRNA_seq/processed/extract.trim")
output_path = Path("/home/emadeldin/HTRIBE_ncRNA/mm10_run/trimmed")
output_path.mkdir(exist_ok=True)

all_files = sorted([f.name for f in input_path.glob("*.fastq")])

read_pairs = {}
for filename in all_files:
    if "_R1" in filename:
        reverse_candidate = filename.replace("_R1", "_R2")
        if reverse_candidate in all_files:
            read_pairs[filename] = reverse_candidate

for forward_read, reverse_read in read_pairs.items():
    print(f"Processing: {forward_read} and {reverse_read}...")
    command = [
        "cutadapt",
        "--trim-n", "--match-read-wildcards", "-n", "4", "-u", "15",
        "-a", "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
        "-a", "CCAAGTCT", "-a", "ATCTCGTATGCCGTCTTCTGCTTG",
        "-a", "GGGGGGGGGG", "-a", "AAAAAAAAAA", "-A", "TTTTTTTTTT",
        "-e", "0.2", "-j", "0", "--nextseq-trim", "20", "-m", "40:40",
        "-o", str(output_path / forward_read),
        "-p", str(output_path / reverse_read),
        str(input_path / forward_read),
        str(input_path / reverse_read)
    ]

    subprocess.run(command, check=True)
