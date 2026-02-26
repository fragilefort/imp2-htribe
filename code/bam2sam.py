from pathlib import Path
import subprocess
from concurrent.futures import ProcessPoolExecutor


def bam2sam(bam_file):
    base_name = bam_file.stem
    out_path = Path("/home/emadeldin/HTRIBE_ncRNA/m38_20Feb2026/dedupped_reads")
    out_sam = out_path / f"{base_name}_dedup.sam"

    command = [
        "samtools", "view",
        "-h", "-o", str(out_sam),
        bam_file
    ]

    print(f"Starting {bam_file.name}")
    subprocess.run(command, check = True)
    print(f"Finished {bam_file.name}")


if __name__ == "__main__":

    input_path = Path("/home/emadeldin/HTRIBE_ncRNA/m38_20Feb2026/dedupped_reads")
    bam_files = sorted(input_path.glob("*.bam"))
    with ProcessPoolExecutor(max_workers=9) as executor:
        results = list(executor.map(bam2sam, bam_files))

