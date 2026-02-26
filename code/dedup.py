from pathlib import Path
import subprocess
from concurrent.futures import ProcessPoolExecutor

def run_dedup(sf):
    base_name = sf.stem
    output_path = Path("/home/emadeldin/HTRIBE_ncRNA/m38_20Feb2026/dedupped_reads")
    out_bam = output_path / f"{base_name}_dedup.bam"
    log_file = output_path / f"{base_name}_dedup.log"
    
    command = [
        "umi_tools", "dedup",
        f"--stdin={sf}",
        f"--log={log_file}",
        "--stdout", str(out_bam)
    ]
    
    print(f"Starting: {sf.name}")
    subprocess.run(command, check=True)
    return f"Finished: {sf.name}"

if __name__ == "__main__":
    input_path = Path("/home/emadeldin/HTRIBE_ncRNA/notrim_mm10/aligned_reads")
    output_path = Path("/home/emadeldin/HTRIBE_ncRNA/m38_20Feb2026/dedupped_reads")
    output_path.mkdir(parents=True, exist_ok=True)

    sam_files = sorted(input_path.glob("*.sam"))

    with ProcessPoolExecutor(max_workers=9) as executor:
        results = list(executor.map(run_dedup, sam_files))

    for result in results:
        print(result)
