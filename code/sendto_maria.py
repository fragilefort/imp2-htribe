import subprocess
import re
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

def sendto_maria(file_info):
    sf, exp, number = file_info
    table_name = "lastHope"
    command = [
        "/home/emadeldin/HTRIBE_ncRNA/m38_20Feb2026/code/load_table.sh",
        str(sf),  
        table_name,
        exp,
        str(number)
    ]

    print(f"Uploading {sf.name} as {exp} (ID: {number})")
    
    result = subprocess.run(command, check=True, capture_output=True, text=True)
    return f"Done: {sf.name}"

if __name__ == "__main__":
    input_path = Path("/home/emadeldin/HTRIBE_ncRNA/m38_20Feb2026/dedupped_reads")
    sam_files = sorted(input_path.glob("*.sam"))
    
    tasks = []

    for sf in sam_files:
        match = re.search(r"DS(\d+)", sf.name)
        
        if match:
            n = int(match.group(1))
            
            if 1 <= n <= 3:
                exp = "wt"
            elif 4 <= n <= 6:
                exp = "mcherry"
            elif 7 <= n <= 9:
                exp = "imp2"
            else:
                continue 
            
            tasks.append((sf, exp, n))

    with ProcessPoolExecutor(max_workers=9) as executor:
        results = list(executor.map(sendto_maria, tasks))

    for r in results:
        print(r)
