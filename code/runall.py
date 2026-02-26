from pathlib import Path
import subprocess
from concurrent.futures import ProcessPoolExecutor

def run_bash_script(script_path):
    print(f"Starting: {script_path.name} ")
    
    try:
        result = subprocess.run(
            ["bash", str(script_path)], 
            check=True, 
            capture_output=True, 
            text=True
        )
        return f"GOOD: {script_path.name}\nOutput: {result.stdout[:100]}"
    except subprocess.CalledProcessError as e:
        return f"ERROR: {script_path.name}\nExit Code: {e.returncode}\nError: {e.stderr}"

if __name__ == "__main__":
    scripts_dir = Path("/home/emadeldin/HTRIBE_ncRNA/m38_20Feb2026/code/find_edit_sites")
    
    bash_scripts = sorted(scripts_dir.glob("*.sh"))

    print(f"Found {len(bash_scripts)} scripts to run.")

    with ProcessPoolExecutor(max_workers=9) as executor:
        results = list(executor.map(run_bash_script, bash_scripts))

    for result in results:
        print(result)
