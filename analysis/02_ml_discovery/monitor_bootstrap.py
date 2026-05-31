"""
Monitor CellChat bootstrap progress and auto-extract stats from new rds files.
Run in background: python monitor_bootstrap.py
"""
import os, glob, pickle, time, json
from datetime import datetime

BOOTSTRAP_DIR = r"D:\Research\tomato\results\cellchat_ko_corrected\bootstrap"
LOG_FILE = os.path.join(BOOTSTRAP_DIR, "bootstrap_monitor_log.json")
STRONG_RDS = r"D:\Research\tomato\results\cellchat_ko_corrected\cellchat_strong.rds"

SEED_MAP = {i: 20250516 + i * 1000 for i in range(7, 27)}

def extract_fib_treg(rds_path):
    try:
        obj = pickle.load(open(rds_path, 'rb'))
        df = obj['net']['df']
        ft = df[(df['source'] == 'Fibroblast') & (df['target'] == 'Treg')]
        return float(ft['prob'].sum()), int(len(ft))
    except Exception as e:
        return None, str(e)

def get_strong_baseline():
    try:
        obj = pickle.load(open(STRONG_RDS, 'rb'))
        df = obj['net']['df']
        ft = df[(df['source'] == 'Fibroblast') & (df['target'] == 'Treg')]
        return float(ft['prob'].sum()), int(len(ft))
    except:
        return None, None

def main():
    print(f"[{datetime.now()}] Bootstrap monitor started")
    strong_prob, strong_n = get_strong_baseline()
    print(f"  Strong baseline: prob={strong_prob}, n={strong_n}")
    
    while True:
        rds_files = sorted(glob.glob(os.path.join(BOOTSTRAP_DIR, "cellchat_weak_bootstrap_*.rds")))
        results = []
        for f in rds_files:
            name = os.path.basename(f)
            # Extract iter number
            try:
                iter_num = int(name.replace("cellchat_weak_bootstrap_", "").replace(".rds", ""))
            except:
                continue
            mtime = datetime.fromtimestamp(os.path.getmtime(f)).isoformat()
            size_mb = os.path.getsize(f) / 1024 / 1024
            
            prob, n_int = extract_fib_treg(f)
            seed = SEED_MAP.get(iter_num, None)
            
            results.append({
                "iter": iter_num,
                "seed": seed,
                "file": name,
                "modified": mtime,
                "size_mb": round(size_mb, 1),
                "fib_treg_prob": prob,
                "fib_treg_n": n_int
            })
        
        # Sort by iter
        results.sort(key=lambda x: x['iter'])
        
        # Compute stats for iter >= 7 (new)
        new_results = [r for r in results if r['iter'] >= 7 and r['fib_treg_prob'] is not None]
        old_results = [r for r in results if r['iter'] <= 6 and r['fib_treg_prob'] is not None]
        
        all_valid = [r['fib_treg_prob'] for r in results if r['fib_treg_prob'] is not None]
        
        summary = {
            "timestamp": datetime.now().isoformat(),
            "total_iterations_found": len(results),
            "old_iterations_01_06": len(old_results),
            "new_iterations_07_26": len(new_results),
            "expected_new": 20,
            "remaining": 20 - len(new_results),
            "all_median_prob": round(float(sorted(all_valid)[len(all_valid)//2]), 10) if all_valid else None,
            "strong_baseline_prob": strong_prob
        }
        
        if all_valid and strong_prob:
            med = sorted(all_valid)[len(all_valid)//2]
            summary["weak_strong_ratio"] = round(med / strong_prob, 4)
        
        log_data = {"summary": summary, "details": results}
        with open(LOG_FILE, 'w') as f:
            json.dump(log_data, f, indent=2, default=str)
        
        print(f"[{datetime.now()}] {len(results)} total rds | {len(new_results)}/20 new done | {summary.get('weak_strong_ratio', 'N/A')}x")
        
        if len(new_results) >= 20:
            print(f"[{datetime.now()}] ALL 20 NEW ITERATIONS COMPLETE!")
            break
        
        time.sleep(300)  # Check every 5 minutes

if __name__ == "__main__":
    main()
