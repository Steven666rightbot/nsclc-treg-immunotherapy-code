import pandas as pd
import os

output_dir = "D:/Research/tomato/results/scTenifoldKnk_ko_merged"

files = {
    "CD44": "ko_merged_receptor_CD44.csv",
    "ITGB1": "ko_merged_receptor_ITGB1.csv",
    "DDR1": "ko_merged_receptor_DDR1.csv",
    "ITGA4": "ko_merged_receptor_ITGA4.csv",
    "MEF2C": "ko_merged_TF_MEF2C.csv",
    "WWTR1": "ko_merged_TF_WWTR1.csv",
}

for name, fname in files.items():
    path = os.path.join(output_dir, fname)
    df = pd.read_csv(path)
    sig = df[df["p.adj"] < 0.05].sort_values("p.adj")
    print(f"\n=== {name} KO ===")
    print(f"Total perturbed: {len(df)}, Significant: {len(sig)}")
    if len(sig) > 0:
        print(sig[["gene", "Z", "p.value", "p.adj"]].head(15).to_string(index=False))
        
        # Count up vs down
        up = sum(sig["Z"] > 0)
        down = sum(sig["Z"] < 0)
        print(f"Direction: {up} up, {down} down")

# Check GO BP files
print("\n\n=== GO BP ITGB1 ===")
go = pd.read_csv(os.path.join(output_dir, "go_bp_receptor_ITGB1.csv"))
print(go[["Description", "pvalue", "p.adjust", "Count"]].head(10).to_string(index=False))

print("\n\n=== GO BP DDR1 ===")
go = pd.read_csv(os.path.join(output_dir, "go_bp_receptor_DDR1.csv"))
print(go[["Description", "pvalue", "p.adjust", "Count"]].head(10).to_string(index=False))

print("\n\n=== GO BP MEF2C ===")
go = pd.read_csv(os.path.join(output_dir, "go_bp_TF_MEF2C.csv"))
print(go[["Description", "pvalue", "p.adjust", "Count"]].head(10).to_string(index=False))
