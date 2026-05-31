path = "D:/research/cucumber/fig1/xgboost_shap_analysis.py"
with open(path, "r", encoding="utf-8") as f:
    lines = f.readlines()
for i, line in enumerate(lines):
    if "color=CMAP_GWP" in line:
        lines[i] = line.replace("color=CMAP_GWP", "cmap='viridis'")
        print(f"Fixed line {i+1}: {lines[i].strip()}")
with open(path, "w", encoding="utf-8") as f:
    f.writelines(lines)
print("Done")
