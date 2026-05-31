import os

path = "D:/research/cucumber/fig1/xgboost_shap_analysis.py"
with open(path, "r", encoding="utf-8", errors="replace") as f:
    content = f.read()

content = content.replace("color=CMAP_SHAP", "cmap='viridis'")

with open(path, "w", encoding="utf-8") as f:
    f.write(content)

print("Fixed cmap parameter")
