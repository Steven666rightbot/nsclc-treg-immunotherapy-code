path = "D:/research/cucumber/fig1/xgboost_shap_analysis.py"
with open(path, "r", encoding="utf-8") as f:
    content = f.read()

# Add SHAP color override after "import shap"
old = "import shap\nfrom matplotlib.colors import LinearSegmentedColormap"
new = """import shap
# Override SHAP default red-blue with viridis (colorblind-safe)
import matplotlib.cm as cm
shap.plots.colors.red_blue = cm.get_cmap('viridis')
shap.plots.colors.red_blue_no_bounds = cm.get_cmap('viridis')
shap.plots.colors.red_blue_circle = cm.get_cmap('viridis')
from matplotlib.colors import LinearSegmentedColormap"""

if old in content:
    content = content.replace(old, new)
    print("Added SHAP color override")
else:
    print("Pattern not found")

with open(path, "w", encoding="utf-8") as f:
    f.write(content)
print("Done")
