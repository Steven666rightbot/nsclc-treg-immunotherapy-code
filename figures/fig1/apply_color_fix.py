import os

path = "D:/research/cucumber/fig1/xgboost_shap_analysis.py"
with open(path, "r", encoding="utf-8") as f:
    lines = f.readlines()

replacements = {
    "COLOR_HIGH = '#8bc98b'   # soft green\nCOLOR_LOW  = '#b08ebf'   # soft purple\nCMAP_GWP = LinearSegmentedColormap.from_list('gwp',\n    [COLOR_LOW, '#d4c2df', '#ffffff', '#c2e4c2', COLOR_HIGH])\nplt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans']\nplt.rcParams['axes.unicode_minus'] = False":
    "# ── Colorblind-safe colormap for SHAP (viridis) ──\n# viridis: perceptually uniform, readable by all color vision types\nimport matplotlib.cm as cm\nCMAP_SHAP = cm.get_cmap('viridis')\nplt.rcParams['font.sans-serif'] = ['Arial', 'Arial Unicode MS', 'DejaVu Sans']\nplt.rcParams['axes.unicode_minus'] = False",

    "        color=CMAP_GWP\n    )":
    "        cmap='viridis'\n    )",

    "colors = CMAP_GWP(np.linspace(0, 1, top_n))":
    "colors = CMAP_SHAP(np.linspace(0, 1, top_n))",

    "# Green-purple gradient for 5 folds\nfold_colors = ['#c2a5cf', '#d4c2df', '#e8e8e8', '#c2e4c2', '#a6dba0']":
    "# Colorblind-safe fold colors (Okabe-Ito inspired, distinguishable)\nfold_colors = ['#56B4E9', '#009E73', '#F0E442', '#E69F00', '#0072B2']",

    "ax.plot(mean_fpr, mean_tpr, color='#5a9e5a', lw=2.5,\n        label=f'Mean (AUC = {mean_auc:.3f} ± {std_auc:.3f})')\nax.fill_between(mean_fpr, np.maximum(mean_tpr - np.std(tprs, axis=0), 0),\n                np.minimum(mean_tpr + np.std(tprs, axis=0), 1), color='#c2e4c2', alpha=0.35)":
    "# Mean curve: project standard green (MPR/Responder)\nax.plot(mean_fpr, mean_tpr, color='#5DAE5D', lw=2.5,\n        label=f'Mean (AUC = {mean_auc:.3f} ± {std_auc:.3f})')\nax.fill_between(mean_fpr, np.maximum(mean_tpr - np.std(tprs, axis=0), 0),\n                np.minimum(mean_tpr + np.std(tprs, axis=0), 1), color='#5DAE5D', alpha=0.15)",
}

content = ''.join(lines)
for old, new in replacements.items():
    if old in content:
        content = content.replace(old, new)
        print(f"Replaced: {old[:50]}...")
    else:
        print(f"NOT FOUND: {old[:50]}...")

with open(path, "w", encoding="utf-8") as f:
    f.write(content)

print("Done")
