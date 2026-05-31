#!/usr/bin/env python3
"""Standalone KO progress monitor."""

import time
import re
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

LOG_PATH = Path(r"C:\Users\qpmzw\.kimi\sessions\7230328cb54a1028d337e03b8acf3189\7224d5fb-dfb0-4504-acb8-c70a1af65246\tasks\bash-v32rcnnw\output.log")
OUTPUT_PATH = Path(r"D:\Research\tomato\figures\ko_progress.png")
OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)

RECEPTORS = ["CD44", "ITGB1", "ITGA4", "DDR1", "DDR2"]
GROUPS = ["weak", "strong"]
TOTAL_KO = 10
KO_ORDER = [(g, r) for g in GROUPS for r in RECEPTORS]


def parse_log(log_path):
    if not log_path.exists():
        return 0, None, []
    text = log_path.read_text(encoding='utf-8', errors='ignore')
    completed = len(re.findall(r'Completed\.\s+\d+\s+perturbed genes identified\.', text))
    current = None
    for m in re.finditer(r'========== KO: (\w+) \((\w+)\) ==========', text):
        current = (m.group(2), m.group(1))
    done_kos = [KO_ORDER[i] for i in range(min(completed, len(KO_ORDER)))]
    return completed, current, done_kos


def draw_progress(n_completed, current, done_kos, elapsed_sec):
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 3.2)
    ax.axis('off')
    fig.suptitle('scTenifoldKnk Receptor KO Progress', fontsize=16, fontweight='bold', y=0.98)

    # Progress bar
    bar_y = 2.4
    bar_height = 0.4
    ax.add_patch(plt.Rectangle((0.5, bar_y), 9, bar_height, fill=False, edgecolor='gray', linewidth=2))
    if n_completed > 0:
        fill_width = 9 * (n_completed / TOTAL_KO)
        ax.add_patch(plt.Rectangle((0.5, bar_y), fill_width, bar_height, facecolor='#4CAF50', edgecolor='none'))
    pct = n_completed / TOTAL_KO * 100
    ax.text(5, bar_y + bar_height/2, f'{n_completed}/{TOTAL_KO} ({pct:.0f}%)',
            ha='center', va='center', fontsize=14, fontweight='bold')

    # Grid
    grid_y = 1.0
    cell_w = 1.6
    cell_h = 0.5
    for i, (group, receptor) in enumerate(KO_ORDER):
        col = i % 5
        row = i // 5
        x = 0.5 + col * cell_w
        y = grid_y - row * (cell_h + 0.3)
        is_done = (group, receptor) in done_kos
        is_current = current == (group, receptor) and not is_done
        if is_done:
            color, text_color, label = '#4CAF50', 'white', '✓'
        elif is_current:
            color, text_color, label = '#2196F3', 'white', '⋯'
        else:
            color, text_color, label = '#E0E0E0', '#666666', '○'
        ax.add_patch(plt.Rectangle((x, y), cell_w - 0.1, cell_h, facecolor=color, edgecolor='black', linewidth=1))
        ax.text(x + (cell_w - 0.1)/2, y + cell_h/2, f'{label}\n{receptor}\n({group})',
                ha='center', va='center', fontsize=9, color=text_color, fontweight='bold')

    # Status text with elapsed time
    if current and current not in done_kos:
        mins = elapsed_sec // 60
        status = f"Running: {current[1]} ({current[0]}) — {mins}m elapsed"
    elif n_completed >= TOTAL_KO:
        status = "All done!"
    else:
        status = "Waiting..."
    ax.text(5, 0.15, status, ha='center', va='center', fontsize=12, style='italic', color='#333333')

    plt.tight_layout()
    plt.savefig(OUTPUT_PATH, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Progress saved: {OUTPUT_PATH}")


start_time = time.time()
while True:
    try:
        n, c, d = parse_log(LOG_PATH)
        elapsed = time.time() - start_time
        draw_progress(n, c, d, elapsed)
        ts = time.strftime("%H:%M:%S")
        print(f"[{ts}] {n}/10 done, current: {c}, elapsed: {elapsed:.0f}s")
        if n >= 10:
            print("All KOs completed!")
            break
    except Exception as e:
        print(f"Error: {e}")
    time.sleep(30)
