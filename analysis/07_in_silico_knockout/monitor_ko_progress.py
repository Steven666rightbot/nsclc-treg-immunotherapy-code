#!/usr/bin/env python3
"""Monitor scTenifoldKnk KO progress and generate visualization."""

import time
import re
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

LOG_PATH = Path(r"C:\Users\qpmzw\.kimi\sessions\7230328cb54a1028d337e03b8acf3189\7224d5fb-dfb0-4504-acb8-c70a1af65246\tasks\bash-0bklj40p\output.log")
OUTPUT_PATH = Path(r"D:\Research\tomato\figures\ko_progress.png")
OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)

RECEPTORS = ["CD44", "ITGB1", "ITGA4", "DDR1", "DDR2"]
GROUPS = ["weak", "strong"]
TOTAL_KO = len(RECEPTORS) * len(GROUPS)  # 10

# Expected order: weak CD44, weak ITGB1, ..., weak DDR2, strong CD44, ...
KO_ORDER = [(g, r) for g in GROUPS for r in RECEPTORS]


def parse_log(log_path):
    """Parse log file to extract progress."""
    if not log_path.exists():
        return 0, None, []
    
    text = log_path.read_text(encoding='utf-8', errors='ignore')
    
    # Count completed KOs
    completed_matches = list(re.finditer(r'Completed\.\s+(\d+)\s+perturbed genes identified\.', text))
    n_completed = len(completed_matches)
    
    # Find last started KO
    current = None
    for m in re.finditer(r'========== KO: (\w+) \((\w+)\) ==========', text):
        current = (m.group(2), m.group(1))  # (group, receptor)
    
    # Check which specific KOs are done
    done_kos = []
    # Each "Completed" corresponds to a KO in order
    for i in range(min(n_completed, len(KO_ORDER))):
        done_kos.append(KO_ORDER[i])
    
    return n_completed, current, done_kos


def draw_progress(n_completed, current, done_kos, output_path):
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 3)
    ax.axis('off')
    
    # Title
    fig.suptitle('scTenifoldKnk Receptor KO Progress', fontsize=16, fontweight='bold', y=0.95)
    
    # Progress bar
    bar_y = 2.2
    bar_height = 0.4
    ax.add_patch(plt.Rectangle((0.5, bar_y), 9, bar_height, fill=False, edgecolor='gray', linewidth=2))
    if n_completed > 0:
        fill_width = 9 * (n_completed / TOTAL_KO)
        ax.add_patch(plt.Rectangle((0.5, bar_y), fill_width, bar_height, facecolor='#4CAF50', edgecolor='none'))
    
    pct = n_completed / TOTAL_KO * 100
    ax.text(5, bar_y + bar_height/2, f'{n_completed}/{TOTAL_KO} ({pct:.0f}%)', 
            ha='center', va='center', fontsize=14, fontweight='bold')
    
    # Grid of KOs
    grid_y = 0.8
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
            color = '#4CAF50'
            text_color = 'white'
            label = '✓'
        elif is_current:
            color = '#2196F3'
            text_color = 'white'
            label = '⋯'
        else:
            color = '#E0E0E0'
            text_color = '#666666'
            label = '○'
        
        ax.add_patch(plt.Rectangle((x, y), cell_w - 0.1, cell_h, facecolor=color, edgecolor='black', linewidth=1))
        ax.text(x + (cell_w - 0.1)/2, y + cell_h/2, f'{label}\n{receptor}\n({group})', 
                ha='center', va='center', fontsize=9, color=text_color, fontweight='bold')
    
    # Status text
    if current and current not in done_kos:
        status = f"Running: {current[1]} ({current[0]})"
    elif n_completed >= TOTAL_KO:
        status = "All done!"
    else:
        status = "Waiting..."
    
    ax.text(5, 0.1, status, ha='center', va='center', fontsize=12, style='italic', color='#333333')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Progress saved: {output_path}")


def main():
    print("Starting KO progress monitor...")
    while True:
        try:
            n_completed, current, done_kos = parse_log(LOG_PATH)
            draw_progress(n_completed, current, done_kos, OUTPUT_PATH)
            print(f"[{time.strftime('%H:%M:%S')}] Completed: {n_completed}/{TOTAL_KO}, Current: {current}")
            
            if n_completed >= TOTAL_KO:
                print("All KOs completed!")
                break
        except Exception as e:
            print(f"Error: {e}")
        
        time.sleep(30)


if __name__ == '__main__':
    main()
