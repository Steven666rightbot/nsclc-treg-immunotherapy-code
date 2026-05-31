#!/usr/bin/env python3
import time
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from monitor_ko_progress import parse_log, draw_progress

LOG_PATH = Path(r"C:\Users\qpmzw\.kimi\sessions\7230328cb54a1028d337e03b8acf3189\7224d5fb-dfb0-4504-acb8-c70a1af65246\tasks\bash-6ilt6y3x\output.log")
OUTPUT_PATH = Path(r"D:\Research\tomato\figures\ko_progress.png")

while True:
    try:
        n, c, d = parse_log(LOG_PATH)
        draw_progress(n, c, d, OUTPUT_PATH)
        ts = time.strftime("%H:%M:%S")
        print(f"[{ts}] {n}/10 done, current: {c}")
        if n >= 10:
            print("All KOs completed!")
            break
    except Exception as e:
        print(f"Error: {e}")
    time.sleep(30)
