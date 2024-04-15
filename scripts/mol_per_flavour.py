import os
import json
import re

results_path = '/hkfs/work/workspace/scratch/mt4181-dft/final_db/results/tasks'

flavours = []
report = []
for item in os.listdir(results_path):
    item_path = os.path.join(results_path, item)
    if os.path.isdir(item_path):
        flavours.append((item_path, item))

for flavour_path, flavour_name in flavours:
    os.chdir(flavour_path)
    with open('task_info.json') as f:
        task_info = json.load(f)
    report.append((task_info['num_molecules'], flavour_name))


def extract_flv_num(item):
    name = item[1]
    match = re.search(r'flv(\d+)_', name)
    if match:
        return int(match.group(1))
    else:
        return float('inf')
    

sorted_report = sorted(report, key=extract_flv_num)
for i in sorted_report:
    print(i)
    
