import os
log_path = '/hkfs/work/workspace/scratch/mt4181-dft/small_db/log'
first_job_name = 'job2346574'


ff = {
    "functionals": ["b-p", "bmk", "tpss", "b3-lyp", "bh-lyp", "pbe0", "tpss", "tpssh", "m06-l", "m06-2x", "b97-d", "b97-3c"],
    "basissets": ["def-SVP", "def2-SV(P)", "def2-TZVP", "def2-TZVPP", "cc-pVDZ", "aug-cc-pVDZ", "cc-pVTZ", "6-31G", "6-311G", "6-311++G**"]
}


dft_flavours = []
for f in ff['functionals']:
    for b in ff['basissets']:
        t = (f, b)
        dft_flavours.append(t)

report = []
for i in range(1, 121):
    file = f'job{2346573+i}_flv{i}.output'
    file_path = os.path.join(log_path, file)

    with open(file_path, 'r') as f:
        for line in f:
            words = line.split()
            if len(words) == 4 and words[1] == 'Wall-clock':
                report.append(f"flv{i}, {dft_flavours[i-1]}, {words[3]}")

print(report)


