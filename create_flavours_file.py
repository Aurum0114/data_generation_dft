import pprint
import json

def create_flavours_file(flavours_source_file, forbidden):

    with open(flavours_source_file, 'r') as f:
        func_and_basis = json.load(f)
    
    functionals = func_and_basis["functionals"]
    basissets = func_and_basis["basissets"]

    dft_flavours = []
    index = 1
    for f in functionals:
        for b in basissets:
            if (f, b) not in forbidden:
                t = {
                    "number": str(index),
                    "functional": f,
                    "basisset": b,
                    "num_molecules": 0
                }
                dft_flavours.append(t)
                index += 1      

    num_flavours = len(dft_flavours)
    print(f"Found {num_flavours} flavours")
    pprint.pprint(dft_flavours)

    with open('final_flavours.json', 'w') as f:
        json.dump(dft_flavours, f)

    return dft_flavours

forbidden = [('bmk', 'aug-cc-pVDZ'), ('b3-lyp', 'aug-cc-pVDZ'), ('b3-lyp', '6-311++G**'), ('pbe0', 'aug-cc-pVDZ'), ('pbe0', '6-311++G**'), ('tpssh', 'aug-cc-pVDZ'), ('tpssh', '6-311++G**'), ('m06-2x', 'aug-cc-pVDZ'), ('m06-2x', '6-311++G**')]
input = 'input_files/flavours.json'
create_flavours_file(input, forbidden)


