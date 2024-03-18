import argparse
import json
import os
import random
import shutil

import utils.xyz_utils as xyz

def create_placeholder_categories(flavour_file, molecule_xyz_file, output_temp_dir, num_molecules, flavour_num=None):
    """
    This function creates placeholder categories that will be used by the next function to create the actual categories.
    :flavour_file: Stores two lists: all functionals and basis sets. Example file: example_files/func_and_base.json
    :molecule_xyz_file: Filtered dataset to sample molecules from. Example file: example_files/inputs.xyz
    :output_temp_dir: The directory that stores the temporary placeholder categories
    :num_molecules: The number of molecules per datatype per category
    """
    if os.path.exists(output_temp_dir):
        shutil.rmtree(output_temp_dir)

    temp_task_dir = os.path.join(output_temp_dir, "tasks")
    os.makedirs(temp_task_dir)

    with open(flavour_file, 'r') as fp:
        func_and_basis = json.load(fp)

    functionals = func_and_basis["functionals"]
    basissets = func_and_basis["basissets"]

    dft_flavours = []
    index = 1
    for f in functionals:
        for b in basissets:
            t = {
                "number": str(index),
                "functional": f,
                "basisset": b,
                "num_molecules": num_molecules
            }
            dft_flavours.append(t)
            index += 1

    num_flavours = len(dft_flavours)
    print("Found number of flavours is: ", num_flavours)

    coords_all, elements_all = xyz.readXYZs(molecule_xyz_file)
    assert len(coords_all) == len(elements_all)

    if flavour_num==None:
        for single_flavour in dft_flavours:
            sample_molecules_for_flavour(single_flavour, dft_flavours, coords_all, elements_all, temp_task_dir)

    elif str(flavour_num).isnumeric() and 1 <= int(flavour_num) < len(dft_flavours):
        single_flavour = dft_flavours[flavour_num]
        sample_molecules_for_flavour(single_flavour, dft_flavours, coords_all, elements_all, temp_task_dir)
        
    else:
        raise ValueError(f"Flavour value should be an integer (index of the flavour to use), maximum value is {dft_flavours[-1]['number']}")
    

def sample_molecules_for_flavour(single_flavour, dft_flavours, coords_all, elements_all, temp_task_dir):
    num_all_mol = len(coords_all)
    sampled_indices = random.sample(range(num_all_mol), int(single_flavour["num_molecules"]))
    coords = [coords_all[i] for i in sampled_indices]
    elements = [elements_all[i] for i in sampled_indices]

    task_coord_filename = f"data_01_{single_flavour['functional']}###{single_flavour['basisset']}.xyz"
    task_dir_path = os.path.join(temp_task_dir,
                                f"flavour{dft_flavours['number']}_{single_flavour['functional']}###{single_flavour['basisset']}")

    os.mkdir(task_dir_path)
    print(f"Creating directory {task_dir_path}")
    xyz.exportXYZs(coords, elements, os.path.join(task_dir_path, task_coord_filename))

    with open(os.path.join(task_dir_path, "task_info.json"), 'w') as fp:
        json.dump(single_flavour, fp)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('flavour_file')
    parser.add_argument('molecule_xyz_file')
    parser.add_argument('output_temp_dir')
    parser.add_argument('num_molecules')
    parser.add_argument('--flavour_num', type=int, default=None)
    
    args = parser.parse_args()
    print("Creating placeholder categories ... ")
    create_placeholder_categories(args.flavour_file, args.molecule_xyz_file, args.output_temp_dir, int(args.num_molecules), args.flavour_num)
    print("Done")