import pandas as pd
import argparse
import json
import os
import random
import shutil
import re

import utils.xyz_utils as xyz

def create_placeholder_categories(flavour_file, xyz_files_dir, output_temp_dir, num_molecules, flavour_idx=None):
    """
    This function creates placeholder categories that will be used by the next function to create the actual categories.
    :flavour_file: Stores two lists: all functionals and basis sets. Example file: example_files/func_and_base.json
    :xyz_files_dir: Filtered dataset to sample molecules from, contained the data in folders
    :output_temp_dir: The directory that stores the temporary placeholder categories
    :num_molecules: The number of molecules per datatype per category
    """
    #if os.path.exists(output_temp_dir):
    #    shutil.rmtree(output_temp_dir)

    temp_task_dir = os.path.join(output_temp_dir, "tasks")
    if not os.path.exists(temp_task_dir):
        os.makedirs(temp_task_dir)

    with open(flavour_file, 'r') as f:
        func_and_basis = json.load(f)

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
    print(f"Found {num_flavours} flavours")

    all_items = os.listdir(xyz_files_dir)
    datatypes_paths = [os.path.join(xyz_files_dir, item) for item in all_items 
                      if os.path.isdir(os.path.join(xyz_files_dir, item))]
    print(f"Found datatypes are {datatypes_paths}")

    if flavour_idx==None:
        for single_flavour in dft_flavours:
            export_molecules_for_flavour(single_flavour, datatypes_paths, temp_task_dir)

    elif str(flavour_idx).isnumeric() and 1 <= int(flavour_idx) <= len(dft_flavours):
        single_flavour = dft_flavours[flavour_idx-1]
        export_molecules_for_flavour(single_flavour, datatypes_paths, temp_task_dir)
        
    else:
        raise ValueError(f"Flavour value should be an integer (index of the flavour to use), maximum value is {dft_flavours[-1]['number']}")

def export_molecules_for_flavour(single_flavour, datatypes_paths, temp_task_dir):
    datapoints_to_sample = int(single_flavour["num_molecules"])
    working_dir = os.getcwd()
    coords_all = []
    elements_all = []
    charges_all = []

    for datatype_path in datatypes_paths:
        os.chdir(datatype_path)
        all_files = os.listdir()
        
        #change to se3 if necessary
        file_pattern_coords = r'.*se1_ID_0_coords.xyz$'
        file_pattern_info = r'.*se1_ID_0_info.csv$'

        pattern_coords = re.compile(file_pattern_coords)
        pattern_info = re.compile(file_pattern_info)

        coords_file = [file for file in all_files if pattern_coords.match(file)][0]
        print(f"Coords file for data type {datatype_path} is {coords_file}")
        coords_all_1type, elements_all_1type = xyz.readXYZs(coords_file)
        assert len(coords_all) == len(elements_all)
        
        info_file_path = [file for file in all_files if pattern_info.match(file)][0]
        info_df = pd.read_csv(info_file_path, usecols=['ID', 'names'])
        ids = info_df['ID'].tolist()
        system_names = info_df['names'].tolist()
        
        num_all_mol = len(coords_all_1type)
        sampled_indices = random.sample(range(num_all_mol), datapoints_to_sample)

        coords_sampled = [coords_all_1type[i] for i in sampled_indices]
        elements_sampled = [elements_all_1type[i] for i in sampled_indices]
        charges_sampled = [find_mol_charge(i, ids, system_names) for i in sampled_indices]

        coords_all.extend(coords_sampled)
        elements_all.extend(elements_sampled)
        charges_all.extend(charges_sampled)

        os.chdir(working_dir)

    assert len(coords_all) == len(elements_all) == len(charges_all)
    print(f"in total {len(coords_all)} molecules were collected")
    print(f"Data was collected for {len(coords_all)/datapoints_to_sample} data types")

    task_coord_filename = f"data_01_{single_flavour['functional']}###{single_flavour['basisset']}.xyz"
    task_dir_path = os.path.join(temp_task_dir,
                                f"FLV_{single_flavour['number']}_{single_flavour['functional']}###{single_flavour['basisset']}")

    if not os.path.exists(task_dir_path):
        os.mkdir(task_dir_path)
        print(f"Creating directory {task_dir_path}")

    xyz.exportXYZs(coords_all, elements_all, charges_all, os.path.join(task_dir_path, task_coord_filename))

    with open(os.path.join(task_dir_path, "task_info.json"), 'w') as fp:
        json.dump(single_flavour, fp)
        print("hell yeah it worked")
        

def find_mol_charge(i, ids, system_names):
    system = system_names[i]
    pos_chr = system.count('+')
    neg_chr = system.count('-')
    tot_chr = pos_chr - neg_chr   
    print(f"For system {system}, ID {ids[i]} = {i}, the charge is {tot_chr}")

    return tot_chr

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('flavour_file')
    parser.add_argument('molecule_xyz_file')
    parser.add_argument('output_temp_dir')
    parser.add_argument('num_molecules')
    parser.add_argument('--flavour_idx', default=None)
    
    args = parser.parse_args()
    print("Creating placeholder categories ... ")
    create_placeholder_categories(args.flavour_file, args.molecule_xyz_file, args.output_temp_dir, int(args.num_molecules), args.flavour_idx)
    print("Done")