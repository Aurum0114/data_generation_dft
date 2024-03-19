import argparse
import json
import os
import random
import shutil

import utils.xyz_utils as xyz

def create_placeholder_categories(flavour_file, xyz_files_dir, output_temp_dir, num_molecules, flavour_idx=None):
    """
    This function creates placeholder categories that will be used by the next function to create the actual categories.
    :flavour_file: Stores two lists: all functionals and basis sets. Example file: example_files/func_and_base.json
    :xyz_files_dir: Filtered dataset to sample molecules from. Example file: example_files/inputs.xyz
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


    datatypes = find_all_subdirs(xyz_files_dir)
    #coords_all, elements_all = xyz.readXYZs(molecule_xyz_file)
    #assert len(coords_all) == len(elements_all)

    if flavour_idx==None:
        for single_flavour in dft_flavours:
            sample_molecules_for_flavour(single_flavour, datatypes, temp_task_dir)

    elif str(flavour_idx).isnumeric() and 1 <= int(flavour_idx) <= len(dft_flavours):
        single_flavour = dft_flavours[flavour_idx-1]
        sample_molecules_for_flavour(single_flavour, datatypes, temp_task_dir)
        
    else:
        raise ValueError(f"Flavour value should be an integer (index of the flavour to use), maximum value is {dft_flavours[-1]['number']}")
    
def find_all_subdirs(xyz_files_dir):
    all_items = os.listdir(xyz_files_dir)
    subdirectories = [item for item in all_items if os.path.isdir(os.path.join(xyz_files_dir, item))]
    return subdirectories

def sample_molecules_for_flavour(single_flavour, datatypes, temp_task_dir):
    #num_all_mol = len(coords_all)
    datapoints_to_sample = int(single_flavour["num_molecules"])
    coords = []
    elements = []

    for datatype in datatypes:
        sample_molecules_for_datatype(coords, elements, datapoints_to_sample)

    #sampled_indices = random.sample(range(num_all_mol), datapoints_to_sample)
    #coords = [coords_all[i] for i in sampled_indices]
    #elements = [elements_all[i] for i in sampled_indices]

    task_coord_filename = f"data_01_{single_flavour['functional']}###{single_flavour['basisset']}.xyz"
    task_dir_path = os.path.join(temp_task_dir,
                                f"FLV_{single_flavour['number']}_{single_flavour['functional']}###{single_flavour['basisset']}")

    os.mkdir(task_dir_path)
    print(f"Creating directory {task_dir_path}")
    xyz.exportXYZs(coords, elements, os.path.join(task_dir_path, task_coord_filename))

    with open(os.path.join(task_dir_path, "task_info.json"), 'w') as fp:
        json.dump(single_flavour, fp)

def sample_molecules_for_datatype(coords, elements, datatype, datapoints_to_sample):
    #go into that directory
    #find the coords file
    #find the info file

    #read out the coords file
    #change reading function so that it adds the charge
    #coords_all, elements_all = xyz.readXYZs(molecule_xyz_file)
    
    #find number of molecules
    #sampled_indices = random.sample(range(num_all_mol), datapoints_to_sample)
    #list compreh for coords and elements
    #append to the total coords and elements lists
    pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('flavour_file')
    parser.add_argument('molecule_xyz_file')
    parser.add_argument('output_temp_dir')
    parser.add_argument('num_molecules')
    parser.add_argument('--flavour_idx', default=None)
    
    args = parser.parse_args()
    print("Creating placeholder categories ... ")
    create_placeholder_categories(args.flavour_file, args.molecule_xyz_file, args.output_temp_dir, int(args.num_molecules), int(args.flavour_idx))
    print("Done")