import pandas as pd
import argparse
import random
import json
import os
import re

import utils.xyz_utils as xyz

def create_placeholder_categories(flavour_file, xyz_files_dir, temp_cat_dir, num_molecules, flavour_idx=None):
    """ Creates the folders (placeholder categories) in which the data and settings for energy calculations are stored.
    Data sampling assumes that the dataset directory contains different datatypes in separate folders and each one of them contains: 
    - one xyz file containing multiple molecules 
    - one csv file containing their corresponding names (and/or other info) 
    To adjust the normal expression for finding those files, go to lines 96 and 97. 
    Such an approach guarantees balanced sampling of each data type (eg. amino acids, dipeptides, tripeptides).

    flavour_file (str): path to the flavour file - json file with a dictionary storing two lists: functinals and basis sets names (example: input_files/flavours.json)
    xyz_files_dir (str): path to a dataset to sample molecules from, containing each datatype in separate folders
    temp_cat_dir (str): path to store the output (placeholder categories) in
    num_molecules (int): number of molecules to sample per datatype per category (flavour)
    flavour_idx (int, optional): if defined, the code creates a placeholder category for only one flavour, specified with this number
    """

    # create the output directory if it doesn't exist yet
    if not os.path.exists(temp_cat_dir):
        os.makedirs(temp_cat_dir)

    # load the functionals and basis sets data
    with open(flavour_file, 'r') as f:
        func_and_basis = json.load(f)
    functionals = func_and_basis["functionals"]
    basissets = func_and_basis["basissets"]

    # define the flavours to skip (e.g. too computationally expensive ones)
    forbidden = [('bmk', 'aug-cc-pVDZ'), ('b3-lyp', 'aug-cc-pVDZ'), ('b3-lyp', '6-311++G**'),
                ('pbe0', 'aug-cc-pVDZ'), ('pbe0', '6-311++G**'), ('tpssh', 'aug-cc-pVDZ'),
                ('tpssh', '6-311++G**'), ('m06-2x', 'aug-cc-pVDZ'), ('m06-2x', '6-311++G**'),
                ('bmk', '6-311++G**'), ('bh-lyp', '6-311++G**'), ('bh-lyp', 'aug-cc-pVDZ')]

    # create the dictionary storing each flavour and its index
    dft_flavours = []
    index = 1
    for f in functionals:
        for b in basissets:
            if (f, b) not in forbidden:
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

    # list the datatypes in the dataset directory
    all_items = os.listdir(xyz_files_dir)
    datatypes_paths = [os.path.join(xyz_files_dir, item) for item in all_items 
                      if os.path.isdir(os.path.join(xyz_files_dir, item))]
    print(f"Found datatypes are {datatypes_paths}")

    # create placeholder categories for each flavour
    if flavour_idx==None:
        for single_flavour in dft_flavours:
            export_molecules_for_flavour(single_flavour, datatypes_paths, temp_cat_dir)

    # create a placeholder category only for flavour number = flavour_idx
    elif str(flavour_idx).isnumeric() and 1 <= int(flavour_idx) <= len(dft_flavours):
        single_flavour = dft_flavours[int(flavour_idx)-1]
        export_molecules_for_flavour(single_flavour, datatypes_paths, temp_cat_dir)
        
    else:
        raise ValueError(f"Flavour value should be an integer (index of the flavour to use), maximum value is {dft_flavours[-1]['number']}")

def export_molecules_for_flavour(single_flavour, datatypes_paths, temp_cat_dir):
    """ Samples the data from the dataset and saves the placeholder category.

    single_flavour (dict): contains the index, functional, basis set, and number of molecules to sample for a single flavour
    datatypes_paths (list): contains the paths to all data types in the dataset
    temp_cat_dir (str): path to store the output (placeholder category) in
    """

    datapoints_to_sample = int(single_flavour["num_molecules"])
    working_dir = os.getcwd()
    coords_all, elements_all, charges_all = [], [], []

    for datatype_path in datatypes_paths:

        # go into the folder containing certain data type and list all files inside
        os.chdir(datatype_path)
        all_files = os.listdir()
        
        # find the xyz and info file based on the normal expression (adjust in case of a different dataset!)
        file_pattern_coords = r'.*se1_ID_0_coords.xyz$' 
        file_pattern_info = r'.*se1_ID_0_info.csv$'         #change to se3 if necessary
        pattern_coords = re.compile(file_pattern_coords)
        pattern_info = re.compile(file_pattern_info)

        # read the coordinates and elements from the xyz file
        coords_file = [file for file in all_files if pattern_coords.match(file)][0]
        coords_all_1type, elements_all_1type = xyz.readXYZs(coords_file)
        assert len(coords_all) == len(elements_all)
        
        # read the ids and names of molecules from the info file
        info_file_path = [file for file in all_files if pattern_info.match(file)][0]
        info_df = pd.read_csv(info_file_path, usecols=['ID', 'names'])
        ids = info_df['ID'].tolist()
        system_names = info_df['names'].tolist()
        
        # sample random molecules
        num_all_mol = len(coords_all_1type)
        sampled_indices = random.sample(range(num_all_mol), datapoints_to_sample)

        # save the data of sampled molecules for that datatype
        coords_sampled = [coords_all_1type[i] for i in sampled_indices]
        elements_sampled = [elements_all_1type[i] for i in sampled_indices]
        charges_sampled = [find_mol_charge(i, ids, system_names) for i in sampled_indices]

        # append the info to all sampled molecules for flavour
        coords_all.extend(coords_sampled)
        elements_all.extend(elements_sampled)
        charges_all.extend(charges_sampled)

        # go back to the dataset directory
        os.chdir(working_dir)

    assert len(coords_all) == len(elements_all) == len(charges_all)
    print(f"Data was collected for {len(coords_all)/datapoints_to_sample} data types, in total: {len(coords_all)} molecules")

    # define the name of the data file (TODO -  meaningful index values)
    cat_coord_filename = f"data_01_{single_flavour['functional']}###{single_flavour['basisset']}.xyz"
     
    # create the directory storing molecules and category info file
    cat_dir_path = os.path.join(temp_cat_dir, f"flv{single_flavour['number']}_{single_flavour['functional']}###{single_flavour['basisset']}")
    if not os.path.exists(cat_dir_path):
        os.mkdir(cat_dir_path)

    # save the sampled molecules in xyz format
    xyz.exportXYZs_with_charges(coords_all, elements_all, charges_all, os.path.join(cat_dir_path, cat_coord_filename))

    #TODO - change the info file name to "category_info.json"
    # save the flavour settings within placeholder category
    with open(os.path.join(cat_dir_path, "task_info.json"), 'w') as fp:
        json.dump(single_flavour, fp)
        

def find_mol_charge(i, ids, system_names):
    ''' Determines the charge of the molecule based on its name containing + or - signs

    i (int): index of the molecule to be calculated
    ids (list): list of indexes from the info csv file
    system_names (list): list of names of the molecules from the info csv file
    '''
    system = system_names[i]
    pos_chr = system.count('+')
    neg_chr = system.count('-')
    tot_chr = pos_chr - neg_chr   

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