import numpy as np
import argparse
import shutil
import json
import os

from parallel_qm import calculate_energies_for_task


def calculate_energies_for_categories(flavours_dir, results_dir, num_workers):
    """ Calculates the energies for provided molecules using parallel Turbomole (DFT) calculations.

    flavours_dir (str): path to the directory storing the temporary placeholder categories
    results_dir (str): path to store the dataset in
    num_workers (int): the number of workers / number of parallel executions 
    """

    # create the results directory if it doesn't exist
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)

    base_settings = {"qm_method": "dft",
                     "turbomole_method": "ridft",
                     "turbomole_basis": "6-311++G**",
                     "turbomole_functional": "bmk",
                     "delete_calculation_dirs": False,
                     "copy_mos": False,
                     "partial_chrg": False,
                     "unp_el": 1,
                     "h20": True,
                     }
    
    # finds all directories (different flavours) in the flavours directory
    all_flavours_todo_paths = find_all_category_dirs(flavours_dir)
    print("All flavours to be calculated are:", all_flavours_todo_paths)

    for flavour_todo in all_flavours_todo_paths:

        #fl_num = flavour_todo[3:4] - to be implemented later, for now the name is fixed with 01
        flavour_todo_path = os.path.join(flavours_dir, flavour_todo)
        flavour_done_path = os.path.join(results_dir, flavour_todo)

        # calculate energies
        print("Calculating energies...")
        energies = calculate_energies_for_task(path_to_task=flavour_todo_path,
                                               base_settings=base_settings,
                                               number_of_workers=num_workers)
        
        # define the results names and paths
        results_file_name = f"final_energies_01.npy"
        results_file_path = os.path.join(flavour_done_path, results_file_name)
        results_file_path_in_flavours_dir = os.path.join(flavour_todo_path, results_file_name)

        # if results file exists, append the results there
        if os.path.exists(results_file_path):
            print(f"File {results_file_name} has been found in {flavour_done_path}, appending the results there...")

            update_molecules_and_task_info(flavour_todo_path, flavour_done_path)
            existing_energies = np.load(results_file_path)
            energies = np.array(energies)
            energies = np.concatenate((existing_energies, energies))
            np.save(results_file_path, energies)
            
            # if there are any leftover files, move them to the results directory
            for item in os.listdir(flavour_todo_path):
                item_path = os.path.join(flavour_todo_path, item)
                shutil.move(item_path, flavour_done_path)
            
            # remove the placeholder category directory
            os.rmdir(flavour_todo_path)

        else:
            # save the results in a new file
            print(f"File {results_file_name} has not been found. Creating a new one...")
            np.save(results_file_path_in_flavours_dir, energies)
            shutil.move(flavour_todo_path, flavour_done_path)
    
    print(f"All calculations finished!")

def update_molecules_and_task_info(source_dir, destination_dir):
    ''' Appends the data (molecular coordinates and the total number of calculated molecules) to the final results files

    source_dir (str): path of the directory from which the data will be copied
    destination_dir (str): path of the directory where the data will be appended
    '''

    # find correct molecules and task files
    existing_molecules_path = find_file_path(destination_dir, '.xyz')
    molecules_to_append_path = find_file_path(source_dir, '.xyz')
    existing_info_path = find_file_path(destination_dir, 'info.json')
    info_to_append_path = find_file_path(source_dir, 'info.json')

    # read in the molecules data and append to the existing file
    with open(molecules_to_append_path, 'r') as src_file:
        molecules_to_append = src_file.read()
    with open(existing_molecules_path, 'a') as dest_file:
        dest_file.write(molecules_to_append)

    # extract the flavour information from info file
    with open(info_to_append_path, 'r') as src_file:
        info_to_append = json.load(src_file)
        functional = info_to_append['functional']
        basisset = info_to_append['basisset']
        num_mol_to_append = info_to_append['num_molecules']
    
    # assert the flavour information is the same in both info files
    with open(existing_info_path, 'r') as dest_file:
        info_to_update = json.load(dest_file)
        assert info_to_update['functional'] == functional, f"Functionals mismatch: source has {function}, whereas desination has {info_to_update['functional']}"
        assert info_to_update['basisset'] == basisset, f"Basis sets mismatch: source has {basisset}, whereas desination has {info_to_update['basisset']}"

        # update the total number of calculated molecules
        info_to_update['num_molecules'] = info_to_update['num_molecules'] + num_mol_to_append

    with open(os.path.join(existing_info_path), 'w') as new_dest_file:
        json.dump(info_to_update, new_dest_file)
    
    # remove both source files
    os.remove(molecules_to_append_path)
    os.remove(info_to_append_path)


def find_file_path(path_to_search, norm_expression):
    ''' Finds the first file with provided normal expression in the provided directory and returns its path

    path_to_search (str): directory where the file is searched for
    norm_expression (str): normal expression to search for
    '''
    target_file_path = None
    for file in os.listdir(path_to_search):
        if file.endswith(norm_expression):
            target_file_path = os.path.join(path_to_search, file)
            break
    if target_file_path is None:
        print(f"No files with normal expression {norm_expression} were found in {path_to_search}")
    return target_file_path


def find_all_category_dirs(path_to_categories):
    ''' Finds all of the folders that start with "flv" in the provided directory
    
    path_to_categories (str): directory where the flavour directories are stored
    '''
    all_task_dirs = []
    for x in os.listdir(path_to_categories):
        if os.path.isdir(os.path.join(path_to_categories, x)) and x.startswith("flv"):
            all_task_dirs.append(x)
    all_task_dirs = sorted(all_task_dirs)
    return all_task_dirs


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('temp_dir')
    parser.add_argument('output_dir')
    parser.add_argument('num_workers')
    args = parser.parse_args()
    calculate_energies_for_categories(args.temp_dir, args.output_dir, int(args.num_workers))