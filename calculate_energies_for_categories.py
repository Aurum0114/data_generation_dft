import argparse
import os
import shutil
import numpy as np
import json
from parallel_qm import calculate_energies_for_task


def calculate_energies_for_categories(flavours_dir, results_dir, num_workers):
    """
    :param temp_dir: Path to the directory storing the temporary placeholder categories
    :param output_dir: Where to store the dataset
    :param num_workers: The number of workers / number of parallel executions
    """

    if not os.path.exists(results_dir):
        os.makedirs(results_dir)

    base_settings = {"qm_method": "dft",
                     "delete_calculation_dirs": False,
                     "copy_mos": False,
                     "dispersion": True,
                     "turbomole_method": "ridft",
                     "turbomole_basis": "6-311++G**",
                     "turbomole_functional": "bmk",
                     "partial_chrg": False,
                     "unp_el": 1,
                     "h20": True,
                     }
    
    all_flavours_todo_paths = find_all_task_dirs(flavours_dir)
    print("All of the task dir to be calculated are:", all_flavours_todo_paths)

    for flavour_todo in all_flavours_todo_paths:
        fl_num = flavour_todo[4:5]
        print("Flavour to be calculated directory ", flavour_todo)
        flavour_todo_path = os.path.join(flavours_dir, flavour_todo)
        flavour_done_path = os.path.join(results_dir, flavour_todo)

        #calculates energies here
        energies = calculate_energies_for_task(path_to_task=flavour_todo_path,
                                               base_settings=base_settings,
                                               number_of_workers=num_workers)
        
        results_file_name = f"final_energies_{fl_num}.npy"
        results_file_path = os.path.join(flavour_done_path, results_file_name)
        results_file_path_in_flavours_dir = os.path.join(flavour_todo_path, results_file_name)

        if os.path.exists(results_file_path):
            print(f"File {results_file_name} has been found in {flavour_done_path}, appending the results there...")

            update_molecules_and_task_info(flavour_todo_path, flavour_done_path)

            existing_energies = np.load(results_file_path)
            energies = np.array(energies)
            energies = np.concatenate((existing_energies, energies))
            np.save(results_file_path, energies)
            
            for item in os.listdir(flavour_todo_path):
                print(f"moving {item} to done path")
                item_path = os.path.join(flavour_todo_path, item)
                shutil.move(item_path, flavour_done_path)
            
            print(f"Removing flavour todo path {flavour_todo_path}")
            os.rmdir(flavour_todo_path)

        else:
            print(f"File {results_file_name} has not been found. Creating a new one...")
            np.save(results_file_path_in_flavours_dir, energies)
            shutil.move(flavour_todo_path, flavour_done_path)
    
    print(f"All calculations finished!")

def update_molecules_and_task_info(source_dir, destination_dir):

    existing_molecules_path = find_file_path(destination_dir, '.xyz')
    print(f"existing molecules path is {existing_molecules_path}")
    molecules_to_append_path = find_file_path(source_dir, '.xyz')
    existing_info_path = find_file_path(destination_dir, 'info.json')
    info_to_append_path = find_file_path(source_dir, 'info.json')

    with open(molecules_to_append_path, 'r') as src_file:
        molecules_to_append = src_file.read()
    with open(existing_molecules_path, 'a') as dest_file:
        dest_file.write(molecules_to_append)

    with open(info_to_append_path, 'r') as src_file:
        info_to_append = json.load(src_file)
        functional = info_to_append['functional']
        basisset = info_to_append['basisset']
        num_mol_to_append = info_to_append['num_molecules']
    
    with open(existing_info_path, 'r') as dest_file:
        info_to_update = json.load(dest_file)
        assert info_to_update['functional'] == functional, f"Functionals mismatch: source has {function}, whereas desination has {info_to_update['functional']}"
        assert info_to_update['basisset'] == basisset, f"Basis sets mismatch: source has {basisset}, whereas desination has {info_to_update['basisset']}"
        info_to_update['num_molecules'] = info_to_update['num_molecules'] + num_mol_to_append

    with open(os.path.join(existing_info_path), 'w') as new_dest_file:
        json.dump(info_to_update, new_dest_file)
    
    os.remove(existing_molecules_path)
    os.remove(info_to_append_path)


def find_file_path(path_to_search, norm_expression):
    target_file_path = None
    for file in os.listdir(path_to_search):
        if file.endswith(norm_expression):
            target_file_path = os.path.join(path_to_search, file)
            print(f"file {target_file_path} was found in {path_to_search}")
            break
    if target_file_path is None:
        print(f"No files with normal expression {norm_expression} were found in {path_to_search}")
    return target_file_path


def find_all_task_dirs(path_to_tasks): 
    all_task_dirs = []
    for x in os.listdir(path_to_tasks):
        if os.path.isdir(os.path.join(path_to_tasks, x)) and x.startswith("FLV_"):
            all_task_dirs.append(x)
    all_task_dirs = sorted(all_task_dirs)
    return all_task_dirs


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('temp_dir')
    parser.add_argument('output_dir')
    parser.add_argument('num_workers')
    args = parser.parse_args()
    print("Calculating energies ...")
    print("main function arguments are: ", args.temp_dir, args.output_dir, args.num_workers)
    calculate_energies_for_categories(args.temp_dir, args.output_dir, int(args.num_workers))
    print("Done")