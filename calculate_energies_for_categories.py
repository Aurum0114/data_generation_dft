import argparse
import os
import shutil
import numpy as np
from parallel_qm import find_all_task_dirs, calculate_energies_for_task


def calculate_energies_for_categories(temp_dir, output_dir, num_workers):
    """
    :param temp_dir: Path to the directory storing the temporary placeholder categories
    :param output_dir: Where to store the dataset
    :param num_workers: The number of workers / number of parallel executions
    """
    path_to_temp_tasks = os.path.join(temp_dir, "tasks")
    path_to_finished_tasks = os.path.join(output_dir, "tasks")

    if not os.path.exists(path_to_finished_tasks):
        os.makedirs(path_to_finished_tasks)

    base_settings = {"qm_method": "dft",
                     "delete_calculation_dirs": False,
                     "copy_mos": False,
                     "use_dispersions": True,
                     "turbomole_method": "ridft",
                     "partial_chrg": False,
                     "unp_el": 1,
                     "h20": False,
                     }

    all_todo_task_dirs = find_all_task_dirs(path_to_temp_tasks)
    print("All of the task dir to be calculated are:", all_todo_task_dirs)

    for task_dir in all_todo_task_dirs:
        print("Task to be calculated directory ", task_dir)
        path_task_todo_dir = os.path.join(path_to_temp_tasks, task_dir)
        path_task_done_dir = os.path.join(path_to_finished_tasks, task_dir)
        energies = calculate_energies_for_task(path_to_task=path_task_todo_dir,
                                               settings=base_settings,
                                               number_of_workers=num_workers)
        
        output_file_name = f"labels_01_energies.npy"
        print("Saving the output files in path: ", path_task_todo_dir, '/', output_file_name)
        np.save(os.path.join(path_task_todo_dir, output_file_name), energies)

        # move task to done
        print("Moving task to done, from ", path_task_todo_dir, " to ", path_task_done_dir)
        shutil.move(path_task_todo_dir, path_task_done_dir)


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