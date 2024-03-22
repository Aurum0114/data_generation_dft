import copy
import json
import os
from multiprocessing.pool import Pool
import numpy as np

import utils.dft_utils as dft
import utils.xtb_utils as xtb
import utils.xyz_utils as xyz

#dft_settings = {"copy_mos": False,
#                "use_dispersions": True,
#                "turbomole_method": "ridft",
#                "turbomole_basis": "6-311++G**", #  def2-SV(P)  6-311++G**
#                "turbomole_functional": "bmk"} #  BMK?? b3-lyp

def create_flavour_setting(base_settings, flavour_def): 
    task_settings = copy.deepcopy(base_settings)
    task_settings["turbomole_functional"] = flavour_def["functional"]
    task_settings["turbomole_basis"] = flavour_def["basisset"]

    return task_settings


def calculate_energies_for_task(path_to_task, base_settings, number_of_workers):
    """
    Function that calculates energies for placeholder categories
    :param path_to_task: path to placeholder category
    :param settings: settings for dft
    :param number_of_workers: number of workers
    :return:
    """
    print("Calculate energies for task in path: ", path_to_task)
    # load flavour information from placeholder category dir
    with open(os.path.join(path_to_task, "task_info.json"), 'r') as fp:
        flavour_def = json.load(fp)

    xyz_files = []
    for x in os.listdir(path_to_task):
        if x.endswith(".xyz"):
            xyz_files.append(x)

    if len(xyz_files) != 1:
        raise NotImplementedError

    xyz_file_path = os.path.join(path_to_task, xyz_files[0])

    coords_all, elements_all, charges_all = xyz.readXYZs_with_charges(xyz_file_path)
    assert len(coords_all) == len(elements_all) == len(charges_all)
    num_calcs = len(coords_all)
    print("The number of calculations to be performed is:", num_calcs)

    task_settings = create_flavour_setting(base_settings=base_settings, flavour_def=flavour_def)
    print("Task settings are equal to:", task_settings)

    items = [(i, [coords_all[i], elements_all[i], charges_all[i], task_settings]) for i in range(num_calcs)]
    print("Provided items are: ", items)

    temp_flavour_folder = f"dft_tmpdirs_{task_settings['turbomole_functional']}_{task_settings['turbomole_basis']}"
    temp_flavour_folder_path = os.path.join(os.getcwd(), temp_flavour_folder)
    if not os.path.exists(temp_flavour_folder_path):
        os.makedirs(temp_flavour_folder_path)

    startdir = os.getcwd()
    os.chdir(temp_flavour_folder_path)
    energies = calc_energies_for_items(items, number_of_workers=number_of_workers, coords_all=coords_all)
    os.chdir(startdir)
    return energies


def calc_energies_for_items(items, number_of_workers, coords_all):
    """

    :param items: Items to calculate the energies for
    :param number_of_workers:
    :param coords_all:
    :return:
    """
    with Pool(number_of_workers) as pool:
        # issues tasks to process pool
        results = pool.starmap_async(qm_task, items).get()

        # iterate results
        energies_all = []
        # gradients_all = []
        for molidx, results_here in enumerate(results):
            print("Got result: {}".format(results_here["energy"]), flush=True)
            # sanity check:
            coords_i = items[molidx][1][0]
            assert coords_all[molidx] == coords_i
            diff = np.array(results_here["coords"]) - np.array(coords_all[molidx])
            if np.max(np.abs(diff)) > 1e-5:
                print("WARNING: the coordinates of molecule {} do not agree with results".format(molidx))
                results_here["energy"] = None
                results_here["gradient"] = None
            energies_all.append(results_here["energy"])
            # gradients_all.append(results_here["gradient"].tolist())
        
        pool.close()
        pool.join()
        print("The pool finished, yielding ", energies_all)
    # process pool is closed automatically
    return energies_all


def qm_task(identifier, data):
    print("Calculating task number: ", identifier)
    coords = data[0]
    print("Provided coordinates for task number", identifier, " are:", coords)
    elements = data[1]
    print("Provided elements for task number", identifier, " are:", elements)
    charge = data[2]
    print(f"Provided charge for task {identifier} is {charge}")
    settings = data[3]
    print("Provided settings for task number", identifier, " are:", settings)
    
    
    if settings["qm_method"] == "xtb":
        results = xtb.xtb_calc(settings, coords, elements, opt=False, grad=False, hess=False, charge=0, freeze=[])
    elif settings["qm_method"] == "dft":
        results = dft.dft_calc(settings, coords, elements, charge, partial_chrg=settings['partial_chrg'], unp_el=settings['unp_el'], h20=settings['h20'])
    else:
        results = {}
    
    print("Qm task number ", identifier, " finished with results: ", results)
    
    return (results)

