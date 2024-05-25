import copy
import json
import os
from multiprocessing.pool import Pool
import numpy as np

import utils.dft_utils as dft
import utils.xtb_utils as xtb
import utils.xyz_utils as xyz

def create_flavour_settings(base_settings, flavour_def): 
    task_settings = copy.deepcopy(base_settings)
    task_settings["turbomole_functional"] = flavour_def["functional"]
    task_settings["turbomole_basis"] = flavour_def["basisset"]

    return task_settings


def calculate_energies_for_category(path_to_flavour, base_settings, number_of_workers):
    """ Calculates energies for one category

    path_to_task (str): path to placeholder category
    settings (dict): settings for dft
    number_of_workers (int): number of workers (parallel calculations)
    """

    # load category (flavour) information from placeholder category dir
    with open(os.path.join(path_to_flavour, "task_info.json"), 'r') as fp:
        flavour_def = json.load(fp)

    # find all xyz files in the category
    xyz_files = []
    for x in os.listdir(path_to_flavour):
        if x.endswith(".xyz"):
            xyz_files.append(x)

    if len(xyz_files) != 1:
        raise NotImplementedError

    xyz_file_path = os.path.join(path_to_flavour, xyz_files[0])

    # read in the data about molecules
    coords_all, elements_all, charges_all = xyz.readXYZs_with_charges(xyz_file_path)
    assert len(coords_all) == len(elements_all) == len(charges_all)
    num_calcs = len(coords_all)

    # create dft settings 
    flavour_settings = create_flavour_settings(base_settings=base_settings, flavour_def=flavour_def)

    items = [(i, [coords_all[i], elements_all[i], charges_all[i], flavour_settings]) for i in range(num_calcs)]

    # create a temporary folder to store the calculations
    temp_flavour_folder = f"dft_tmpdirs_{flavour_settings['turbomole_functional']}_{flavour_settings['turbomole_basis']}"
    temp_flavour_folder_path = os.path.join(os.getcwd(), temp_flavour_folder)
    if not os.path.exists(temp_flavour_folder_path):
        os.makedirs(temp_flavour_folder_path)

    startdir = os.getcwd()
    os.chdir(temp_flavour_folder_path)
    energies, gradients = calc_energies_for_items(items, number_of_workers=number_of_workers, coords_all=coords_all)
    os.chdir(startdir)
    return energies, gradients


def calc_energies_for_items(items, number_of_workers, coords_all):
    """ Calculates the energies for items (molecules) in category

    items: items to calculate the energies for
    number_of_workers: number of workers (parallel calculations)
    coords_all: coordinates of all items (molecules)
    """
    with Pool(number_of_workers) as pool:
        # issues tasks to process pool
        results = pool.starmap_async(qm_task, items).get()

        # iterate results
        energies_all = []
        gradients_all = []
        for molidx, results_here in enumerate(results):
            print("Got result: {}".format(results_here["energy"]), flush=True)
            print(f"gradients have length {len(results_here['gradient'])}")
            # sanity check:
            coords_i = items[molidx][1][0]
            assert coords_all[molidx] == coords_i
            diff = np.array(results_here["coords"]) - np.array(coords_all[molidx])
            if np.max(np.abs(diff)) > 1e-5:
                print("WARNING: the coordinates of molecule {} do not agree with results".format(molidx))
                results_here["energy"] = None
                results_here["gradient"] = None
            energies_all.append(results_here["energy"])
            gradients_all.append(results_here["gradient"].tolist())
        
        pool.close()
        pool.join()
    # process pool is closed automatically
    return energies_all, gradients_all


def qm_task(identifier, data):

    # read it the data for one task
    coords = data[0]
    elements = data[1]
    charge = data[2]
    settings = data[3]
    
    if settings["qm_method"] == "xtb":
        results = xtb.xtb_calc(settings, coords, elements, opt=False, grad=False, hess=False, charge=0, freeze=[])
    elif settings["qm_method"] == "dft":
        results = dft.dft_calc(settings, coords, elements, charge, partial_chrg=settings['partial_chrg'], unp_el=settings['unp_el'], h20=settings['h20'])
    else:
        results = {}
    
    return (results)

