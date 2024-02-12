import os
import sys
import numpy as np
import scipy.spatial as scsp
import resource
import time
import random
import subprocess
import shutil
import uuid
import getpass
import socket
import yaml
import warnings

import utils.dft_utils as dft

kcal_to_eV = 0.0433641153
kB = 8.6173303e-5  # eV/K
T = 298.15
kBT = kB * T
AToBohr = 1.889725989
HToeV = 27.211399

def check_basis_and_func(basis_todo, func_todo, path_to_control):
    func = None
    basis = None
    for line in open(path_to_control, "r"):
        if "functional" in line:
            func = line.split()[1]
        elif "basis =" in line and basis is None:
            basis = line.split()[2]
    if basis is None:
        warnings.warn(f"Warning: No basis found in control file! (Path to control: {path_to_control})")
        #exit(f"No basis found in control file! (Path to control: {path_to_control})")
    if func is None:
        warnings.warn(f"Warning: No functional found in control file! (Path to control: {path_to_control})")
        #exit(f"No functional found in control file! (Path to control: {path_to_control})")
    if basis != basis_todo:
        warnings.warn(f"Warning: Wrong basis in control file: Expected {basis_todo} but found {basis}!"
                      f" (Path to control: {path_to_control})")
        #exit(f"Wrong basis in control file: Expected {basis_todo} but found {basis}! (Path to control: {path_to_control})")
    if func != func_todo:
        warnings.warn(f"Warning: Wrong functional in control file: Expected {func_todo} but found {func}! "
                      f"(Path to control: {path_to_control})")
        #exit("Wrong functional is control file: Expected %s but found %s" % (func_todo, func))

    print(f"basis and functional seems to be correct for {path_to_control}")
   

def try_mkdir(dirname):
    if not os.path.exists(dirname):
        try:
            os.makedirs(dirname)
        except:
            pass


def copy_dir_contents_to_dir(in_directory, out_directory, dir_to_copy):
    try:
        # We copy all files in in_directory to out_directory without following directories recursively:
        file_with_dir = "%s/%s" % (in_directory, dir_to_copy)
        if os.path.exists(file_with_dir):
            shutil.copytree(file_with_dir, "%s/%s" % (out_directory, dir_to_copy))
        else:
            print("did not find the directory %s to copy back from scratch" % (file_with_dir))
            exit()
    except Exception as exc:
        print("Moving files to from %s to %s has failed. Reraising Exception:" % (in_directory, out_directory))
        print(exc)
        raise

def GoToScratch():
    oldcwd = os.getcwd()
    try:
        SCRATCH_BASE = os.environ["SCRATCH"]
        username = getpass.getuser()
        randstring = uuid.uuid4()
        scratch_directory = "%s/%s/%s" % (SCRATCH_BASE, username, randstring)
        os.makedirs(scratch_directory)
    except KeyError as exc:
        print(
            "A KeyError occured, when querying the Scratch Directory. Check the environment settings. Exception was: %s. Turning off scratch handling." % (
                exc))
        scratch_directory = oldcwd
    except Exception as exc:
        # In case there was something unforseen, we reraise to bomb out of the application.
        print("An unexpected exception as occured of type %s. Exception was: %s. Reraising." % (type(exc), exc))
        raise
    os.chdir(scratch_directory)
    return ([oldcwd, scratch_directory])

def ComeBachFromScratch(oldcwd, scratch_directory, dir_to_copy):
    if oldcwd != scratch_directory:
        # copy result back to oldcwd, change back, remove scratch
        copy_dir_contents_to_dir(scratch_directory, oldcwd, dir_to_copy)
        os.chdir(oldcwd)
        shutil.rmtree(scratch_directory)
    else:
        print("Warning, oldcwd ", oldcwd, "was equal to scratch_directory", scratch_directory,
              "review log for exceptions.")

def read_dft_grad():
    if not os.path.exists("gradient"):
        return (None)
    grad = []
    for line in open("gradient", "r"):
        if len(line.split()) == 3 and "grad" not in line:
            line = line.replace("D", "E")
            grad.append([float(line.split()[0]), float(line.split()[1]), float(line.split()[2])])
    if len(grad) == 0:
        grad = None
    else:
        grad = np.array(grad) * HToeV * AToBohr
    return (grad)

def read_dft_hess():
    hess = None
    if not os.path.exists("hessian"):
        return (None, None, None)
    hess = []
    for line in open("hessian", "r"):
        if "hess" not in line:
            for x in line.split():
                hess.append(float(x))
    if len(hess) == 0:
        hess = None
    else:
        hess = np.array(hess)

    vibspectrum = None
    if not os.path.exists("vibspectrum"):
        return (None, None, None)
    vibspectrum = []
    read = False
    for line in open("vibspectrum", "r"):
        if "end" in line:
            read = False

        if read:
            if len(line.split()) == 5:
                vibspectrum.append(float(line.split()[1]))
            elif len(line.split()) == 6:
                vibspectrum.append(float(line.split()[2]))
            else:
                print("WARNING: weird line length: %s" % (line))
        if "RAMAN" in line:
            read = True

    reduced_masses = None
    if not os.path.exists("g98.out"):
        print("g98.out not found")
        return (None, None, None)
    reduced_masses = []
    read = False
    for line in open("g98.out", "r"):
        if "Red. masses" in line:
            for x in line.split()[3:]:
                try:
                    reduced_masses.append(float(x))
                except:
                    pass

    if len(vibspectrum) == 0:
        vibspectrum = None
        print("no vibspectrum found")
    else:
        vibspectrum = np.array(vibspectrum)

    if len(reduced_masses) == 0:
        reduced_masses = None
        print("no reduced masses found")
    else:
        reduced_masses = np.array(reduced_masses)

    return (hess, vibspectrum, reduced_masses)

def RunTMRelaxation(moldir, dft_settings):
    startdir = os.getcwd()
    os.chdir(moldir)

    instring = dft.prep_define_file(dft_settings, 0)
    dft.ExecuteDefineString(instring)

    if dft_settings["use_dispersions"]:
        dft.AddStatementToControl("control", "$disp3")

    if dft_settings["turbomole_method"] == "ridft":
        os.system("jobex -ri -c 200 > jobex.out")
    elif dft_settings["turbomole_method"] == "dscf":
        os.system("jobex -c 200 > jobex.out")
    else:
        exit("ERROR in turbomole_method: %s" % (dft_settings["turbomole_method"]))

    if os.path.exists("GEO_OPT_CONVERGED"):
        finished = True
        os.system("eiger > eiger.out")
    else:
        finished = False

    dft.AddStatementToControl("control", "$esp_fit kollman")
    if dft_settings["turbomole_method"] == "ridft":
        os.system("ridft -proper > TM_proper.out")
    elif dft_settings["turbomole_method"] == "dscf":
        os.system("dscf -proper > TM_proper.out")

    os.chdir(startdir)
    return (finished)

def getTMpartialcharges(outfilename, noOfAtoms):
    TMfile = open(outfilename, "r")
    lines = TMfile.readlines()
    TMfile.close()
    # finding position of ESP partial charges in file
    idx = 0
    for line in lines:
        spl = line.split()
        if len(spl) != 0 and len(spl) != 1:
            if spl[0] == "atom" and spl[1] == "radius/au":
                index = idx + 1
                break
        idx += 1
    partialCharges = []
    for i in range(noOfAtoms):
        partialCharges.append(float(lines[idx + 1 + i].split()[3]))
    return (partialCharges)

def getTMCoordinates(moldir, startOrEnd):
    infile = open("%s/gradient" % (moldir))
    lines = infile.readlines()
    infile.close()
    coordinates_all = []
    eles = []
    for idx, line in enumerate(lines):
        if "cycle =      2" in line:
            number_of_atoms = (idx - 2) / 2
            break
    print("number of atoms: %i" % (number_of_atoms))

    if len(lines) % (number_of_atoms * 2 + 1) != 2:
        print("WARNING")
        exit()

    number_of_steps = (len(lines) - 2) / (number_of_atoms * 2 + 1)
    print("found %i gradient steps" % (number_of_steps))

    for step in range(0, number_of_steps):
        coordinates_all.append([])
        eles.append([])
        startline = 1 + step * (number_of_atoms * 2 + 1) + 1
        for line in lines[startline:startline + number_of_atoms]:
            coordinates_all[step].append(
                [float(line.split()[0]) / AToBohr, float(line.split()[1]) / AToBohr, float(line.split()[2]) / AToBohr])
            eles[step].append(line.split()[3])

    if startOrEnd == "end":
        coordinates = coordinates_all[-1]
        elements = eles[-1]
    elif startOrEnd == "start":
        coordinates = coordinates_all[0]
        elements = eles[0]

    return (coordinates, elements)

def PrepTMInput(moldir, coords, elements, dihedral, dft_settings):
    frozen_atoms = dihedral[0]
    coordfile = open("%s/coord" % (moldir), 'w')
    coordfile.write("$coord \n")
    for idx, atom in enumerate(coords):
        if idx in frozen_atoms:
            frozen = "  f"
        else:
            frozen = ""
        coordfile.write(
            "%f  %f  %f  %s%s\n" % (atom[0] * AToBohr, atom[1] * AToBohr, atom[2] * AToBohr, elements[idx], frozen))
    coordfile.write("$end\n")
    coordfile.close()

    return ()