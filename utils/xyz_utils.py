import numpy as np
import subprocess
import os

kcal_to_eV=0.0433641153
kB=8.6173303e-5 #eV/K
AToBohr=1.889725989
HToeV = 27.211399
T=298.15
kBT=kB*T

#converts XYZ files into turbomole coordinates
def x2t_command(infile, outfile, moldir):
    startdir = os.getcwd()
    os.chdir(moldir) 

    # run the Turbomole the conversion command
    subprocess.run("x2t {} > {}".format(infile, outfile), shell=True) 
    os.chdir(startdir) 

# for normal xyz files
def readXYZs(filename):
    ''' Reads multiple molecules from a xyz file. '''
    infile = open(filename,"r")
    coords = [[]]
    elements = [[]]
    for line in infile.readlines():
        words = line.split()

        if len(words)==1 and str(words[0]).isnumeric() and len(coords[-1])!=0:
            coords.append([])
            elements.append([])
        elif len(line.split())==4:
            elements[-1].append(line.split()[0].capitalize())
            coords[-1].append([float(line.split()[1]),float(line.split()[2]),float(line.split()[3])])

    infile.close()
    return coords,elements

def readXYZs_with_charges(filename):
    ''' Reads multiple molecules from a xyz file containing the charge of the molecules next to the atoms count. '''

    infile=open(filename,"r")
    coords, elements, charges = [], [], []

    for line in infile.readlines():
        words = line.split()
        
        if len(words)==2 and str(words[0]).isnumeric():
            coords.append([])
            elements.append([])
            charges.append(words[1])
        elif len(words)==4:
            elements[-1].append(words[0].capitalize())
            coords[-1].append([float(line.split()[1]), float(line.split()[2]), float(line.split()[3])])

    infile.close()
    return coords, elements, charges

def exportXYZ(coords, elements, filename, mask=[]):
    ''' Saves the information about molecules in one xyz file'''

    outfile = open(filename, "w")

    if len(mask)==0:
        outfile.write("%i\n\n"%(len(elements)))
        for atomidx,atom in enumerate(coords):
            outfile.write("%s %f %f %f\n"%(elements[atomidx].capitalize(),atom[0],atom[1],atom[2]))
    else:
        outfile.write("%i\n\n"%(len(mask)))
        for atomidx in mask:
            atom = coords[atomidx]
            outfile.write("%s %f %f %f\n"%(elements[atomidx].capitalize(),atom[0],atom[1],atom[2]))
    outfile.close()


def exportXYZs_with_charges(coords, elements, charges, filename):
    ''' Saves the information about molecules in one xyz file containing charge value next to the atoms count'''

    outfile = open(filename, "a")
    
    for i in range(len(coords)):
        outfile.write(f"{len(elements[i])} {charges[i]}\n\n")

        for atomidx, atom in enumerate(coords[i]):
            outfile.write(f"{elements[i][atomidx].capitalize()} {atom[0]} {atom[1]} {atom[2]}\n")

    outfile.close()

def export_forces(forces, elements, filename):
    ''' Saves the information about forces in one xyz file '''
    print(f"exporting {len(forces)} forces for {len(elements)} elements")

    outfile = open(filename, "a")
    assert len(forces) == len(elements), "more forces than molecules are returned"

    for i in range(len(forces)):
        outfile.write(f"{len(elements[i])} \n")

        for atomidx, force in enumerate(forces[i]):
            outfile.write(f"{elements[i][atomidx].capitalize()} {force[0]} {force[1]} {force[2]}\n")
        outfile.write("\n")

    outfile.close()

#for the xtb calculations
def readXYZ(filename):
    infile = open(filename,"r")
    coords, elements = [], []
    lines = infile.readlines()

    if len(lines)<3:
        exit("ERROR: no coordinates found in %s/%s"%(os.getcwd(), filename))
    for line in lines[2:]:
        elements.append(line.split()[0].capitalize())
        coords.append([float(line.split()[1]),float(line.split()[2]),float(line.split()[3])])

    infile.close()
    coords = np.array(coords)
    return coords,elements

def a2bohr_exportXYZ(coords, elements, filename):
    #coords = coords*AToBohr
    outfile = open(filename, "w") 
    outfile.write("%i\n\n"%(len(elements)))
    for atomidx,atom in enumerate(coords):
        outfile.write("%s %f %f %f\n"%(elements[atomidx].capitalize(), atom[0]* AToBohr, atom[1]* AToBohr, atom[2]* AToBohr))
    outfile.close()