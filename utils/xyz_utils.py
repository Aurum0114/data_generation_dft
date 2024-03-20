import numpy as np
import os
import subprocess

kcal_to_eV=0.0433641153
kB=8.6173303e-5 #eV/K
T=298.15
kBT=kB*T
AToBohr=1.889725989
HToeV = 27.211399

#converts XYZ files into turbomole coordinates
def x2t_command(infile, outfile, moldir):
    startdir = os.getcwd() #gets the current working directory
    os.chdir(moldir) #changes the directory to moldir
    #accesses Turbomole, running the conversion command
    subprocess.run("x2t {} > {}".format(infile, outfile), shell=True) 
    os.chdir(startdir) #goes back to the firs working directory

def readXYZ(filename):
    infile=open(filename,"r")
    coords=[]
    elements=[]
    lines=infile.readlines()
    if len(lines)<3:
        exit("ERROR: no coordinates found in %s/%s"%(os.getcwd(), filename))
    for line in lines[2:]:
        elements.append(line.split()[0].capitalize())
        coords.append([float(line.split()[1]),float(line.split()[2]),float(line.split()[3])])
    infile.close()
    coords=np.array(coords)
    return coords,elements


def readXYZs(filename):
    infile=open(filename,"r")
    coords=[[]]
    elements=[[]]
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

#Adrian's one
def exportXYZ(coords, elements, filename, mask=[]):
    outfile=open(filename, "w")

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

#Marlen's one
def exportXYZs(coords, elements, charges, filename):
    outfile = open(filename, "a")
    
    for i in range(len(coords)):
        outfile.write(f"{len(elements[i])} {charges[i]}\n\n")

        for atomidx, atom in enumerate(coords[i]):
            print(f"total coords list for i molecule is is {coords[i]}")
            print(f"list of elements for i is {elements[i]}")
            print(f'finally, the atomidx element of that list is {elements[i][atomidx]}')
            outfile.write(f"{elements[i][atomidx].capitalize()} {atom[0]} {atom[1]} {atom[2]}\n")

    outfile.close()

#filename = name of XYZ output file
def a2bohr_exportXYZ(coords, elements, filename):
    #coords = coords*AToBohr
    outfile = open(filename, "w") #creates a new file or truncates an existing one
    outfile.write("%i\n\n"%(len(elements))) #first the number of elements, then a blank line
    for atomidx,atom in enumerate(coords):
        outfile.write("%s %f %f %f\n"%(elements[atomidx].capitalize(), atom[0]* AToBohr, atom[1]* AToBohr, atom[2]* AToBohr))
    outfile.close()