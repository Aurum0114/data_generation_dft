#!/usr/bin/env python
import os
import numpy as np
from io import StringIO
import resource
import subprocess
import datetime
import uuid

import utils.xyz_utils as xyz
import utils.xtb_utils as xtb
import utils.miscall as misc

kcal_to_eV = 0.0433641153
kB = 8.6173303e-5  # eV/K
T = 298.15
kBT = kB * T
AToBohr = 1.889725989
HToeV = 27.211399

def dft_calc(settings, coords, elements, charge, opt=False, grad=False, hess=False, freeze=[], dirname = None, partial_chrg = False, unp_el = 1, dispersion= False, h20=False):
    print("The dft settings in dft_calc are:", settings)

    if opt and grad:
        exit("opt and grad are exclusive")
    if hess and grad:
        exit("hess and grad are exclusive")

    if hess or grad:
        if len(freeze)!=0:
            print("WARNING: please test the combination of hess/grad and freeze carefully")

    if dirname is None:
        #might add an ordinal number passed somewhere
        rundir=f"dft_tmpdir_{settings['turbomole_functional']}_{settings['turbomole_basis']}_{uuid.uuid4().hex}_"
    else:
        rundir = dirname
    
    if not os.path.exists(rundir):
        os.makedirs(rundir)
    else:
        if len(os.listdir(rundir))>0:
            os.system("rm %s/*"%(rundir)) #removes all files in rundir

    startdir = os.getcwd() #stores current directory
    os.chdir(rundir) #goes to the rundir
    
    print("Starting to prepare input")
    PrepTMInputNormal(".", coords, elements)
    
    #if unp_el != None and unp_el != 0:
    # run calculation
    print("Starting to run TM calculation")
    RunTMCalculation(".", settings, charge, uhf = unp_el, disp = dispersion, pop = partial_chrg, water = h20)
    #else:
    #    RunTMCalculation(".", dft_settings, disp = dispersion, pop = partial_chrg)
    
    # read out results  
    print("TM calculation converged, reading out the results")  
    if opt:
        os.system("t2x coord > opt.xyz")
        coords_new, elements_new = xtb.readXYZ("opt.xyz")
    else:
        coords_new, elements_new = coords, elements

    if grad:
        grad = read_dft_grad()
    else:
        grad = None

    if hess:
        hess, vibspectrum, reduced_masses = read_dft_hess()
    else:
        hess, vibspectrum, reduced_masses = None, None, None

    e = getTMEnergies(".")[2]

    if partial_chrg:
        partialcharges = getMullikans(outfilename = 'TM.out', noOfAtoms=len(coords))
    else:
        partialcharges = None
    # read mull
    
    os.chdir(startdir)

    path_to_control = os.path.join(rundir, "control")
    misc.check_basis_and_func(path_to_control=path_to_control, basis_todo=settings["turbomole_basis"],
                         func_todo=settings["turbomole_functional"])

    #os.system("rm -r %s"%(rundir))

    if settings["delete_calculation_dirs"]:
        os.system("rm -r %s" % (rundir))

    results = {"energy": e, "coords": coords_new, "elements": elements_new, "gradient": grad, "hessian": hess, "vibspectrum": vibspectrum, "reduced_masses": reduced_masses, 'partial_charges': partialcharges}
    print("dft_calc results are: ", results)

    return(results)


def RunTMCalculation(moldir, settings, charge, uhf = None, disp=False, pop = False, water = False):
    startdir = os.getcwd()
    os.chdir(moldir)
    
    #create define string
    if uhf == None or uhf == 1:
        instring = prep_define_file_uhf_1(settings, charge)
        
    if uhf == 3:
        instring = prep_define_file_uhf_3(settings, charge)
    
    print("Starting to execute define string")
    ExecuteDefineString(instring)
    
    # add functional to control file
    func = settings['turbomole_functional']
    add_functional_to_control('control', func)
    
    # add other options to control file like dispersion, solution in water
    if disp:
        AddStatementToControl("control", "$disp3")
    if water:
        AddStatementToControl('control', '$cosmo')
        AddStatementToControl('control', '   epsilon=78.3')  ## adapt?
    if pop:
        AddStatementToControl('control', '$pop')
    
    if settings["copy_mos"]:
        if os.path.exists("%s/pre_optimization/mos"%(settings["main_directory"])):
            print("   ---   Copy the old mos file from precalculation")
            os.system("cp %s/pre_optimization/mos ."%(settings["main_directory"]))
        else:
            print("WARNING: Did not find old mos file in %s/pre_optimization"%(settings["main_directory"]))
    
    # do calculation  
    print("Got to the terminal interaction part of TM!")   
    if settings["turbomole_method"]=="ridft":
        os.system("ridft > TM.out")
        #os.system("rdgrad > rdgrad.out")    ###removed grad
    elif settings["turbomole_method"]=="dscf":
        os.system("dscf > TM.out")
        #os.system("rdgrad > rdgrad.out")
    else:
        exit("ERROR in turbomole_method: %s"%(settings["turbomole_method"]))

    finished=False   
    number_of_iterations=None
    for line in open("TM.out","r"):
        if "convergence criteria satisfied after" in line:
            number_of_iterations=int(line.split()[4])
        if "all done" in line:
            finished=True
            break
    if number_of_iterations!=None:
        print("   ---   converged after %i iterations" %(number_of_iterations))
    else:
        pass

    if finished:
        os.system("eiger > eiger.out")

    os.chdir(startdir)
    return(finished)

#------------------------------------------------------------------preparation functions

def PrepTMInputNormal(moldir, coords, elements):
    coordfile = open("%s/coord"%(moldir), 'w')
    coordfile.write("$coord \n")
    for idx, atom in enumerate(coords):
        coordfile.write("%f  %f  %f  %s\n" % (atom[0] * AToBohr, atom[1] * AToBohr, atom[2] * AToBohr, elements[idx]))
    coordfile.write("$end\n")
    coordfile.close()
    return ()

# define file preperation
def prep_define_file_uhf_1(dft_settings, charge):

    basisset = dft_settings["turbomole_basis"]
    functional = dft_settings["turbomole_functional"]

    try:
        from StringIO import StringIO as mStringIO
    except ImportError:
        from io import StringIO as mStringIO

    outfile = mStringIO()
    
    outfile.write("\n\n\n\n\na coord\n*\nno\n")
    outfile.write("\n\nb all %s\n\n\n" % basisset)
    #if charge !=0: #== +1 or charge == -1
    outfile.write("*\neht\n\n%i\nn\nu 1 \n" % int(charge)) 
    outfile.write("*\nn\n\n\n") 
    #outfile.write("*\neht\n\n%i\ny\nn\n\n\n\n\n\n\n" % int(charge)) ## hier uhf commands \nu int write out natural orbitals? n oder y?
        #outfile.write("\n\n\n\n") ?
    #else:
    #   outfile.write("*\neht\n\n\n\n\n\n\n") # uhf change oder nur 0
    
    if functional != "HF":
        #outfile.write("dft\non\nfunc %s\n\n\n" % functional)
        outfile.write("dft\non\n\n\n")

    if dft_settings["turbomole_method"]=="ridft":
        outfile.write("ri\non\nm1500\n\n\n\n" )

    outfile.write("scf\niter\n100\n\n\n") 
    #outfile.write("scf\ndsp\non\n\n\n") 
    outfile.write("scf\nconv\n5\n\n\n") 
    outfile.write("\n\n\n\n*\n")

    returnstring = outfile.getvalue()
    outfile.close()

    print("Define file is: ", returnstring)
    return returnstring

def prep_define_file_uhf_3(dft_settings, charge):

    basisset = dft_settings["turbomole_basis"]
    functional = dft_settings["turbomole_functional"]

    try:
        from StringIO import StringIO as mStringIO
    except ImportError:
        from io import StringIO as mStringIO

    outfile = mStringIO()
    outfile.write("\n\n\n\n\na coord\n*\nno\n")
    outfile.write("\n\nb all %s\n\n\n" % basisset)
    #if charge !=0: #== +1 or charge == -1
    outfile.write("*\neht\n\n%i\nn\nu 3 \n" % int(charge)) ## hier uhf commands \nu int write out natural orbitals? n oder y?
    outfile.write("*\nn\n\n\n") 
    #else:
    #    outfile.write("*\neht\n\n\n\n\n\n\n") # uhf change oder nur 0
    
    if functional != "HF":
        #outfile.write("dft\non\nfunc %s\n\n\n" % functional) 
        outfile.write("dft\non\n\n\n")

    if dft_settings["turbomole_method"]=="ridft":
        outfile.write("ri\non\nm1500\n\n\n\n" )

    outfile.write("scf\niter\n500\n\n\n") #1000
    #outfile.write("scf\ndsp\non\n\n\n") ## statement nicht erkannt ? 
    outfile.write("scf\nconv\n5\n\n\n") # 7?
    outfile.write("\n\n\n\n*\n")

    returnstring = outfile.getvalue()
    outfile.close()
    return returnstring

def add_functional_to_control(file_path, func):
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()

        with open(file_path, 'w') as file:
            for line in lines:
                if 'functional' in line and 'b-p' in line:
                    # Replace 'b-p' with the value of 'func'
                    line = line.replace('b-p', func)
                file.write(line)

        #print(f"Text 'b-p' replaced with '{func}' in '{file_path}'")

    except FileNotFoundError:
        print(f"File '{file_path}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

def AddStatementToControl(controlfilename, statement):
    inf = open(controlfilename, 'r')
    lines = inf.readlines()
    inf.close()
    already_in=False
    outf = open(controlfilename, 'w')
    for line in lines:
        if statement.split()[0] in line:
            already_in=True
        if len(line.split()) > 0:
            if line.split()[0] == "$end" and not already_in:
                outf.write("%s\n" % (statement))
        outf.write(line)
    outf.close()

def RemoveStatementFromControl(controlfilename, statement):
    inf = open(controlfilename, 'r')
    lines = inf.readlines()
    inf.close()
    outf = open(controlfilename, 'w')
    writeOutput = True
    for line in lines:
        if len(line.split()) > 0:
            if line.split()[0] == statement.split()[0]:
                writeOutput = False
            else:
                writeOutput = True
        if writeOutput:
            outf.write(line)
    outf.close()

def ExecuteDefineString(instring):
    instring = instring + "\n\n\n\n"
    out = ""
    err = ""

    process = subprocess.Popen(["define"], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, preexec_fn=setulimit, encoding='utf8')
    out, err = process.communicate(input = instring)

    if "normally" in err.split():
        return
    if "normally" not in err.split():
        print("ERROR in define")
        print("STDOUT was: %s"%(out))
        print("STDERR was: %s"%(err))
        print("Now printing define input:")
        print("--------------------------")
        print(instring)
        print("--------------------------")
        with open("define.input",'w') as defineinput:
            defineinput.write(instring)
        exit()

def setulimit():
    resource.setrlimit(resource.RLIMIT_STACK,(-1,resource.RLIM_INFINITY))

#----------------------------------------------------read out calculation results
    
def getTMEnergies(moldir):

    print("current directory at the start of getTMEnergies function is ", os.getcwd())
    print("moldir is defined as: ", moldir)
    try:
        eigerfile=open("%s/eiger.out"%(moldir),"r")
    except FileNotFoundError as err:
        raise FileNotFoundError("the eiger.out file is not found, was searching in directory", moldir)

    eigerlines=eigerfile.readlines()
    eigerfile.close()
    total_energy=0.0
    energy_homo=0.0
    energy_lumo=0.0
    for eigerline in eigerlines:
        if len(eigerline.split())!=0:
            if eigerline.split()[0]=="Total":
                total_energy=eigerline.split()[6]
            elif eigerline.split()[0]=="HOMO:":
                energy_homo=eigerline.split()[8]
            elif eigerline.split()[0]=="LUMO:":
                energy_lumo=eigerline.split()[8]
                break
    return([float(energy_homo),float(energy_lumo),float(total_energy)])

def read_dft_grad():
    if not os.path.exists("gradient"):
        return(None)
    grad = []
    for line in open("gradient","r"):
        if len(line.split())==3 and "grad" not in line:
            line = line.replace("D","E")
            grad.append([float(line.split()[0]), float(line.split()[1]), float(line.split()[2])])
    if len(grad)==0:
        grad=None
    else:
        grad = np.array(grad)*HToeV*AToBohr
    return(grad)

def read_dft_hess():
    hess = None
    if not os.path.exists("hessian"):
        return(None, None, None)
    hess = []
    for line in open("hessian","r"):
        if "hess" not in line:
            for x in line.split():
                hess.append(float(x))
    if len(hess)==0:
        hess=None
    else:
        hess = np.array(hess)

    vibspectrum = None
    if not os.path.exists("vibspectrum"):
        return(None, None, None)
    vibspectrum = []
    read=False
    for line in open("vibspectrum","r"):
        if "end" in line:
            read=False

        if read:
            if len(line.split())==5:
                vibspectrum.append(float(line.split()[1]))
            elif len(line.split())==6:
                vibspectrum.append(float(line.split()[2]))
            else:
                print("WARNING: weird line length: %s"%(line))
        if "RAMAN" in line:
            read=True
    
    reduced_masses = None
    if not os.path.exists("g98.out"):
        print("g98.out not found")
        return(None, None, None)
    reduced_masses = []
    read=False
    for line in open("g98.out","r"):
        if "Red. masses" in line:
            for x in line.split()[3:]:
                try:
                    reduced_masses.append(float(x))
                except:
                    pass

    if len(vibspectrum)==0:
        vibspectrum=None
        print("no vibspectrum found")
    else:
        vibspectrum = np.array(vibspectrum)

    if len(reduced_masses)==0:
        reduced_masses = None
        print("no reduced masses found")
    else:
        reduced_masses = np.array(reduced_masses)

    return(hess, vibspectrum, reduced_masses)


def getMullikans(outfilename="ridft.log", noOfAtoms=0):
    TMfile=open(outfilename,"r")
    lines=TMfile.readlines()
    TMfile.close()
    # finding position of mullikan partial charges in file
    idx=0
    for line in lines:
        spl=line.split()
        if len(spl)!=0 and len(spl)!=1:
            if spl[0]=="atom" and spl[1]=="charge":
                index=idx+1
                break
        idx+=1
    partialCharges=[]
    for i in range(noOfAtoms):
        partialCharges.append( float(lines[idx+1+i][10:18]) )
    return(partialCharges)