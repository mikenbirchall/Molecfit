#!/usr/bin/env python3


import os
import shutil
import sys
import numpy as np
import math

# GLOBAL VARIABLES:
LBLRTM_BIN     ="lblrtm"
LBLRTM_MRG_BIN ="ODintMerge"
OD_MERGE_FILE  ="OD_merged.dat"
TOPL_DIR=os.getcwd()
TAPE5_FILENAME='TAPE5'
TAPE1_FILENAME='TAPE1'
TAPE3_FILENAME='TAPE3'
TAPE6_FILENAME='TAPE6'
TAPE7_FILENAME='TAPE7'
TAPE10_FILENAME='TAPE10'
TAPE3_DIRNAME="TAPE3_DIR"

# -----------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------

# ==========
# Functions:
# ==========


# -----------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------

# -------------
# MAIN ROUTINE:
# -------------
print(sys.argv[0])

# Get the list of moleculrs in the TAPE3 directory
folder=os.path.join(TOPL_DIR,TAPE3_DIRNAME)
tape5 =os.path.join(TOPL_DIR,TAPE5_FILENAME)
mol_lst = [ f.name for f in os.scandir(folder) if f.is_dir() ]
print("Molecules of Interest = ", mol_lst);

# Iterate through each molecule
od_file_lst=[]
for mol in mol_lst:

    # Make clean directory for this molecule
    if (os.path.exists(mol)):
            print("exists")
            shutil.rmtree(mol)
    os.mkdir(mol)

    # Soft link the associated TAPE3 file into this directory
    source=os.path.join(folder,mol,TAPE3_FILENAME)
    target=os.path.join(mol,TAPE3_FILENAME)
    os.symlink(source,target)
    wdir=os.path.join("./",mol)
    target=os.path.join(wdir,TAPE5_FILENAME)
    print("Work in path ", wdir)

    #  Create a softlink to the LBLRTM TAPE5 into this directory
    if (os.path.exists(target)):
            os.unlink(target)
    os.symlink(tape5,target)

    # Invoke lblrtm in this working directory
    cmd_str="cd " + wdir + " ; " + LBLRTM_BIN
    print(cmd_str)
    os.system(cmd_str)

    # Invoke the lblrtio merge application
    cmd_str="cd " + wdir + " ; " + LBLRTM_MRG_BIN
    print(cmd_str)
    #os.system(cmd_str)

    od_file=os.path.join(wdir,OD_MERGE_FILE)
    od_file_lst.append(od_file)
sys.exit()
cnt=0
for od_file in od_file_lst:
    print ("Merging ", od_file)

    # Read in file data
    fid=open(od_file,'r')
    data=fid.read()
    lines=data.split("\n")
    fid.close()

    # If this is the first file then allocate od vector of appropriate size
    if (cnt==0):
        n=len(lines)-1
        odA=np.zeros(n)
        fqA=np.zeros(n)
    cnt=cnt+1

    #Iterate through the lines
    line_idx=0
    for line in lines:
        lst=line.split()
        if (line_idx<5 or line_idx > n-5) :
            print (line_idx+1, "line=",line, " list=", lst)
        fqval=lst[0]
        odval=lst[1]
        odA[line_idx]=odA[line_idx]+float(odval)
        if (cnt==1):
            fqA[line_idx]=fqval
        line_idx=line_idx+1
        if (line_idx==n):
            break

# Now dump the od vector
fid=open(OD_MERGE_FILE,'w')
for i in range(0,n-1):
    fid.write( str( fqA[i] ) + " " + str( math.exp(-1.0*odA[i] ) ) + "\n")
fid.close()


