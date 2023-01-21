#!/usr/bin/env python3

import os
import shutil
import sys
import tarfile

# WIP RELATED HACKS
LNFL_BIN="lnfl"
TOPL_DIR=os.getcwd()
ALL_FLAG=True


# GLOBAL VARIABLES:
TAPE5_FILENAME  = 'TAPE5'
TAPE1_FILENAME  = 'TAPE1'
TAPE3_FILENAME  = 'TAPE3'
TAPE6_FILENAME  = 'TAPE6'
TAPE7_FILENAME  = 'TAPE7'
TAPE10_FILENAME = 'TAPE10'
TAPE3_DIRNAME   = "TAPE3_DIR"
ALLMOLS_DIRNAME = "ALL"
N_MOL_TYPES=47
MOL_STR="\
( 1)  H2O  ( 2)  CO2  ( 3)    O3 ( 4)   N2O ( 5)    CO ( 6)   CH4 ( 7)    O2 \
( 8)   NO  ( 9)  SO2  (10)   NO2 (11)   NH3 (12)  HNO3 (13)    OH (14)    HF \
(15)  HCL  (16)  HBR  (17)    HI (18)   CLO (19)   OCS (20)  H2CO (21)  HOCL \
(22)   N2  (23)  HCN  (24) CH3CL (25)  H2O2 (26)  C2H2 (27)  C2H6 (28)   PH3 \
(29) COF2  (30)  SF6  (31)   H2S (32) HCOOH (33)   HO2 (34)     O (35)CLONO2 \
(36)  NO+  (37) HOBR  (38)  C2H4 (39) C3HOH (40) CH3Br (41) CH3CN (42)   CF4 \
(43) C4H2  (44) HC3N  (45)    H2 (46)    CS (47)   SO3                       \
"
# Parse the MOLECULAR set string into a workable list
def ParseMOLSTR(MOL_STR):

    ret_lst=[]
    for i in range(1,N_MOL_TYPES+1):
        if (i<10):
            substr="( "+str(i)+")"
        else:
            substr="("+str(i)+")"
        idx =MOL_STR.find(substr)
        mols=MOL_STR[idx+4:idx+11].strip()
        ret_lst.append(mols)
    # Add an all moles dirname to the list
    ret_lst.append(ALLMOLS_DIRNAME)
    return ret_lst
MOLECULES_FULL_LST=ParseMOLSTR(MOL_STR)
# Store the ALLMLS idx (which will be at the end of the list)
ALLMOLS_IDX=len(MOLECULES_FULL_LST)-1

# -----------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------

# ==========
# Functions:
# ==========

def CleanOldMolecularDirectories(dir_names):

    if (os.path.exists(TAPE3_FILENAME)):
            os.remove(TAPE3_FILENAME)

    filename_lst=[TAPE1_FILENAME,TAPE5_FILENAME,TAPE10_FILENAME,TAPE3_FILENAME,TAPE6_FILENAME,TAPE7_FILENAME]

    for dirname in dir_names:
        if (os.path.exists(dirname)):
            print("Directory ", dirname, " exists")
            for filename in filename_lst:
                pathname=os.path.join(dirname,filename)
                if (os.path.exists(pathname)):
                    print("File ", pathname, " exists")
                    os.remove(pathname)
                if (os.path.islink(pathname)):
                    print("File ", pathname, " exists as link")
                    os.unlink(pathname)

            os.rmdir(dirname)


# -----------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------

# -------------
# MAIN ROUTINE:
# -------------
print(sys.argv[0])

# If the output directory already exists then remove it
if (os.path.exists(TAPE3_DIRNAME)):
    shutil.rmtree(TAPE3_DIRNAME)

# Clean out any leftover Molecular Specific Directories
CleanOldMolecularDirectories(MOLECULES_FULL_LST)

# Open TAPE 5 and read the 5 records that it contains
fid=open(TAPE5_FILENAME,'r')
record1=fid.readline()
record2=fid.readline()
record3=fid.readline()
record4=fid.readline()
record5=fid.readline()
fid.close()

# Parse record3 of the TAPE5 file into a list of molecules and define
# sel_lst as the list ofindicies of the molecules flagged in record3
lst=record3.split()
field=lst[0]
print ("TAPE5 RECORD 3          =",record3)
print ("MOLECULE SELECTION FIELD=",field)
sel_lst=[]
for idx in range(0,N_MOL_TYPES):
    flag=field[idx]
    if (flag=="1"):
        sel_lst.append(idx)
        print("Selected:", MOLECULES_FULL_LST[idx])


# Create a subdirectory for each selected molecule and populate with
# a copy of TAPE5 but with record3 only having flag for that molecule
for idx in sel_lst:
    mol_name=MOLECULES_FULL_LST[idx]
    field=""
    for idx2 in range(0,N_MOL_TYPES):
        if (idx2==idx):
            field=field+"1"
        else:
            field=field+"0"
    new_record=field+record3[N_MOL_TYPES:len(record3)]

    os.mkdir(mol_name)
    newfile=os.path.join(mol_name,TAPE5_FILENAME)
    fid=open(newfile,'w')
    fid.write(record1)
    fid.write(record2)
    fid.write(new_record)
    fid.write(record4)
    fid.write(record5)
    fid.close()

    # Now make a softlink for TAPE1
    newfile=os.path.join(mol_name,TAPE1_FILENAME)
    srcfile=os.path.join("../",TAPE1_FILENAME)
    os.symlink(srcfile,newfile)

# If ALL flag is true then create a "Molecule" subdirectoy
# for all mols, which will conatin a Softlink TAPE5 to the
# full TAPE5
if ALL_FLAG:
    # Make directory ALL
    os.mkdir(ALLMOLS_DIRNAME)
    # Softlink TAPE5
    newfile=os.path.join(ALLMOLS_DIRNAME,TAPE5_FILENAME)
    srcfile=os.path.join("../",TAPE5_FILENAME)
    os.symlink(srcfile,newfile)
    #Softlink TAPE2
    newfile=os.path.join(ALLMOLS_DIRNAME,TAPE1_FILENAME)
    srcfile=os.path.join("../",TAPE1_FILENAME)
    os.symlink(srcfile,newfile)
    # Now add the ALL molecule idx to the selection list
    sel_lst.append(ALLMOLS_IDX)

# Now execute lnfl in each molecule directory
for idx in sel_lst:
    mol_name=MOLECULES_FULL_LST[idx]
    cmd_str="cd " + mol_name + " ; " + LNFL_BIN
    print(cmd_str)
    os.system(cmd_str)

# Now Create a Directory Substructure suitable for lblrtm
os.mkdir(TAPE3_DIRNAME)
for idx in sel_lst:
    mol_name=MOLECULES_FULL_LST[idx]
    pathname=os.path.join(TAPE3_DIRNAME,mol_name)
    os.mkdir(pathname)
    target=os.path.join(pathname,TAPE3_FILENAME)
    print("TOPL_DIR=",TOPL_DIR)
    source=os.path.join(mol_name,TAPE3_FILENAME)
    source=os.path.join(TOPL_DIR,source)
    os.symlink(source,target)



