#!/usr/bin/env python3

import os
import shutil
import sys
import tarfile
import numpy as np

# GLOBAL VARIABLES:
LBLRTM_BIN      = 'lblrtm'
MOLPRFSDIR      = 'INDIVIDUAL_MOLECULE_PROFS'
TAPE5_FILENAME  = 'TAPE5'
TAPE1_FILENAME  = 'TAPE1'
TAPE3_FILENAME  = 'TAPE3'
TAPE6_FILENAME  = 'TAPE6'
TAPE7_FILENAME  = 'TAPE7'
TAPE10_FILENAME = 'TAPE10'
TAPE28_FILENAME = 'TAPE28'
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

def Mol2Idx(mol_str):
    idx=-1
    n_molecules=len(MOLECULES_FULL_LST)
    for i in range(n_molecules):
        if(MOLECULES_FULL_LST[i]==mol_str):
            idx=i+1
            return idx
    return idx
# ==================================================================================
#                     CAST TAPE5 UTILITY ROUTINES BEGIN
# ==================================================================================
def setzero(str, idx):
    idx1=(idx-1)*10+1
    idx2=idx1+9
    nstr=str[0:idx1] + "0.000e+00" + str[idx2:]
    return nstr

def setzeros(str0, idx0):
    lst=str0.split()
    n=len(lst)
    nstr=str0
    for i in range(n):
        idx=i+1
        if idx!=idx0:
            nstr=setzero(nstr,idx)
    return nstr

def setallzeros(str0):
    lst=str0.split()
    n=len(lst)
    nstr=str0
    for i in range(n):
        idx=i+1
        nstr=setzero(nstr,idx)
    return nstr

def READ_TAPE5(dirname):

    verbose=False
    filename=os.path.join(dirname,'TAPE5')

    # Read TAPE5 file and dump int a data list
    fid=open(filename,'r')
    data=fid.readlines()
    fid.close()
    header=[]
    for idx in range(7):
        header.append(data[idx].strip("\n"))
    prerecords=data[7].strip("\n")
    nlayers=int(data[7])
    print(nlayers)
    records=[]
    line_idx=8
    if (nlayers>0): # There will be at least one sample line folowing
        sample_line=data[8]
        lst=sample_line.split() # Get the multiple AAA string as the 5th
        astr=lst[4]             # element in this line
        nas=len(astr)           # Count the number of A's
        n_lines_per_record=(nas-1)//8+1 # Each record is split into lines of 8 molecules
        n_total_lines=nlayers*(1+n_lines_per_record)
        line_idx=8 + n_total_lines
        for i in range(nlayers):
            beg_idx=8+i*(1+n_lines_per_record)
            end_idx=beg_idx+1+n_lines_per_record
            lst=[]
            for idx in range(beg_idx,end_idx):
                lst.append(data[idx].strip("\n"))
            records.append(lst)
    footer=[]
    for idx in range(line_idx,len(data)):
        footer.append(data[idx].strip("\n"))
    print(astr,nas, n_lines_per_record)
    return header, prerecords, records, footer

def writeTAPE5(filename,header, prerecords, records, footer):
    fid=open(filename,'w')
    for line in header:
        fid.write(line+"\n")
    fid.write(prerecords+"\n")
    for record in records:
        for line in record:
            fid.write(line+"\n")
    for line in footer:
        fid.write(line+"\n")
    fid.close()

def modRecord2SingleMol(record,mol_idx):

    new_record=[]
    new_record.append(record[0])
    n=len(record)
    nparts=n-1
    significant_part=(mol_idx-1)//8+1
    significant_idx = mol_idx-(significant_part-1)*8
    for part_idx in range(1,nparts+1):
        record_part=record[part_idx]
        if (part_idx==significant_part):
            new_part=setzeros(record_part,significant_idx)
        else:
            new_part=setallzeros(record_part)
        new_record.append(new_part)

    return new_record

def modRecords4SingleMol(records,mol_idx):
    new_records=[]
    for record in records:
        new_record=modRecord2SingleMol(record,mol_idx)
        new_records.append(new_record)
    return new_records

def test():
    str=" 1.312e+03 0.000e+00 2.903e-02 0.000e+00 0.000e+00 1.768e+00"
    str1 = setzero(str,1)
    str3 = setzero(str,3)
    str6  = setzero(str,6)

    print(str)
    print(str1)
    print(str3)
    print(str6)
    strx=setzeros(str,3)
    print(strx)

    header, prerecords, records, footer = READ_TAPE5("./")
    print ("HEADER")
    for line in header:
        print (line)
    print ("PRERECORD")

    records1=modRecords4SingleMol(records,1)
    records3=modRecords4SingleMol(records,3)
    records6=modRecords4SingleMol(records,6)
    writeTAPE5("test0.txt",header, prerecords, records, footer)
    writeTAPE5("test1.txt",header, prerecords, records1, footer)
    writeTAPE5("test3.txt",header, prerecords, records3, footer)
    writeTAPE5("test6.txt",header, prerecords, records6, footer)
    print(setallzeros(str))
    rec=records[0]
    print(rec)
    nrec=modRecord2SingleMol(rec,3)
    print(nrec)

def ModifyTAPE5Header4dv(header,dv):

    # We need to express the dv value as a fixed string format
    f = format("%.3e" % dv)
    # Replace the in indecies [31:40] of the 3rd line in the header
    elem=header[2]
    new_elem=elem[0:31]+str(f)+elem[40:]

    # Create modified header as copy of the header with
    # new_elem replacing the 3rd line
    mod_header=[]
    nlines=len(header)
    for i in range(nlines):
        if (i==2):
            mod_header.append(new_elem)
        else :
            mod_header.append(header[i])

    # Return the modified header
    return mod_header

def cast_TAPE5 (TAPE5source,wdir,molecule_list,dv):

    # Parse the TAPE source file into sections
    header, prerecords, records, footer = READ_TAPE5(os.path.dirname(TAPE5source))

    # Modify the header to have the specified dv value
    mod_header=ModifyTAPE5Header4dv(header,dv)

    # For each molecule in given list generate a TAPE5 file with the dv value
    # and specification that only this moecule is to be used
    for molecule in molecule_list:

        # Specify the new TAPE5 file to create in this directory
        target=os.path.join(wdir,molecule,TAPE5_FILENAME)

        # Special case for psuedo molecule 'ALL'
        if (molecule=="ALL"):

            snglmol_records=records

        else:
            # Determine the idx for this molecule
            mol_idx=Mol2Idx(molecule)

            # Modify the records do that only the flas for this molecule
            # is non zero
            snglmol_records=modRecords4SingleMol(records,mol_idx)

        # Write the new TAPE5
        writeTAPE5(target,mod_header, prerecords, snglmol_records,footer)

# ==================================================================================
#                     CAST TAPE5 UTILITY ROUTINES END
# ==================================================================================
# ==================================================================================

def ReadTAPE28(dirname):

    # Read all lines in asci file TAPE28
    file='TAPE28'
    filename=os.path.join(dirname,file)

    # Check if present if not assume it has been moved to ../TAPE28_1
    if (not os.path.exists(filename)):
        target=os.path.join(dirname,'../TAPE28_1')
        os.symlink(target,filename)

    fid=open(filename,'r')
    lines=fid.readlines()
    fid.close()

    # Iterate through all lines and extract Wavenumbers and Transmission values
    cnt=0;            # No of lines with values extracted
    extracting=False; # Boolean if in extracting mode
    wave_lst=[];      # The wavenumber extracted values list
    trans_lst=[]       # The transmission extracted values list
    for line in lines:

        # Skip empty lines
        if (len(line.strip())==0):
            continue

        # If this line contains "WAVENUMBER" and "TRANSMISSION" then start extracting
        # on the next iteration
        if (line.find("WAVENUMBER")!=-1 and line.find("TRANSMISSION")!=-1):
                extracting=True
                continue

        # Skip if we are not extracting
        if (not extracting):
            continue

        # Real Values are to be extracted from this line
        cnt+=1
        lst=line.split()
        wave_lst.append (lst[0])
        trans_lst.append(lst[1])

    # Convert to numpy vectors
    waveV=np.zeros(cnt)
    transV=np.zeros(cnt)
    for i in range (cnt):
        waveV [i] = wave_lst[i]
        transV[i] =trans_lst[i]

    return waveV, transV, cnt

def GetdvFromTAPE28(dirname):

    waveV, transV, cnt = ReadTAPE28(dirname)
    stepV=np.zeros(cnt-1)
    for i in range(cnt-1):
        stepV[i]=waveV[i+1]-waveV[i]
    dv=np.median(stepV)
    return dv

def GetTAPE3_DIRPath(dirname):
    tape3=os.path.join(dirname,TAPE3_FILENAME)
    realpath=os.path.realpath(tape3)
    realdir=os.path.dirname(realpath)
    tape3dir=os.path.join(realdir,TAPE3_DIRNAME)
    return tape3dir

def MakeMOLPRFSDir(dirname,tape5,tape28,tape3dir):
    molprfs_dir=os.path.join(dirname,MOLPRFSDIR)
    print ("HERE dirname=",dirname, "molprfs=",molprfs_dir)
    if (os.path.exists(molprfs_dir)):
        print("Molecule Profile Directory allready exists. Will delete.")
        shutil.rmtree(molprfs_dir)
    os.mkdir(molprfs_dir)
    target=os.path.join(molprfs_dir,TAPE5_FILENAME)
    shutil.copy(tape5,target)
    target=os.path.join(molprfs_dir,TAPE28_FILENAME)
    shutil.copy(tape28,target)
    target=os.path.join(molprfs_dir,TAPE3_DIRNAME)
    os.symlink(tape3dir,target)
    return molprfs_dir


def CastTAPE3(dirname):
    # Get the list of moleculrs in the TAPE3 directory
    folder=os.path.join(dirname,TAPE3_DIRNAME)
    tape5 =os.path.join(dirname,TAPE5_FILENAME)
    mol_lst = [ f.name for f in os.scandir(folder) if f.is_dir() ]
    print("Molecules of Interest = ", mol_lst);

    # Iterate through each molecule
    od_file_lst=[]
    for mol in mol_lst:
        moldir=os.path.join(dirname,mol)
        # Make clean directory for this molecule
        if (os.path.exists(moldir)):
            print("exists")
            shutil.rmtree(mol)
        os.mkdir(moldir)

        # Soft link the associated TAPE3 file into this directory
        source=os.path.join(folder,mol,TAPE3_FILENAME)
        source=os.path.realpath(source)
        target=os.path.join(moldir,TAPE3_FILENAME)
        os.symlink(source,target)

def CastTAPE5(dirname,dv):
    # Get the list of molecules in the TAPE3 directory
    folder=os.path.join(dirname,TAPE3_DIRNAME)
    tape5 =os.path.join(dirname,TAPE5_FILENAME)
    mol_lst = [ f.name for f in os.scandir(folder) if f.is_dir() ]
    print("Molecules of Interest = ", mol_lst);

    cast_TAPE5(tape5,dirname,mol_lst,dv)
    return

    # Iterate through each molecule
    od_file_lst=[]
    for mol in mol_lst:
        moldir=os.path.join(dirname,mol)

        # Soft link the associated TAPE3 file into this directory
        source=os.path.join(dirname,TAPE5_FILENAME)
        target=os.path.join(moldir,TAPE5_FILENAME)
        shutil.copy(source,target)

def RunLBLRTM(dirname):
    # Get the list of molecules in the TAPE3 directory
    folder=os.path.join(dirname,TAPE3_DIRNAME)
    mol_lst = [ f.name for f in os.scandir(folder) if f.is_dir() ]
    for mol in mol_lst:
        wdir=os.path.join(dirname,mol)
        cmd_str="cd " + wdir + " ; " + LBLRTM_BIN
        print(cmd_str)
        os.system(cmd_str)


print("Molecule by Molecule Profile Extraction")
dv= GetdvFromTAPE28(".")
tape3dir_path=GetTAPE3_DIRPath(".")
outdir=os.path.dirname(tape3dir_path)
print ("dv estimates at", dv)
print ("TAPE3_DIR found at ", tape3dir_path)
molprfs_dir=MakeMOLPRFSDir(outdir,TAPE5_FILENAME,TAPE28_FILENAME,tape3dir_path)
CastTAPE3(molprfs_dir)
CastTAPE5(molprfs_dir,dv)
RunLBLRTM(molprfs_dir)