
# WIP RELATED HACKS
TOPL_DIR=os.getcwd()
ALL_FLAG=True



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


def ForceMakeDir(dir):

    # Forecably make a directory, ie if it already exists then delete first.

    if os.path.exists(dir):
            shutil.rmtree(dir)
    os.mkdir(dir)

# -----------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------

# -------------
# MAIN ROUTINE:
# -------------
# Clean out any leftover Molecular Specific Directories
#CleanOldMolecularDirectories(MOLECULES_FULL_LST)
print("lnflOda ODA_OPTION=", ODA_OPTION)

# Check if option is to run the default method
if (ODA_OPTION=="NONE" or ODA_OPTION=="STD" or ODA_OPTION=="BOTH"):
    cmd_str=LNFL_BIN
    print("Running Standard:", cmd_str)
    os.system(cmd_str)
#    shutil.copyfile(TAPE3_FILENAME,'TAPE3_BACKUP')

# Check if option is to run the Oda method
if (ODA_OPTION!="ODA" and ODA_OPTION!="BOTH"):
    # No request to run ODA so exist
    sys.exit()


# If the output directory already exists then remove it
if (os.path.exists(TAPE3_DIRNAME)):
    shutil.rmtree(TAPE3_DIRNAME)

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
#print ("TAPE5 RECORD 3          =",record3)
#print ("MOLECULE SELECTION FIELD=",field)
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

    ForceMakeDir(mol_name)
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
    ForceMakeDir(ALLMOLS_DIRNAME)
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



