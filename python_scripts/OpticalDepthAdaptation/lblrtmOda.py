# GLOBAL VARIABLES:
LBLRTM_MRG_BIN ="ODintMerge"
OD_MERGE_FILE  ="OD_merged.dat"
TOPL_DIR=os.getcwd()

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
print("lblrtmOda ODA_OPTION=", ODA_OPTION)

# Check if option is to run the default method
if (ODA_OPTION=="NONE" or ODA_OPTION=="STD" or ODA_OPTION=="BOTH"):
    cmd_str=LBLRTM_BIN
    print(cmd_str)
    os.system(cmd_str)
    print("Finished here")


# Check if option is to run the Oda method. If not then exit
if (ODA_OPTION!="ODA" and ODA_OPTION!="BOTH"):
    # No request to run ODA so exist
    sys.exit()

# Check if the Molecular Profile Exists and if not then  create it with mbmOda
if (not MolPrfDirExists(".")):
    cmd_str="mbmOda"
    print(cmd_str)
    os.system(cmd_str)

# Create the Oda version of TAPE28 via mrgmOda
cmd_str="mrgOda"
print(cmd_str)
os.system(cmd_str)
