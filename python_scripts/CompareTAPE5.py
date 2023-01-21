import numpy as np
import os

MAX_NMOLS=50

def READ_TAPE5(dirname):

    verbose=False
    filename=os.path.join(dirname,'TAPE5')

    # Read TAPE5 file and dump int a data list
    fid=open(filename,'r')
    data=fid.readlines()
    fid.close()


    # Iterate through each line and parse
    cnt=0

    max_nmols=0
    record_section=False
    new_record    =False
    for line in data:

        cnt=cnt+1
        if (cnt<10 and verbose) :
            print(cnt, ") ", line)

        if (record_section and new_record):

            # Parse the line as the header record in the form "val val val AA AAAA??"
            lst=line.split()
            hgth=lst[0]
            pres=lst[1]
            temp=lst[2]
            mols=lst[4]         # A string of 'A' characters one for each molecules
            nmols=len(mols)     # No of molecules = no of 'A' characters
            max_nmols=max(max_nmols,nmols)
            nmol_lines=(nmols-1)//8+1 # Record is broken up into groups of 8 molecules
            hvec[layer]=hgth
            pvec[layer]=pres
            tvec[layer]=temp
            new_record=False
            record_part=0
            if (verbose):
                print (layer+1, hgth, pres, temp, mols,nmol_lines)

        elif (record_section and layer<nlayers) :
            record_part=record_part+1
            # Parse the line as a list of abundencies
            lst=line.split()
            for i in range(len(lst)):
                j=i
                abun[layer,j]=lst[i]
            if (record_part==nmol_lines):
                new_record=True
                layer=layer+1
            if (layer==nlayers):
                new_record=False
                record_section=False

        if (cnt==8) :
            # This is a datum line containing the no of layers
            nlayers=int(line)
            # Now allocate the data vactors and arrays
            abun=np.zeros((nlayers,MAX_NMOLS))
            hvec=np.zeros((nlayers))
            tvec=np.zeros((nlayers))
            pvec=np.zeros((nlayers))

            # Decalre we are now in the record section and we are reading a new record
            new_record    =True
            record_section=True
            layer=0     # Start the layer index at 0

    if (verbose) :
        print('N lines=',cnt)
        print('n_layers=',nlayers)

    return hvec,pvec,tvec,abun,nlayers,max_nmols

# end def READ_TAPE5(dirname)

def printData(hvec,pvec,tvec,abunA,nlayers,nmols):
    print("Number of Layers=",nlayers)
    for layer in range(nlayers):
        lst=[]
        for mol_idx in range(nmols):
            lst.append(abunA[layer,mol_idx])
        print (layer+1,hvec[layer],pvec[layer],tvec[layer],lst)

#end def printData(hvec,pvec,tvec,abunA,nlayers,nmols):

def CompareTAPE5(dirname1,dirname2):

    hvec1,pvec1,tvec1,abunA1,nlayers1,nmols1=READ_TAPE5(dirname1)
    hvec2,pvec2,tvec2,abunA2,nlayers2,nmols2=READ_TAPE5(dirname2)
    different_nlayers=False
    different_nmols  =False
    if (nlayers1!=nlayers2):
        different_nlayers=True
    if (nmols1!=nmols2):
        different_nmols  =True
    print("Check layers: ",nlayers1,nlayers2, "Diff=", different_nlayers)
    print("Check nmols: ", nmols1,nmols2, "Diff=", different_nmols)
    if (different_nlayers or different_nmols):
        return
    nlayers=nlayers1
    nmols=nmols1
    sA=np.zeros((nlayers,nmols))
    for layer in range(nlayers):
        lst=[]
        for mol_idx in range(nmols):
            a1=abunA1[layer][mol_idx]
            a2=abunA2[layer][mol_idx]
            if (a2==0.0):
                scale=0.0
            else:
                scale=a1/a2
            lst.append(scale)
            sA[layer][mol_idx]=scale
        print (layer,lst)
    lst=[]
    avgV=np.median(sA,axis=0)
    for mol_idx in range(nmols):
        avg=avgV[mol_idx]
        lst.append(avg)
    return lst

# end def CompareTAPE5(dirname1,dirname2):
if __name__=="__main__":
    hvec,pvec,tvec,abunA,nlayers,nmols=READ_TAPE5(".")

#printData(hvec,pvec,tvec,abunA,nlayers,nmols)

    wrkdir="/home/mnb/workspace/ESO/molecfit/work_dir/crires1"
    tmpdir="telluriccorr_tmp_folder_LJ8w1z"
    run1_1 ="run_1_wdir_lblrtm_range_1_lblrtm_call_1"
    run1_2 ="run_6_wdir_lblrtm_range_1_lblrtm_call_91"
    dir1=os.path.join(wrkdir,tmpdir,run1_1,"wv_number_0")
    dir2=os.path.join(wrkdir,tmpdir,run1_2,"wv_number_0")
    lst=CompareTAPE5(dir1,dir2)
    print("AVG:",lst)
    print(dir1)
    print (dir2)
