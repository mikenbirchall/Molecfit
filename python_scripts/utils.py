import os

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

def cast_TAPE5 (TAPE5source,wdir,molecule_list):

    header, prerecords, records, footer = READ_TAPE5(TAPE5source)

    for molecule in molecule_list:
        # WIP CHEATS HERE
        if (molecule=="H2O"):
            mol_idx=1
        if (molecule=="O3"):
            mol_idx=3
        if (molecule=="CH4"):
            mol_idx=6
        target=os.path.join(wdir,molecule,"TAPE5")
        snglmol_records=modRecords4SingleMol(records,mol_idx)
        writeTAPE5(target,header, prerecords, snglmol_records,footer)

path="/home/mnb/workspace/ESO/molecfit/work_dir/crires1/telluriccorr_tmp_folder_LJ8w1z/run_1_wdir_lblrtm_range_1_lblrtm_call_1/wv_number_0"
molecules=["H2O","O3","CH4"]
wdir=path
cast_TAPE5(wdir,wdir,molecules)
