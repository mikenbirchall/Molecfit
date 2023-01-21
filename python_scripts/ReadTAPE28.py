import itertools
import os
import numpy as np
import math
import CompareTAPE5

def GenCombTrans(mol_lst,r_lst,npts):
    n=len(r_lst)
    print("Combining ", n, " molecule transitivity profiles")
    tV=np.zeros(npts)
    for idx in range(n):
        sV=mol_lst[idx]
        r =r_lst[idx]
        #r=1.0/r
        for j in range(npts):
            trans=max(sV[j],1.0e-08)
            tau=-1.0*math.log(trans)
            tV[j]=tV[j]+r*tau

    for j in range(npts):
        tau=tV[j]
        tV[j]=math.exp(-tau)

    return tV


def dump4plot(fname,xV,yV,n):
    fid=open(fname,'w')
    for i in range(n):
        lst = []
        lst.append(xV[i])
        lst.append(yV[i])
        #line=' '.join(lst)
        #fid.write(line)
        fid.write(str(xV[i])+' '+str(yV[i])+'\n')
        #fid.write(str(yV[i]))
        #fid.write('\n')
        #fid.write("%d %d \n" % xV[i] % yV[i])
    fid.close()



def ReadTAPE28(dirname):

    # Read all lines in asci file TAPE28
    file='TAPE28'
    filename=os.path.join(dirname,file)
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


def AnalyseData(waveV,transV,n):
    stepV=np.zeros(n-1)
    for i in range(n-1):
        stepV[i]=waveV[i+1]-waveV[i]
    mean1=stepV.mean()
    median=np.median(stepV)
    mean2 =np.mean(stepV)
    std   =np.std(stepV)
    print("mean1  step=",mean1)
    print("mean2  step=",mean2)
    print("median step=",median)
    print("SD     step=",std, "(",std/mean2*100,"% of mean")

def CompV(v1,v2,n):
    dv=np.zeros(n)
    dv=abs(v2-v1)
    mn=np.mean(dv)
    sd=np.std(dv)
    print("Diff Stats=",mn,sd)



bdir="/home/mnb/workspace/ESO/molecfit/work_dir/crires1/telluriccorr_tmp_folder_LJ8w1z"
h0dir="run_1_wdir_lblrtm_range_1_lblrtm_call_1"
f0dir="run_6_wdir_lblrtm_range_1_lblrtm_call_91"
f0dir="run_4_wdir_lblrtm_range_1_lblrtm_call_40"

hdir=os.path.join(bdir,h0dir,"wv_number_0")
fdir=os.path.join(bdir,f0dir,"wv_number_0")


print("ANALYSING DEFAULT")
waveDV,transDV,nptsD= ReadTAPE28(hdir)
print("Finished parsing. Extracted ", nptsD, " pairs of values")


print("ANALYSING ALL")
wdir=os.path.join(hdir,'ALL')
waveAV,transAV,nptsA= ReadTAPE28(wdir)
print("Finished parsing. Extracted ", nptsA, " pairs of values")

print("ANALYSING H2O")
wdir=os.path.join(hdir,'H2O')
waveHV,transHV,nptsH= ReadTAPE28(wdir)
print("Finished parsing. Extracted ", nptsH, " pairs of values")

print("ANALYSING O3")
wdir=os.path.join(hdir,'O3')
waveOV,transOV,nptsO= ReadTAPE28(wdir)
print("Finished parsing. Extracted ", nptsO, " pairs of values")

print("ANALYSING CH4")
wdir=os.path.join(hdir,'CH4')
waveCV,transCV,nptsC= ReadTAPE28(wdir)
print("Finished parsing. Extracted ", nptsC, " pairs of values")

print("Comparing wavenumbers:")
print("H2O with ALL")
CompV(waveHV,waveAV,nptsA)
print("O3 with ALL")
CompV(waveOV,waveAV,nptsA)
print("CH4 with ALL")
CompV(waveCV,waveAV,nptsA)

# Create a product:
transPV=np.zeros(nptsA)
for i in range(nptsA):
    prod=1.0
    prod=prod*transHV[i]
    prod=prod*transOV[i]
    prod=prod*transCV[i]
    transPV[i]=prod

dump4plot("prod.dat",waveAV,transPV,nptsA)
dump4plot("ALL.dat",waveAV,transAV,nptsA)
dump4plot("H2O.dat",waveAV,transAV,nptsA)
#AnalyseData(wave2V,trans2V,npts2)


waveFV,transFV,nptsF= ReadTAPE28(fdir)
dump4plot("F.dat",waveFV,transFV,nptsF)

mol_lst=[]
mol_lst.append(transHV)
mol_lst.append(transOV)
mol_lst.append(transCV)
r_lst=[1.0,1.0,1.0]

tV=GenCombTrans(mol_lst,r_lst,nptsA)
dump4plot("G0.dat",waveHV,tV,nptsH)

lst=CompareTAPE5.CompareTAPE5(fdir,hdir)
r_lst=[]
r_lst.append(lst[0]) # H2O
r_lst.append(lst[2]) # O3
r_lst.append(lst[5]) # CH4
print(r_lst)

tV=GenCombTrans(mol_lst,r_lst,nptsA)
dump4plot("GF.dat",waveHV,tV,nptsH)

exit()


print("ANALYSING CH4")
wdir=os.path.join(hdir,'CH4')
wave2V,trans2V,npts2= ReadTAPE28(wdir)
print("Finished parsing. Extracted ", npts2, " pairs of values")
AnalyseData(wave2V,trans2V,npts2)


exit()

if (npts2==npts):
    CompV(wave2V,waveV,npts)
    CompV(trans2V,transV,npts)
nu1  =3.028e+03
nu2  =3.052e+03
dvset=6.464e-04
m=(nu2-nu1)/dvset
print ("m=",m)
for i in range(10):
    print (i, wave2V[i],waveV[i],wave2V[i]-waveV[i])

#for i in itertools.chain(range(1,5), range(cnt-4,cnt)):
#    print(i, wave_lst[i], trans_lst[i])

