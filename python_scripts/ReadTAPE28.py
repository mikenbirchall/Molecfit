import itertools
# Read all lines in asci file TAPE28
filename='TAPE28'
fid=open(filename,'r')
lines=fid.readlines()
fid.close()

# Iterate through all lines and extract Wavenumbers and Transmission values
cnt=0;         # No of lines with values extracted
extracting=0;  # Boolean if in extracting mode
wave_lst=[];   # The wavenumber extracted values list
trans_lst=[]   # The transmission extracted values list
for line in lines:

   # Skip empty lines
    if (len(line.strip())==0):
        continue

    # If this line contains "WAVENUMBER" and "TRANSMISSION" then start extracting
    # on the next iteration
    if (line.find("WAVENUMBER")!=-1 and line.find("TRANSMISSION")!=-1):
            extracting=1
            continue

    # Skip if we are not extracting
    if (not extracting):
        continue

    # Real Values are to be extracted from this line
    cnt+=1
    lst=line.split()
    wave_lst.append (lst[0])
    trans_lst.append(lst[1])

print("Finished parsing. Extracted ",cnt, " pairs of values")

for i in itertools.chain(range(1,5), range(cnt-4,cnt)):
    print(i, wave_lst[i], trans_lst[i])
