filename='TAPE5'

fid=open(filename,'r')
data=fid.readlines()
cnt=0
for line in data:
    cnt=cnt+1
    if (cnt<6) :
        print(cnt, ") ", line)
print('N lines=',cnt)

#record=fid.readline()
#print('RECODRD=,',record)

fid.close()

