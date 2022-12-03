filename='TAPE12'

fid=open(filename,'rb')

cnt=0
for record in fid:
    if (cnt==0): print('RECODRD=,',record)
    cnt=cnt+1
print('N lines=',cnt)

#record=fid.readline()
#print('RECODRD=,',record)

fid.close()

n=data.length()
