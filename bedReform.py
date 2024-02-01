
import os
import sys

infile = sys.argv[1]
prefix = infile.split('.')[0]
outfile = prefix + '.reform.bed'

###############################################################################
with open(infile,'r') as fi, open(outfile,'w') as fo :
	for line in fi :
		data = line.strip().split('\t')
		ch,start,end = data[0:3]
		start = int(start)
		end = int(end)
		#print(111,ch,start,end,end-start+1)	
		for i in range(0,end-start,200):
			#print(i)
			sub_start = start + i 
			sub_end = sub_start + 200 -1 if sub_start + 200 -1 <= end else end
			#print(222,ch,sub_start,sub_end,sub_end-sub_start+1) 
			fo.write('\t'.join([ch,str(sub_start),str(sub_end)]) + '\n')
