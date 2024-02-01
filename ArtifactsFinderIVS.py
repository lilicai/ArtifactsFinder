 

import os
import sys
import re
from pyfasta import Fasta
import pysam
import argparse
from collections import Counter
from line_profiler import LineProfiler


################################################################################
__author = 'LiaoRui'
__date = '2023/02/14'
__email = 'ruiliao@chosenmedtech.com'


################################################################################
S_LEN = 2
D_LEN = 8
STEM_LEN = 5


################################################################################
def getRegionSequence(ch,start,end,genome):
	### the input start must be 1-base
    seq_dic = Fasta(genome)
    seq = seq_dic[ch][start-1:end]
    return(seq.upper())


################################################################################
def getRegionSequence1(ch,start,end,genome):
	fasta_open = pysam.Fastafile(genome)
	seq = fasta_open.fetch(ch,start-1,end)
	fasta_open.close()
	return(seq.upper())
	
	
################################################################################
def baseTrans(base) :
	BASE_DIC = {'A':'T','T':'A','C':'G','G':'C'}
	return(BASE_DIC[base])
	
	
###############################################################################
def reverSeq(seq) :
	seq = seq.translate(str.maketrans('ACGTacgt','TGCAtgca'))
	seq = seq[::-1]
	return(seq)
	
	
###############################################################################
def transPos(kmerseq_dic,inseq,CH,inseq_start,inseq_end) :
	add_kmerseq_dic = {}
	for k in kmerseq_dic :
		subseq_start,subseq_end,revseq = kmerseq_dic[k]
		map_index = [[r.start(),r.end()] for r in re.finditer(revseq,inseq)]
		re_pos = []
		for r in re.finditer(revseq,inseq) :
			map_index_start,map_index_end = r.start(),r.end()
			map_start = inseq_start + map_index_start
			map_end = inseq_start + map_index_end - 1
			re_pos.append([map_start,map_end])
		key = '_'.join([CH,str(subseq_start),str(subseq_end)])
		add_kmerseq_dic[key] = [CH,subseq_start,subseq_end,re_pos]
	return(add_kmerseq_dic)
	
		
###############################################################################
### step1 : kmer inseq 
def cycleKmerSplit(inseq,CH,inseq_start,inseq_end) :
	SL = len(inseq)
	kmerseq_dic = {}
	for K in range(2, int(SL/2) + 1) :
		kmers = []
		for i in range(0,len(inseq)- K + 1):
			subseq = inseq[i:i + K]
			kmers.append(subseq)
			subseq_start = inseq_start + i
			subseq_end = subseq_start + K - 1
			subseq_len = K
			revseq = reverSeq(subseq)
			if revseq not in inseq : continue
			#if inseq.find(revseq) == -1 : continue
			key = str(i) + '-' + str(K)
			kmerseq_dic[key] = [subseq_start,subseq_end,revseq]
	return(kmerseq_dic)

		
###############################################################################
### step2 : merge overlap : for >=2 map postion
def split2One(kmerseq_dic):
	pos_list = kmerseq_dic.keys()
	onebyone_list = []
	del_list = []
	for k in kmerseq_dic :
		Ch,Start,End,MapPos = kmerseq_dic[k]
		for map_pos in MapPos :
			map_k = Ch + '_' + str(map_pos[0]) + '_' + str(map_pos[1])
			if map_pos[0] == Start and map_pos[1] == End : continue
			if map_pos[0] <= End <= map_pos[1] : continue
			if map_pos[0] <= Start <= map_pos[1] : continue
			onebyone_list.append([Ch,Start,End,map_pos[0],map_pos[1]])

	for K in onebyone_list :
		Ch,Start,End,mStart,mEnd = K
		for k in onebyone_list :
			ch,start,end,mstart,mend = k
			if K != k and Start >= start and End <= end \
			and Start - start == mend - mEnd and End - end == mstart - mStart \
			and K not in del_list :
				del_list.append(K)
	for d in del_list : onebyone_list.remove(d)
	return(onebyone_list)


###############################################################################
### step3 : del short seq, include in long sequence
def delShortSeq(onebyone_list) : 
	del_list = []
	for S in onebyone_list :
		Ch,Start,End,Pos = S
		for s in onebyone_list :
			ch,start,end,pos = s
			if S != s and Start >= start and End <= end \
				and S not in del_list and len(Pos) == 1 :
				del_list.append(S)
	for d in del_list : onebyone_list.remove(d)
	return(onebyone_list)
	
	
###############################################################################
### step4 : merge sub to long 
def checkPoly(k1,k2,inseq,inseq_start):
	ch1,start1,end1,mstart1,mend1 = k1
	ch2,start2,end2,mstart2,mend2 = k2
	seq1 = inseq[start1-inseq_start:end1-inseq_start+1]
	seq2 = inseq[start2-inseq_start:end2-inseq_start+1]
	polyBase = seq1[-1]
	polySeq = ''
	#CTCCC CCCGCAG
	#for i in range(len(seq1),1,-1):
	#	polySeq = polyBase * i
	#	if polySeq in seq1 and seq1.endswith(polySeq) : break
	
	for i in range(len(seq1)-1,1,-1) :
		polySeq = seq1[i:-1]
		if seq1[i] != polyBase: polySeq = seq1[i+1:-1];break
		
	if len(polySeq) >= 2 and seq1.endswith(polySeq):
		start2 += len(polySeq)
		mend2 -= len(polySeq)
	poly = False
	if len(polySeq) == len(seq1) : poly = True
	return(poly,[ch1,start1,end1,mstart1,mend1],[ch2,start2,end2,mstart2,mend2])

	
###############################################################################
def mergeNearby(onebyone_list,inseq,inseq_start) :
	map_list = []
	tmp_list = []
	del_list = []
	gap_dic = {}
	for k1 in onebyone_list :
		ch1,start1,end1,mstart1,mend1 = k1		
		for k2 in onebyone_list :
			ch2,start2,end2,mstart2,mend2 = k2
			
			if k1 == k2 : continue
			if [ch1,start1,end1] == [ch2,mstart2,mend2] : continue
			if [ch1,mstart1,mend1] == [ch1,start1,end1] : continue
			if [ch1,start1,end1] ==  [ch1,mstart1,mend1] : continue
			if [ch1,mstart1,mend1] == [ch2,mstart2,mend2] : continue
			
			### ploy
			poly,k1,k2 = checkPoly(k1,k2,inseq,inseq_start)
			if poly : continue	#all is a base,for example AAAAAA
			ch1,start1,end1,mstart1,mend1 = k1
			ch2,start2,end2,mstart2,mend2 = k2
			if mstart1 <= mend2 <= mend1 and mstart1 <= mstart1 <= mend1 : continue
			if mstart2 < mend1 < mend2 and mstart2 < mstart1 < mend2 : continue
			if start2 < end1 <= end2 : continue
			#['chr1', 22157625, 22157631, 22157367, 22157373] ['chr1', 22157632, 22157630, 22157367, 22157365]
			if mstart1 > mend1 or mstart2 > mend2 : continue
			#['chr1', 66037795, 66037800, 66037986, 66037991] ['chr1', 66037798, 66037799, 66037984, 66037985]
			if start1 < start2 < end1 and start1 < end2 < end1 : continue 
			
			dis1 = start2 - end1
			dis2 = mstart1 - mend2
			if abs(dis1) <= 2 and abs(dis2) <= 2 :
				### del raw pos
				del_list.append(k1)
				del_list.append(k2)
				
				### creat new merge pos
				l_start = start1;l_end = end2
				r_start = mstart2;r_end = mend1
				new_k1 = [ch1,l_start,l_end,r_start,r_end]
				new_k2 = [ch1,r_start,r_end,l_start,l_end]
				if new_k1 in map_list or new_k2 in map_list : continue
				map_list.append(new_k1)
				
				### indel 
				if (abs(dis1) == 1 and abs(dis2) == 2) or (abs(dis1) == 2 and abs(dis2)== 1):
					### find gap pos
					gap_info = {'long':[],'short':[]}
					### gap in right
					if end1 + 1 < start2 :
						lg_start = end1 + 1;lg_end = start2 - 1
						#lg_base = getRegionSequence(ch2,lg_start,lg_end,genome)
						lg_base = inseq[lg_start-inseq_start:lg_end-inseq_start+1]
						gap_info['long'] = [ch2,lg_start,lg_end,lg_base,'-']
						#print(11111,lg_base)
						gap_info['short'] = [ch2,mend2,mend2,'-',baseTrans(lg_base)]
					### gap in left
					elif mend2 + 1 < mstart1 :
						rg_start = mend2 + 1
						rg_end = mstart1 - 1 
						#print(22222222,k1,k2)
						#rg_base = getRegionSequence(ch2,rg_start,rg_end,genome)
						rg_base = inseq[rg_start-inseq_start:rg_end-inseq_start+1]
						gap_info['long'] = [ch2,rg_start,rg_end,rg_base,'-']
						#print(222222222,rg_base,rg_start,inseq_start,rg_end,inseq_start+1)
						gap_info['short'] = [ch2,end1,end1,'-',baseTrans(rg_base)]
					gap_dic['_'.join(list(map(str,new_k1)))] = gap_info
				### snp
				elif abs(dis1) == 2 and abs(dis2) == 2 :
					gap_info = {'left':[],'right':[]}
					#gap_base1 = getRegionSequence(ch1,end1+1,start2-1,genome)
					#gap_base2 = getRegionSequence(ch2,mend2+1,mstart1-1,genome)
					gap_base1 = inseq[end1+1-inseq_start:start2-1-inseq_start+1]
					gap_base2 = inseq[mend2+1-inseq_start:mstart1-1-inseq_start+1]
					gap_info['left'] = [ch1,end1+1,start2-1,gap_base1,reverSeq(gap_base2)]
					gap_info['right'] = [ch2,mend2+1,mstart1-1,gap_base2,reverSeq(gap_base1)]
					gap_dic['_'.join(list(map(str,new_k1)))] = gap_info
				else :
					new_k1 = [ch1,start1,end1,mstart1,mend1]
					new_k2 = [ch1,mstart1,mend1,start1,end1]
					if new_k1 in map_list or new_k2 in map_list : continue
					map_list.append(new_k1)
	return(map_list,gap_dic,del_list)

	
###############################################################################
def mergeLong(onebyone_list,genome,inseq,inseq_start) :
	onebyone_list.sort(key = lambda x : x[1])
	onebyone_list,gap_dic,del_list = mergeNearby(onebyone_list,inseq,inseq_start)
	
	def dealInDel():
		base_info = ''
		min_s_len = 0
		if left_map_len < right_map_len :
			left_info = gap_info['short']
			right_info = gap_info['long']
			#left_base1 = getRegionSequence(ch,start,left_info[1],genome)
			#left_base2 = getRegionSequence(ch,left_info[1]+1,end,genome)
			#right_base = getRegionSequence(ch,mstart,mend,genome)
			left_base1 = inseq[start-inseq_start:left_info[1]-inseq_start+1]
			left_base2 = inseq[left_info[1]+1-inseq_start:end-inseq_start+1]
			right_base = inseq[mstart-inseq_start:mend-inseq_start+1]
			left_base = left_base1 + '-' + left_base2
			min_s_len = min(len(left_base1),len(left_base2))
		else :
			left_info = gap_info['long']
			right_info = gap_info['short']
			#left_base = getRegionSequence(ch,start,end,genome)
			#right_base1 = getRegionSequence(ch,mstart,right_info[1],genome)
			#right_base2 = getRegionSequence(ch,right_info[1]+1,mend,genome)
			left_base = inseq[start-inseq_start:end-inseq_start+1]
			right_base1 = inseq[mstart-inseq_start:right_info[1]-inseq_start+1]
			right_base2 = inseq[right_info[1]+1-inseq_start:mend-inseq_start+1]
			
			right_base = right_base1 + '-' + right_base2
			min_s_len = min(len(right_base1),len(right_base1))
		return(left_base,right_base,min_s_len,left_info,right_info)
	
	def dealSNP():
		base_info = ''
		min_s_len = 0
		left_info = gap_info['left']
		right_info = gap_info['right']
		
		#left_base = getRegionSequence(ch,start,end,genome)
		#left_base1 = getRegionSequence(ch,start,left_info[1]-1,genome)
		#left_base2 = getRegionSequence(ch,left_info[1]+1,end,genome)
		left_base = inseq[start-inseq_start:end-inseq_start+1]
		left_base1 = inseq[start-inseq_start:left_info[1]-1-inseq_start+1]
		left_base2 = inseq[left_info[1]+1-inseq_start:end-inseq_start+1]
		
		#right_base = getRegionSequence(ch,mstart,mend,genome)
		#right_base1 = getRegionSequence(ch,mstart,right_info[1]-1,genome)
		#right_base2 = getRegionSequence(ch,right_info[1]+1,mend,genome)
		right_base = inseq[mstart-inseq_start:mend-inseq_start+1]
		right_base1 = inseq[mstart-inseq_start:right_info[1]-1+inseq_start+1]
		right_base2 = inseq[right_info[1]+1-inseq_start:mend-inseq_start+1]
			
		#base_info = [left_base,right_base]
		min_s_len = min(len(left_base1),len(left_base2),len(right_base1),len(right_base1))
		return(left_base,right_base,min_s_len,left_info,right_info)
		
		
	###########################################################################
	last_seq = {}
	for k in gap_dic :
		#print(k,gap_dic[k])
		gap_info = gap_dic[k]
		ch,start,end,mstart,mend = k.split('_')
		start = int(start);end = int(end)
		mstart = int(mstart);mend = int(mend)
		left_map_len = end - start + 1
		right_map_len = mend - mstart + 1
		map_len = max(left_map_len,right_map_len)
		
		out_string = []
		info1 = '_'.join([ch,str(start),str(end),str(left_map_len)])
		info2 = '_'.join([ch,str(mstart),str(mend),str(right_map_len)])
		out_string = [info1,info2]
		
		### filter stem length
		stem_len = mstart - end if end  <= mstart else start - mend
		if stem_len < STEM_LEN : continue
		
		### filter double length
		if map_len < D_LEN : continue		
		### 		
		if 'short' in gap_info and 'long' in gap_info : 
			left_base,right_base,min_s_len,left_info,right_info = dealInDel()
		elif 'left' in gap_info and 'right' in gap_info :
			left_base,right_base,min_s_len,left_info,right_info = dealSNP()
		else :
			continue
		
		### filter single length
		if min_s_len < S_LEN : continue	
		
		if start > mend : 
			out_string = [str(map_len),info2,info1,'_'.join(list(map(str,right_info))),
					'_'.join(list(map(str,left_info))),right_base,left_base]
		else :
			out_string = [str(map_len),info1,info2,'_'.join(list(map(str,left_info))),
					'_'.join(list(map(str,right_info))),left_base,right_base]
		
		### add at 20230508
		compare = '|'.join(out_string[3:5])
		if compare not in last_seq :
			last_seq[compare] = out_string 
		else : 
			if int(last_seq[compare][0]) < map_len :
				last_seq[compare] = out_string 
		#print(11111,compare)
		#last_seq.append(out_string)
	#last_seq.sort()
	return(last_seq)
		
	
################################################################################		
def basePrect(inseq):
	base_dic = Counter(inseq)
	basebias = False
	base_ratio_dic = {}
	for base in base_dic :
		base_ratio = base_dic[base]/len(inseq)
		base_ratio_dic[base] = base_ratio
		if abs(base_ratio - 0.25) >= 0.1 :
			basebias = True
			#break
	return(basebias,base_ratio_dic)

	
################################################################################		
def choseLong(CH,rawseq_start,rawseq_end,last_seq) :
	#print(1111,last_seq)
	rawseq_info = CH + '_' + str(rawseq_start) + '_' + str(rawseq_end)
	out_string_list = []
	for k in last_seq :
		s = last_seq[k]
		out_string = rawseq_info + '\t' + '\t'.join(s)
		out_string_list.append(out_string)
		#print(out_string)
	return(out_string_list)
	
	
################################################################################		
def run(arg) :
	inbed = arg['bed']
	outfile = os.path.basename(inbed) + '.out'
	genome = arg['genome']
	extend = arg['length']	
	###
	with open(inbed,'r') as f, open(outfile,'w') as fo :
		for line in f :
			if line.startswith('#') : continue
			if line.startswith('\n') : continue
			CH,START,END = line.strip().split('\t')[0:3]
			if '_' in CH : continue 
			
			###
			rawseq_start = int(START) + 1; rawseq_end = int(END) 
			inseq_start = rawseq_start - extend; inseq_end = rawseq_end + extend
			inseq = getRegionSequence(CH,inseq_start,inseq_end ,genome)
			#basebias,base_ratio_dic = basePrect(inseq)
			#print(CH,rawseq_start,rawseq_end)
			#if basebias : continue			
	
			### step1 : kmer seq 
			kmerseq_dic = cycleKmerSplit(inseq,CH,inseq_start,inseq_end)
			kmerseq_dic = transPos(kmerseq_dic,inseq,CH,inseq_start,inseq_end)
			#for k in kmerseq_dic : print(k,kmerseq_dic[k])
			
			### step2 : split >=2 map position to one-by-one format
			onebyone_list = split2One(kmerseq_dic)
			#for k in onebyone_list : print(k)
			
			### step3 : merge sub to long 
			last_seq = mergeLong(onebyone_list,genome,inseq,inseq_start)

			### step4 : print out the best 
			out_string_list = choseLong(CH,rawseq_start,rawseq_end,last_seq)
			if out_string_list : fo.write('\n'.join(out_string_list) + '\n')
	
	
################################################################################
if __name__ == '__main__' :
	parser = argparse.ArgumentParser()
	parser.add_argument('-b', '--bed', help = 'The bed file,0-base',required = True)
	parser.add_argument('-g', '--genome', help = 'The genome file',required = True) 
	parser.add_argument('-l', '--length', help = 'The extension length',
		type = int,default = 50
	) 

	arg = vars(parser.parse_args())
	run(arg)
