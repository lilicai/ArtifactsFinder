
import os
import sys
import re
import collections 
import argparse
from pyfasta import Fasta


__author = 'Liaorui'
__date = '2023/02/14'
__email = 'ruiliao@chosenmedtech.com'


###############################################################################
GAP = 1
DIS = 3 	#two gap length
#LEN = 17 	#subseq length
BASE_DIC = {'#':'#','A':'T','T':'A','C':'G','G':'C'}


###############################################################################
def getSeq(string,i) :
	t = 1
	gap = 0
	length = 0
	gap_list = []
	while 1 :
		if i-t < 0 or i+t >= len(string) : return (gap_list,string[i-t+1:i+t])
		# gap = 0
		if string[i-t] != BASE_DIC[string[i+t]] and gap < GAP : 
			gap_list.append([i-t,i+t])
			gap+=1
			t+=1
			continue
		#gap = 1
		if string[i-t] != BASE_DIC[string[i+t]] and gap == GAP : 
			gap_list.append([i-t,i+t])
			if string[i]=='#' :
				return (gap_list,string[i-t:i+t+1])
			else :
				return (gap_list,string[i-t+1:i+t])
		t+=1


###############################################################################
def longestPalinString(string):
	seq_dic = {}
	j = 0 
	for i in range(1,len(string)-1) :
		gap_list,sequence = getSeq(string,i)
		length = len(sequence.replace('#',''))
		if length in seq_dic :
			seq_dic[length].update({i:[sequence,gap_list]})
		else :
			seq_dic.update({length:{i:[sequence,gap_list]}})
	seq_dic = sorted(seq_dic.items(),key = lambda x : x, reverse = True)
	return(seq_dic)


###############################################################################
def list2Str(ls,j = '\t'):
	return(j.join(list(map(str,ls))))


###############################################################################
def mergeOut(seq_dic,ch,start,end,seq) :
	res_dic = {}
	for i in seq_dic :
		subSeq = seq_dic[i][0]
		gap_list = seq_dic[i][1]
		cleanSeq = subSeq.replace('#','')
		cleanSeqS = start + int(( i -  len(cleanSeq) ) / 2)
		cleanSeqE = cleanSeqS + len(cleanSeq) - 1 
		cleanSeqInfo = list2Str([ch,cleanSeqS,cleanSeqE,cleanSeq],j='-')
		snp_list = []

		#######################################################################
		def sNP_odd(gap):
			gap_pos = start + int(gap / 2 )
			gap_ref = cleanSeq[gap_pos - cleanSeqS]
			gap_alt = BASE_DIC[gap_ref]
			snp = list2Str([ch,gap_pos,gap_pos,gap_ref,gap_alt],j='|')
			snp_list.append([snp])
			

		#######################################################################
		def sNP_even(gap_list):
			gap_1_pos = start + int((gap_list[0]) / 2)	
			gap_2_pos = start + int((gap_list[1]) / 2)	
			gap_1_ref = cleanSeq[gap_1_pos - cleanSeqS]
			gap_2_ref = cleanSeq[gap_2_pos - cleanSeqS]
			gap_1_alt = BASE_DIC[gap_2_ref]
			gap_2_alt = BASE_DIC[gap_1_ref]

			# two snp nearby merge to a mnp
			if abs(gap_1_pos - gap_2_pos) == 1 :
				gap_ref = gap_1_ref + gap_2_ref
				gap_alt = gap_1_alt + gap_2_alt
				snp = list2Str([ch,gap_1_pos,gap_2_pos,gap_ref,gap_alt],j='|')
				snp_list.append([snp])
			else :
				snp1 = list2Str([ch,gap_1_pos,gap_1_pos,gap_1_ref,gap_1_alt],j='|')
				snp2 = list2Str([ch,gap_2_pos,gap_2_pos,gap_2_ref,gap_2_alt],j='|')	
				snp_list.append([snp1,snp2])

	
		#######################################################################
		def sNP_merge(gap_list):
			gap1,gap2 = gap_list
			gap_1_1_pos = start + int((gap_list[0][0]) / 2)	
			gap_1_2_pos = start + int((gap_list[0][1]) / 2)	
			gap_2_1_pos = start + int((gap_list[1][0]) / 2)	
			gap_2_2_pos = start + int((gap_list[1][1]) / 2)	
			
			gap_1_1_ref = cleanSeq[gap_1_1_pos - cleanSeqS]
			gap_1_2_ref = cleanSeq[gap_1_2_pos - cleanSeqS]
			gap_2_1_ref = cleanSeq[gap_2_1_pos - cleanSeqS]
			gap_2_2_ref = cleanSeq[gap_2_2_pos - cleanSeqS]
			
			gap_1_1_alt = BASE_DIC[gap_1_2_ref]
			gap_1_2_alt = BASE_DIC[gap_1_1_ref]
			gap_2_1_alt = BASE_DIC[gap_2_2_ref]
			gap_2_2_alt = BASE_DIC[gap_2_1_ref]

			# two snp nearby merge to a mnp
			gap_1_ref = gap_2_1_ref + gap_1_1_ref
			gap_1_alt = gap_2_1_alt + gap_1_1_alt
			snp1 = list2Str([ch,gap_2_1_pos,gap_1_1_pos,gap_1_ref,gap_1_alt],j='|')

			gap_2_ref = gap_1_2_ref + gap_2_2_ref
			gap_2_alt = gap_1_2_alt + gap_2_2_alt 
			snp2 = list2Str([ch,gap_1_2_pos,gap_2_2_pos,gap_2_ref,gap_2_alt],j='|')
			snp_list.append([snp1,snp2])
	

 		#######################################################################
		### search snps
		### add filter
		### snp-----snp-----snp
		if len(cleanSeq) %2 == 1 :
			sNP_odd(i)
			if len(gap_list) >= 1 and BASE_DIC[cleanSeq[0]] == cleanSeq[-1] :
				sNP_even(gap_list[0])
		else :
			### snpsnp------snpsnp
			if len(gap_list) > 1 and abs(gap_list[0][0] - gap_list[1][0]) == 2 :
				sNP_merge(gap_list)
			else :
				if len(gap_list) >= 1 :
					sNP_even(gap_list[0])
				if len(gap_list) == 2 and BASE_DIC[cleanSeq[0]] == cleanSeq[-1] :
					sNP_even(gap_list[1])
		
		
		#######################################################################
		res_dic[i] = [str(len(cleanSeq)),cleanSeqInfo,snp_list]
	return(res_dic)	

	
################################################################################
def getRegionSequence(ch,start,end,genome):
	### the input start must be 1-base
    seq_dic = Fasta(genome)
    seq = seq_dic[ch][start-1:end]
    return(seq.upper())
	
	
###############################################################################
def reverSeq(seq) :
	seq = seq.translate(str.maketrans('ACGTacgt', 'TGCAtgca'))
	seq = seq[::-1]
	return(seq)
				



###############################################################################
def printOut(data,res_dic,out_info,out_black,black_list) :
	ch,start,end = data
	
	for i in res_dic :
		cleanSeqL,cleanSeqInfo,snp_list = res_dic[i]
		cleanSeq = cleanSeqInfo.split('-')[3]
		resOut = [cleanSeqL,cleanSeqInfo]
		#if int(cleanSeqL) < LEN : continue
		ps_c,ps_s,ps_e,ps_seq = cleanSeqInfo.split('-')
		if int(ps_s) > int(end) or int(ps_e) < int(start) : continue
		#print(res_dic[i])
		for snp in snp_list :
			for ssnp in snp :
				c,s,e,r,a = ssnp.split('|')
				if c == ch and int(s) >= int(start) and int(e) <= int(end) \
					and ssnp not in black_list:
					out_black.write('\t'.join([c,s,e,r,a]) + '\n')
				black_list[ssnp] = 0	
			if len(snp) == 1 :
				resOut+=[snp[0],'-']
			else :
				resOut+=snp
		out_info.write('\t'.join(data + resOut) + '\n')
	return(black_list)

		
################################################################################		
def run(arg) :
	inbed = arg['bed']
	genome = arg['genome']
	extend = arg['length']
	prefix = arg['prefix']
	
	out_info = open(prefix + '.info.txt','w')
	out_black = open(prefix + '.backlist.txt','w')
	black_dic = {}
	
	with open(inbed,'r') as f :
		for line in f :
			if line.startswith('#') : continue
			if line.startswith('\n') : continue			
			ch,start,end = line.strip().split('\t')[0:3]
			#if '_' in ch : continue #chr6_cox_hap2

			#start = int(start) + 1 			# bed is 0-base
			ex_start = int(start) + 1 - extend	# extend 50
			ex_end = int(end) + extend			# extend 
			seq = getRegionSequence(ch,ex_start,ex_end,genome)
			
			seq = '#' + '#'.join(list(seq)) + '#'	#2n+1
			seq_dic = longestPalinString(seq)
			
			# return the longest one seq
			#for s in seq_dic : print(s)
			res_dic = mergeOut(seq_dic[0][1],ch,ex_start,ex_end,seq)		
			black_dic = printOut([ch,start,end],res_dic,out_info,out_black,black_dic)
	out_info.close()
	out_black.close()		
			
			
################################################################################
if __name__ == '__main__' :
	parser = argparse.ArgumentParser()
	parser.add_argument('-b', '--bed', help = 'The bed file,0-base',required = True) 
	parser.add_argument('-g', '--genome', help = 'The genome file',required = True) 
	parser.add_argument('-l', '--length', help = 'The extend length',type = int,default = 50) 
	parser.add_argument('-p', '--prefix', help = 'The outfile prefix',default='out')

	arg = vars(parser.parse_args())
	run(arg)	
