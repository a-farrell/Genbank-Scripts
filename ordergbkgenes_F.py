## This program takes an unordered GBK file and orders the genes, both CDS and RNA, in order
## of their chromosomal locations starting from the ORI

import os
import re
from optparse import OptionParser

options = OptionParser(usage='%prog input output ',
                       description="Specify input gbk file and output file")

options.add_option("-i","--infile",dest="inputfile",
                   help="Input file (.gbk)")
options.add_option("-o","--outfile",dest="outputfile",
                   help="output file")

# Read file and create a dictionary with the starting BP for each gene as the key and the 
# corresponding gene info as its value 
# Returns a dictionary and the sequence as a string 
def readfile(oldgbk):
	f = open(oldgbk)
	lines = f.readlines()
	f.close()
	order = {}
	i=0
	geneinfo=""
	sequence=""
	for line in lines:
		if "ORIGIN" in line:
			order[startloc] = geneinfo 
			i=2
			sequence +=line
		elif i==2 and "ORIGIN" not in line:
			sequence +=line	
		if "     CDS             " in line and i==0:
			startloc = int(re.search('\d+',line).group())
			geneinfo += line
			i=1
		elif "     CDS             " in line or "     tRNA            " in line or "     rRNA            " in line and i==1:
			order[startloc] = geneinfo
			startloc = int(re.search('\d+',line).group())
			geneinfo = ""
			geneinfo += line
		elif "     CDS             " not in line or "     tRNA            " not in line or "     rRNA            " not in line and i==1:
			geneinfo += line
			 
		
	return order, sequence

# Creates ordered gbk file by using a counter to take genes out in ascending order and writes them
def writenewgbk(oldgbk,order,sequence,newgbk):
	f = open(oldgbk)
	lines = f.readlines()
	f.close()
	a = open(newgbk,'w+')
	counter = 0
	for line in lines:
		##if "     CDS             " not in line:
		##	a.write(line)
		if "     source" in line:
			startend = re.findall('\d+',line)
			endbp = int(startend[1])
		if "     CDS             " in line:
			break 
	for i in range(endbp):
		if i in order:
			a.write(order[i])
	a.write(sequence)
	a.close()

def main():
	opts,args = options.parse_args()
	orderdict, NTsequence = readfile(opts.inputfile)
	writenewgbk(opts.inputfile,orderdict,NTsequence,opts.outputfile)

if __name__ == '__main__':
	main()
	
				
					