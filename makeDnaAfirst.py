## Script to make DnaA the first gene in a GBK file
## Author: Alexander Farrell
## July 25, 2016

import os
import re
from optparse import OptionParser

options = OptionParser(usage='%prog input output ',
                       description="Specify input gbk file and output file")

options.add_option("-i","--infile",dest="inputfile",
                   help="Input file (.gbk)")
#options.add_option("-o","--outfile",dest="outputfile",
#                   help="Output file (.gbk)")



## First order the GBK genes



## Removes new line characters, spaces, numbers, and ORIGIN from sequence, so you can index freely
def removeextra(string):
	strin = string.replace('ORIGIN',"")
	stri = strin.replace("\n","")
	str = stri.replace(" ","")
	result = ""
	for ch in str:
		if ch.isdigit():
			pass
		else:
			result += ch
	return result	

## Returns strings that contain the header for the final GBK file as well as a full-sequence 		
def getheaderandsequence(gbk):
	f = open(gbk)
	lines = f.readlines()
	f.close()
	header = ""
	sequence = ""
	a=0 # indicates you are in header
	b=0 # indicates you are in sequence
	for line in lines:
		if "     CDS             " not in line and a==0:
			header += line
		else:
			a=1
		if b==1 and "//" not in line:
			sequence += line
		if "ORIGIN" in line:
			b=1
			sequence += line
	result = removeextra(sequence)
	return header, result

## Create a dictionary with keys = gene # in ordered GBK and values are CDS portions of genes
def creategenedict(file):
	f = open(file)
	lines = f.readlines()
	f.close()
	genes = {}
	key = 1
	i=0
	CDS = ""
	for line in lines:
		if i==0 and ( "     CDS             " in line or line[6:9] == "RNA"):
			i=1
			CDS += line
		elif			
	
	


def main():
	opts, args = options.parse_args()
	head, seq = getheaderandsequence(opts.inputfile)
	genes = creategenedict(opts.inputfile)

main()
	
	
		
			
		  
		
		
	
                