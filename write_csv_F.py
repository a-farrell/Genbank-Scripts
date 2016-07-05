## This program takes an input GBK file and writes an output CSV containing the locus tag 
## and gene product for referencing

import os
import csv 
from optparse import OptionParser

options = OptionParser(usage='%prog input output ',
                       description="Specify input gbk file and output csv file")

options.add_option("-i","--infile",dest="inputfile",
                   help="Input file (.gbk)")
options.add_option("-o","--outfile",dest="outputfile",
                   help="output file (.csv)")
                   
                   
# Read gbk file and return a dictionary of locus tag, product entries                   
def readgbk(filename):
	f=open(filename)
	lines=f.readlines()
	f.close()
	tags={}
	i=0
	prod=""
	SP=""
	for line in lines:
		if i==1 and line[21]!= '/':
			prod = prod + line[20:]
		if i==1 and line[21]== '/':
			prodnew = prod.replace("\"","")
			tags[SP] = prodnew
			i=0
		if "/locus_tag" in line:
			SP = line[33:44]
		if "/product" in line:
			prod = line[31:]
			i=1
	return tags

# Print dictionary in two columns of out csv file
def makecsv(tags,otherfilename):
	writer = csv.writer(open(otherfilename,'wb'))
	writer.writerow(["Locus Tag","Product"])
	for key, value in tags.items():
		writer.writerow([key,value])



def main():
	opts, args = options.parse_args()
	locus = readgbk(opts.inputfile)
	makecsv(locus,opts.outputfile)

if __name__ == '__main__':
    main()
		