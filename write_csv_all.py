## Alexander Farrell
## July 11 2016

## This program takes an input GBK file and writes an output CSV containing the locus tag 
## and gene product for referencing

## Defne Surujon - July 21, 2016
## Adding option to output fasta proteome instead of csv

import os
import csv 
from optparse import OptionParser

options = OptionParser(usage='%prog input output --fasta',
                       description="Specify input gbk file and output file, and whether the output will be in fasta format")

options.add_option("-i","--infile",dest="inputfile",
                   help="Input file (.gbk)")
options.add_option("-o","--outfile",dest="outputfile",
                   help="output file (.csv)")
options.add_option("--fasta",dest="write_fasta",
			action="store_true", default=False,
			help="specify if fasta output is desired.")                   

# Read gbk file and return a dictionary of locus tag, product entries                   
def readgbkprod(filename):
	f=open(filename)
	lines=f.readlines()
	f.close()
	tags={}
	i=0
	prod=""
	SP=""
	count=0
	for line in lines:
		if i==1 and line[21]!= '/':
			if "repeat_region" in line or "     gene" in line:
				prodnew = prod.replace("\"","")
				prodfin = prodnew.replace("\n","")
				count+=1
				tags[SP] = prodfin
				i=0
			else:
				prod = prod + line[20:]
		if i==1 and line[21] == '/':
			prodnew = prod.replace("\"","")
			prodfin = prodnew.replace("\n","")
			count+=1
			tags[SP] = prodfin
			i=0
		if "/locus_tag" in line:
			SPna = line[33:] 
			SPn = SPna.replace("\"","")
			SP = SPn.replace("\n","")
		if "/product" in line:
			prod = line[31:]
			i=1
	print "There are " + str(count) + " genes."
	return tags

# Read gbk file and return a dictionary of locus tag, protein translation entries 
def readgbkprot(filename):
	f=open(filename)
	lines=f.readlines()
	f.close()
	prot={}
	i=0
	a=0
	SP=""
	seq=""
	for line in lines[:-1]:
		if "     CDS        " in line:
			i=1
		if line[6:9] == "RNA":
			i=0
		if a==1 and i==1 and line[21] != '/':
			seq += line[21:]
		if a==1 and i==1 and line[21] == '/':
			seqnew = seq.replace("\"","")
			seqrealnew = seqnew.replace("\n","")
			prot[SP] = seqrealnew
			i=0
			a=0
		if "/locus_tag" in line and i==1:
			SP = line[33:44]
		if "/translation" in line and i==1:
			seq = line[35:] 
			a=1		
	return prot

# Print dictionary in two columns of out csv file
def makecsv(prot,tags,otherfilename):
	writer = csv.writer(open(otherfilename,'wb'))
	writer.writerow(["Locus Tag","Product","Protein Translation"])
	for key, value in tags.items():
		if key in prot:
			writer.writerow([key,value,prot[key]])
		else:
			writer.writerow([key,value,"X"])

def makefasta(prot,tags,otherfilename):
	writer = open(otherfilename,"w")

	for key, value in tags.items():
		fastaheader=">"+key+"|"+value+"\n"
		writer.write(fastaheader)
		if key in prot:
			writer.write(prot[key]+"\n")
		else:
			writer.write("X\n")
			
	writer.close()

def main():
	opts, args = options.parse_args()
	locus = readgbkprod(opts.inputfile)
	trans = readgbkprot(opts.inputfile)
	if opts.write_fasta==False:
		makecsv(trans,locus,opts.outputfile)
	else:
		makefasta(trans,locus,opts.outputfile)

if __name__ == '__main__':
    main()
		
