## This program inserts locus tags in a previously-ordered GBK file, according 
## to the spacing of N's between the initial contigs

import os
import re
from optparse import OptionParser

options = OptionParser(usage='%prog input output ',
                       description="Specify an input GBK file (.gbk)(-i) and a new output GBK file (.gbk)(-o)")

options.add_option("-i","--infile",dest="inputfile",
                   help="Input file (.gbk)")
options.add_option("-o","--outfile",dest="outputfile",
                   help="output file (.gbk)")
                   

## Gets a single list of the entire sequence, removing spaces and new line characters
def getseqlist(file):
	f=open(file)
	lines = f.readlines()
	f.close()
	i = 0
	sequence = ""
	for line in lines:
		if i == 1:
			sequence += line[10:]
		if "ORIGIN" in line:
			i = 1
		if "//" in line:
			i = 0
	a = sequence.replace(" ","")
	b = a.replace("\n","")
	return b
	

## Creates a dictionary with the start and end locations of groups on N's in the sequence
def obtaindict(list):
	N = {}
	i=0
	charcount = 1
	start = 0
	end = 0
	for char in list:
		if char == 'n' and i == 0:
			start = charcount
			i=1
		if char != 'n' and i == 1:
			i=0
			end = charcount
			N[start] = end
		charcount += 1
	return N


## Takes an input dictionary and creates a list of the start locations of groups of N's in ascending order
def getincreasinglist(dict):
	listofkeys=[]
	list1 = list(dict.keys())
	for x in list1:
		num = int(x)
		listofkeys.append(num)
	sortedlist = sorted(listofkeys,key=int)
	return sortedlist


## Creates look-up dictionary with the number of genes that could possibly be in the sequence of N's
def makelookupdict(dict):
	newdict={}
	for x in dict:
		startnum = int(x)
		endnum = int(dict[x])
		distance = endnum - startnum
		## Ensure to account for the fact that some group of N's may be less than 100 bp 
		if distance < 100:
			distancedividedby100 = 1
		else:
			distancedividedby100 = distance/100
		newdict[startnum] = distancedividedby100
	return newdict

	
## Takes previously made ascending list and look-up dictionary to add locus tags into a user-specified 
## new GBK file according to the groups of N's in the sequence
def insertlocustag(oldgbk,ascendinglist,lookupdict,newgbk):
	f=open(oldgbk)
	lines = f.readlines()
	f.close()
	genecount = 5
	locus = "SPS01"     #  Must be altered with each new genome	
	a = open(newgbk,'w+')
	for line in lines[:-1]:
		if "     CDS             " not in line and line[6:9] != "RNA":
			a.write(line)
		else:
			if len(ascendinglist) != 0:
				loclist = re.findall('\d+',line)
				startloc = int(loclist[0])
				if startloc > ascendinglist[0]:
					startN = ascendinglist.pop(0)
					numberofgenes = lookupdict[startN]
					increment = numberofgenes*5
					genecount += increment
			a.write(line)
			if genecount < 10:
				a.write("                     /locus_tag=\"" + locus + "_0000" + str(genecount) + "\"")
			elif genecount < 100:
				a.write("                     /locus_tag=\"" + locus + "_000" + str(genecount) + "\"")
			elif genecount < 1000:
				a.write("                     /locus_tag=\"" + locus + "_00" + str(genecount) + "\"")
			elif genecount < 10000:
				a.write("                     /locus_tag=\"" + locus + "_0" + str(genecount) + "\"")
			else:
				a.write("                     /locus_tag=\"" + locus + "_" + str(genecount) + "\"")
			a.write("\n")
			genecount += 5
	a.close()
	      

## Main function that calls all previously written functions			
def main():
	opts, args = options.parse_args()
	seq = getseqlist(opts.inputfile)
	Ndict = obtaindict(seq)
	##print(Ndict)
	asclist = getincreasinglist(Ndict)
	##print(asclist)
	##print(len(asclist))
	lookup = makelookupdict(Ndict)
	##print(lookup)
	insertlocustag(opts.inputfile,asclist,lookup,opts.outputfile)

if __name__ == '__main__':
	main()
	
