## Streamline Process of Ordering, Tagging, and Outputting a Product-Locus Tag-Translation Chart for each Genome
## Author: Alexander Farrell
## July 20, 2016

import os
import re
from optparse import OptionParser
import csv

options = OptionParser(usage='%prog input output ',
                       description="Specify input gbk file and output file")

options.add_option("-i","--infile",dest="inputfile",
                   help="Input file (.gbk)")
options.add_option("-o","--outfile",dest="outputfile",
                   help="output file (.gbk)")
                

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


## Gets a single list of the entire sequence, removing spaces and new line characters
def getseqlist(file):
	with open(file, 'r') as f:
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
	a = open(newgbk,'w')
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



## The main function will create an ordered GBK file with RNA genes placed correctly and return it as Intermediate.gbk
## It then incorporates the locus tags into the Intermediate.gbk and returns the user-inputted GBK file
## The main function then creates a csv spreadsheet of a Product-Locus Tag-Translation
## MUST ALTER LINE 143 FOR SCRIPT
def main():
	opts,args = options.parse_args()
	correctorder,seq = readfile(opts.inputfile)
	writenewgbk(opts.inputfile,correctorder,seq,"Intermediate.gbk")
	seq = getseqlist("Intermediate.gbk")
	Ndict = obtaindict(seq)
	asclist = getincreasinglist(Ndict)
	lookup = makelookupdict(Ndict)
	insertlocustag("Intermediate.gbk",asclist,lookup,opts.outputfile)
	CSVOUT = ""
	CSVOUT = opts.outputfile[:-3]
	CSVOUT += "csv"
	locus = readgbkprod(opts.outputfile)
	trans = readgbkprot(opts.outputfile)
	makecsv(trans,locus,CSVOUT)

if __name__ == '__main__':
	main()
	
                   
                   