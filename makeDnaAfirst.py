## Script to make DnaA the first gene in a GBK file
## Author: Alexander Farrell
## July 25, 2016

## Edited on July 28th 

import os
import re
from optparse import OptionParser
import sys
import tempfile


options = OptionParser(usage='%prog input output ',
                       description="Specify input gbk file and output file")

options.add_option("-i","--infile",dest="inputfile",
                   help="Input file (.gbk)")
options.add_option("-o","--outfile",dest="outputfile",
                   help="Output file (.gbk)")



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
	headers = header.replace("linear","circular")
	return headers, result

## Return the full length of the NT sequence
def findlastNT(header):
	headerarray = header.split("\n")
	for line in headerarray:
		if "     source          " in line:
			locarray = re.findall(r'\d+', line)
	return locarray[1]

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
		elif i==1 and ("     CDS             " in line or line[6:9] == 'RNA'):
			genes[key] = CDS
			CDS = ""
			CDS += line
			key += 1
		elif "ORIGIN" in line or "BASE COUNT" in line:
			genes[key] = CDS
			i = 0
		elif i==1 and "     CDS           " not in line:
			CDS += line
		elif i==1 and line[6:9] != 'RNA':
			CDS += line
	return genes, key	

## Returns the gene number that contains DnaA	
def findDnaA(dictionary):
	for key, value in dictionary.items():
		if "Chromosomal replication initiator protein DnaA" in value:
			return key

## Create new dictionary where DnaA is the first gene
def reorderdictionary(dictionary,num,highest):
	newdict = {}
	CDS = " "
	DnaA = ""
	key = 1
	DnaA = dictionary[num]
	for i in range(1,highest+1):
		CDS = dictionary[num]
		newdict[key] = CDS
		key +=1
		num += 1
		if num > highest:
			num = 1
	return newdict, DnaA

## Determine what is the NT location of DnaA to reorder the sequence
def findNTlocationofDnaA(CDS):
	StartNT = 0
	EndNT = 0
	array = []
	array = CDS.split("\n")
	NTarray = re.findall(r'\d+', array[0])
	StartNT = NTarray[0]
	EndNT = NTarray[1]
	return StartNT


## Reorder NT sequence
def createfinalsequence(prevseq,DnaAloc):
	first = prevseq[DnaAloc-1:]
	end = prevseq[:DnaAloc]
	finalseq = first + end
	return finalseq
			
## At this point, we have: a) put RNA genes in correct position, b) reorganize the sequence and c) reorganize CDS portion so DnaA is the first gene
## All that is left to do is change the locations of the NTs in the CDS and print the sequence

## Look through previous dictionary and create new one with correct NT locations 
def changeNTlocations(dictionary,highest,sequencelength):
	finaldict = {}
	index = 1
	pick1 = "     CDS             "
	pick2 = "     tRNA            "
	pick3 = "     rRNA            "
	for i in range(highest):
		if index == 1:
			CDS = dictionary[index]
			CDSarray = CDS.split("\n")
			NTarray = re.findall(r'\d+', CDSarray[0])
			StartNT = NTarray[0]
			EndNT = NTarray[1]
			distance = int(EndNT) - int(StartNT)
			distance +=1
			if "CDS" in CDSarray[0]:
				firstline = pick1 + str('1') + ".." + str(distance) + "\n" 
			elif "tRNA" in CDSarray[0]:
				firstline = pick2 + str('1') + ".." + str(distance) + "\n" 
			elif "rRNA" in CDSarray[0]:
				firstline = pick3 + str('1') + ".." + str(distance) + "\n" 
			finalCDS = ""
			finalCDS += firstline
			for line in CDSarray[1:]:
				finalCDS += line + "\n"
			finaldict[index] = finalCDS	
		else:	
			previousCDS = dictionary[index-1]
			prevCDSarray = previousCDS.split("\n")
			prevNTarray = re.findall(r'\d+', prevCDSarray[0])
			prevSNT = int(prevNTarray[0])
			prevENT = int(prevNTarray[1])
			currentCDS = dictionary[index]
			currentCDSarray = currentCDS.split("\n")
			currentNTarray = re.findall(r'\d+', currentCDSarray[0])
			currentSNT = int(currentNTarray[0])
			currentENT = int(currentNTarray[1])
			## Handle the issue of overlapping from end of genome to start of genome 
			currentgenelength = currentENT - currentSNT
			CDSfromfinaldict = finaldict[index-1]
			CDSarrayfromfinaldict = CDSfromfinaldict.split("\n")
			NTarrayfromfinaldict = re.findall(r'\d+', CDSarrayfromfinaldict[0])
			endingNTfromfinaldict = int(NTarrayfromfinaldict[1])
			if currentSNT < 1000000 and prevENT > 1000000:   # Means that you have overlapped
				newSNT = (sequencelength - prevENT) + currentSNT + endingNTfromfinaldict
				newENT = newSNT + currentgenelength
			else:
				newSNT = endingNTfromfinaldict + (currentSNT - prevENT)
				newENT = newSNT + currentgenelength
			if "CDS" in currentCDSarray[0]:
				if "(" in currentCDSarray[0]:
					firstline = pick1 + "complement(" + str(newSNT) + ".." + str(newENT) + ")" + " \n" 
				else:
					firstline = pick1 + str(newSNT) + ".." + str(newENT) + " \n" 
			elif "tRNA" in currentCDSarray[0]:
				if "(" in currentCDSarray[0]:
					firstline = pick2 + "complement(" + str(newSNT) + ".." + str(newENT) + ")" + " \n"
				else:
					firstline = pick2 + str(newSNT) + ".." + str(newENT) + " \n" 
			elif "rRNA" in currentCDSarray[0]:
				if "(" in currentCDSarray[0]:
					firstline = pick3 + "complement(" + str(newSNT) + ".." + str(newENT) + ")" + " \n"
				else:
					firstline = pick3 + str(newSNT) + ".." + str(newENT) + " \n" 
			finalCDS = ""
			finalCDS += firstline
			for line in currentCDSarray[1:]:
				finalCDS += line + "\n"
			finaldict[index] = finalCDS
		index += 1
	return finaldict

# Write header, CDS information, then sequence into output gbk file
def writegbk(file,finaldict,header,sequence):
	a = open(file,'w+')
	a.write(header)
	for key, value in finaldict.items():
		CDS = value
		PartCDS = CDS.split('\n')
		for line in PartCDS:
			if "/transl_table" in line:
				liner = line.replace("\n","")
				a.write(liner)
			else:
				a.write(line)
				a.write("\n")
	a.write("ORIGIN")
	a.write("\n")
	linetowrite = ""
	numbertowrite = 1
	linecount = 0 
	for char in sequence:
		if numbertowrite % 60 == 1:
			numstring = str(numbertowrite)
			lenstring = len(numstring)
			spaces = 9-lenstring
			linetowrite += spaces*" " + numstring + " " + char
		elif numbertowrite % 60 == 0:
			linetowrite += char
			a.write(linetowrite)
			a.write("\n")
			linetowrite = ""
			linecount = -1
		elif linecount % 10 == 0:
			linetowrite += char + " "
		else:
			linetowrite += char
		linecount += 1
		numbertowrite += 1
	a.write(linetowrite)
	a.write("\n")
	a.write("//")
	a.close()

## Remove any blank lines that appear in the GBK file
def removeemptylines(file,otherfile):
	a = open(file)
	lines = a.readlines()
	a.close()
	f = open(otherfile,'w+')
	for line in lines:
		if line.strip():
			f.write(line)
	f.close()
		
## Main Function

def main():
	opts, args = options.parse_args()
	head, seq = getheaderandsequence(opts.inputfile)
	totalseqlen = int(findlastNT(head))
	genes, highestindex = creategenedict(opts.inputfile)
	number = findDnaA(genes)
	newdict, DnaACDS = reorderdictionary(genes,number,highestindex)
	firstNTofDnaA = int(findNTlocationofDnaA(DnaACDS))
	finalsequence = createfinalsequence(seq,firstNTofDnaA)
	finaldict = changeNTlocations(newdict,highestindex,totalseqlen)
	writegbk("Intermediate.gbk",finaldict,head,finalsequence)
	removeemptylines("Intermediate.gbk",opts.outputfile)
	
	
if __name__ == '__main__':
    main()
	
	
		
			
		  
		
		
	
                