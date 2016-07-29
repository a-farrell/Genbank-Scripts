## Post RAST output, this script will Reorder, Tag, and Output either FASTA or CSV for a specific genome
## author: Alexander Farrell
## July 28, 2016

import os
import re
from optparse import OptionParser
import csv
import makeDnaAfirst
import OrderandTagGBK_OutputCSV
import ordergbkgenes

options = OptionParser(usage='%prog -i input -o output -p locus_prefix',
                       description="Specify input gbk file, output file and a prefix for locus tags")

options.add_option("-i","--infile",dest="inputfile",
                   help="Input file (.gbk)")
options.add_option("-o","--outfile",dest="outputfile",
                   help="output file (.gbk)")
options.add_option("-p","--prefix",dest="locus_prefix",help="prefix for all locus tags")
options.add_option("--fasta",dest="write_fasta",
					action="store_true", default=False,
					help="specify if fasta output is desired.")
					
#### ORDER GBK GENES FUNCTIONS

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

#### MAKE DNAA FIRST FUNCTIONS

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
		
### LOCUS TAG AND OUTPUT FUNCTIONS

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
def insertlocustag(oldgbk,ascendinglist,lookupdict,newgbk,locus):
	f=open(oldgbk)
	lines = f.readlines()
	f.close()
	genecount = 5
	
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
	a.write("//")
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

# Print dictionary in two columns of out csv file with protein translations
def makecsv(prot,tags,otherfilename):
	writer = csv.writer(open(otherfilename,'wb'))
	writer.writerow(["Locus Tag","Product","Protein Translation"])
	for key, value in tags.items():
		if key in prot:
			writer.writerow([key,value,prot[key]])
		else:
			writer.writerow([key,value,"X"])

## if the --fasta option is used via command line, the script will output a fasta file with all proteins, not a CSV file with a table
def makefasta(prot,tags,otherfilename):
	writer = open(otherfilename,"wb")

	for key, value in tags.items():
		fastaheader=">"+key+"|"+value+"\n"
		writer.write(fastaheader)
		if key in prot:
			writer.write(prot[key]+"\n")
		else:
			writer.write("X\n")
			
	writer.close()
	
					
#### MAIN FUNCTION					


def main():
	opts,args = options.parse_args()
	orderdict, NTsequence = readfile(opts.inputfile)
	writenewgbk(opts.inputfile,orderdict,NTsequence,"I1.gbk")
	head, seq = getheaderandsequence("I1.gbk")
	totalseqlen = int(findlastNT(head))
	genes, highestindex = creategenedict("I1.gbk")
	DnaAlocation = findDnaA(genes)
	newdict, DnaACDS = reorderdictionary(genes,DnaAlocation,highestindex)
	firstNTofDnaA = int(findNTlocationofDnaA(DnaACDS))
	finalsequence = createfinalsequence(seq,firstNTofDnaA)
	finaldict = changeNTlocations(newdict,highestindex,totalseqlen)
	writegbk("I2.gbk",finaldict,head,finalsequence)
	removeemptylines("I2.gbk","I3.gbk")
	locustags_prefix=str(opts.locus_prefix)
	seqlist = getseqlist("I3.gbk")
	Ndict = obtaindict(seqlist)
	asclist = getincreasinglist(Ndict)
	lookup = makelookupdict(Ndict)
	insertlocustag("I3.gbk",asclist,lookup,opts.outputfile,locustags_prefix)
	CSVOUT = ""
	CSVOUT = opts.outputfile[:-3]
	CSVOUT += "csv"
	locus = readgbkprod(opts.outputfile)
	trans = readgbkprot(opts.outputfile)
	if opts.write_fasta==False:
		makecsv(trans,locus,CSVOUT)
	else:
		FASTA = CSVOUT[:-3] + "fasta"
		makefasta(trans,locus,FASTA)

if __name__ == '__main__':
	main()
	