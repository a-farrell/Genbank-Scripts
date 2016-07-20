## Alexander Farrell
## July 12 2016

## In some of our GBK files, there are two SPS numbers with the same products
## and therefore, the products need to be altered


import os 
from optparse import OptionParser
import csv
import sys

options = OptionParser(usage='%prog input output ', description = "Specify input CSV directory (-i) and output CSV directory (-o)")

options.add_option("-i","--infile",dest="inputfile",help="Input csv directory (.csv)")
options.add_option("-o","--outfile",dest="outputfile",help="Output directory CSV file (.csv)")

def changecsv(file):
	f = open(file,'rt')
	d = {}
	count = 2
	try:
		reader = csv.reader(f)
		for row in reader:
			SP = row[0]
			product = row[1]
			if product in d.values():
				product = product + ".1"
				while product in d.values():
					if count <= 10:
						product = product[:-2] + "." + str(count)
					if count > 10:
						product = product[:-3] + "." + str(count)
					count +=1
					
			d[SP] = product 
			count = 2
	finally:
		f.close()
	return d

def main():
	opts, args = options.parse_args()
	d = changecsv(opts.inputfile)
	with open(opts.outputfile,'wb') as f:
		w = csv.writer(f)
		w.writerow(["Locus Tag","Product"])
		for key, value in d.items():
			w.writerow([key,value])
		
	
	
if __name__ == '__main__':
    main()


			
			