## Alexander Farrell
## July 12 2016

#Edited 12 July 2016
#Defne Surujon

## This program takes all the previously created .csv files from the program write_csv_withproteins.csv
## and creates one large aggregate gene presence/absence chart 

import os 
from optparse import OptionParser
import csv
import sys

options = OptionParser(usage='%prog input output ', description = "Specify input CSV directory (-i) and output CSV file (-o)")

options.add_option("-i","--infile",dest="inputdir",help="Input csv directory (.csv)")
options.add_option("-o","--outfile",dest="outputfile",help="Output TSV file (.tsv)")
options.add_option("-d","--debug", dest="isDebug", action="store_true", help="Debug mode - verbose")


# Read through CSV file and create a product dictionary
#keys are product names, values are SP numbers
def makedicts(file):
        ProdDict={}
        f = open(file,'rt')
        try:
                reader = csv.reader(f)
                for row in reader:
                        if "hypothetical protein" not in row[1]:
                                ProdDict[row[1]] = row[0]
        finally:
                f.close()	
        return ProdDict	


def main(): 
        opts, args = options.parse_args()
        finaldict = {}
        allProducts=[]
        mystrains=[i[:-4] for i in os.listdir(opts.inputdir)]

        if opts.isDebug == True: print("Starting in Debug Mode\n")
        for f in os.listdir(opts.inputdir):
                f = os.path.join(opts.inputdir,f)
                if f.endswith(".csv"):
                        products = makedicts(f)
                        if opts.isDebug == True: print("#products in "+f+": "+str(len(products.keys())))
                        finaldict[f[:-4]]=products
                        allProducts=allProducts+list(products.keys())
        allProductsUnique=list(set(allProducts))
        for product in allProductsUnique:
        	if product == 'Product':
        		allProductsUnique.remove(product)
        if opts.isDebug == True: print(str(len(allProductsUnique))+" unique gene products found (excluding hypotheticals)\n")
        #print(finaldict["CSV_DEF\Strain1"].values())
        #fout=open(opts.outputfile,"w+")
        with open(opts.outputfile,"w+") as fout:
        	header = ['Product']
        	count = 1 
        	while count <=30:
        		header.append("Strain " + str(count))
        		count+=1
        	writer = csv.DictWriter(fout,fieldnames=header)
        	writer.writeheader()
        	for thisproduct in allProductsUnique:
        		forcsv={}
        		forcsv['Product'] = thisproduct
        		for strain in finaldict:
        			if thisproduct in finaldict[strain].keys():
        				thisSP = finaldict[strain][thisproduct] 
        				strainnum = int(thisSP[3:5])
        				forcsv['Strain ' + str(strainnum)] = thisSP 
        		count = 1
        		while count <= 30:
        			if 'Strain ' + str(count) not in forcsv.keys():
        				forcsv['Strain ' + str(count)] = " "
        			count += 1
        		writer.writerow(forcsv)
        				    
        
        
        
        
if __name__ == '__main__':
        main()
