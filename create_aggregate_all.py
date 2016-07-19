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
options.add_option("-o","--outfile",dest="outputfile",help="Output CSV file (.csv)")
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
                                ProdDict[row[1].upper()] = row[0]
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
        	if product == 'PRODUCT':
        		allProductsUnique.remove(product)
        if opts.isDebug == True: print(str(len(allProductsUnique))+" unique gene products found (excluding hypotheticals)\n")
        #print(finaldict["CSV_DEF\Strain1"].values())
        #fout=open(opts.outputfile,"w+")
        refdict = {'SPS01_':'Strain 1','SPS02_':'Strain 2','SPS03_':'Strain 3','SPS04_':'Strain 4','SPS05_':'Strain 5','SPS06_':'Strain 6','SPS07_':'Strain 7','SPS08_':'Strain 8',
        'SPS09_':'Strain 9','SPS10_':'Strain 10','SPS11_':'Strain 11','SPS12_':'Strain 12','SPS13_':'Strain 13','SPS14_':'Strain 14','SPS15_':'Strain 15','SPS16_':'Strain 16','SPS17_':'Strain 17',
        'SPS18_':'Strain 18','SPS19_':'Strain 19','SPS20_':'Strain 20','SPS21_':'Strain 21','SPS22_':'Strain 22','SPS23_':'Strain 23','SPS24_':'Strain 24','SPS25_':'Strain 25','SPS26_':'Strain 26',
        'SPS27_':'Strain 27','SPS28_':'Strain 28','SPS29_':'Strain 29','SPS30_':'Strain 30','SP_RS':'TIGR4','SPH_RS':'19A','SPT_RS':'19F','HMPREF0':'TCH84331/19A','MYY_RS':'ST556','SPNA45':'SPNA45',
        'SPP_RS':'P1031','SPNOXC':'OXC141','SPJ_RS':'JJA','SPNINV':'INV200','INV104':'INV104','HMPREF1':'gamPNI0373','SPG_RS':'G54','SPD_RS':'D39','SPCG_R':'CGSP14','SPN23F':'ATCC','SPAP_R':'AP200',
        'SP7058':'70585','SP670_':'670-6B'}
        
        ## Special case for HMPREF between gamPNI0373 and TCH84331/19A       
        ## R6 is an exception because its locus tags are still old spr       
        ## refdict looks for first six characters of the locus tag to find the appropriate column to write in
                 
        with open(opts.outputfile,"w+") as fout:
        	header = ['Product']
        	for key,value in refdict.items():
        		header.append(value)
        	header.append('R6') 
        	writer = csv.DictWriter(fout,fieldnames=header)
        	writer.writeheader()
        	for thisproduct in allProductsUnique: 
        		forcsv={}
        		forcsv['Product'] = thisproduct
        		for strain in finaldict:
        			if thisproduct in finaldict[strain].keys():
        				thisSP = finaldict[strain][thisproduct] 
        				if 'spr' in thisSP:
        					forcsv['R6'] = thisSP
        				elif 'SP_RS' in thisSP:
        					forcsv['TIGR4'] = thisSP
        				else:
        					smallSP = thisSP[:6]
        					if smallSP == 'HMPREF':
        						smallSP = thisSP[:7]
        						Strain = refdict[smallSP]
        						forcsv[Strain] = thisSP
        					else:
        						Strain = refdict[smallSP]
        						forcsv[Strain] = thisSP 
        		for item in header:
        			if item not in forcsv.keys():
        				forcsv[item] = " " 
        		writer.writerow(forcsv)
        				    
        
        
        
        
if __name__ == '__main__':
        main()
