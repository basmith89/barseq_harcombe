#!/usr/bin/python

import csv
import argparse, textwrap


parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
    description = textwrap.dedent('''\
                                 Author: Brian A. Smith
                                 University of Minnesota
                                 smithba@umn.edu
                                 
                 Tie Barseq information to gene information
                 The BarSeq ratio file input file should be formatted and tab delimitted as follows:
                 barcode position ratio1 ratio2 ...

                 The gene info file is output by the Feba Barseq code found at:
                 https://bitbucket.org/berkeleylab/feba/src/master/'''))

parser.add_argument("-b", "--barseq", required = True,
           help = "Barseq ratio File")
parser.add_argument("-g", "--gene_info", required = True,
           help = "File containing gene and position information")
parser.add_argument("-o", "--output", required = True,
           help = "Desired output file name")


args = parser.parse_args()

out_file = open(args.output, 'w')

prev_end = 0
counter = 0

with open(args.barseq, 'r') as f:
    readcsv_bar = csv.reader(f, delimiter='\t')
    headers = next(readcsv_bar)
    headerlen = len(headers)
    ratio_amt = headerlen - 2
    #print(headers, headers[-ratio_amt:])
    out_file.write(str(headers).strip("[]").replace(',', '\t').replace("'", "") + "\tstart" + "\tend" + "\tgene_desc\n")
    for row in readcsv_bar:
        ###Note this will need to be changed based on how many ratio value are used
        ###The data var will need to be adjusted as well
        ###Could think of a way to automate this or make it user defined
        barcode, pos, ratio1, ratio2, ratio3, ratio4 = row[0], row[1], row[2], row[3], row[4], row[5]
        #print(barcode, pos, ratio1, ratio2)


        with open(args.gene_info, 'r') as gene:
            readcsv_gene = csv.reader(gene, delimiter='\t')
            #next(readcsv_bar) #skip header line
            for gene_row in readcsv_gene:
                curr_start = gene_row[4]
                curr_end   = gene_row[5]
                gene_desc = gene_row[8]
                
                if curr_start == "begin" and curr_end == "end":  #skip header
                    continue
                #print("current position: " + str(pos))
                #print(curr_start, curr_end)

                if int(pos) > int(curr_start) and int(pos) < int(curr_end): #not sure why but sometimes values were not ints
                    #print(gene_desc)
                    #print(curr_start, curr_end, prev_end)
                    data = barcode, pos, ratio1, ratio2, ratio3, ratio4, curr_start, curr_end, gene_desc
                    out_file.write('\t'.join(data) + '\n')
                    #print("Debug" + barcode + " "+ str(curr_start) + " "+ str(curr_end) + " "+ str(prev_end))
                    break
                
                if int(prev_end) < int(pos) < int(curr_start):
                    #print("Intergenic insertion")
                    #print(curr_start, curr_end, prev_end)
                    data = barcode, pos, ratio1, ratio2, ratio3, ratio4, curr_start, curr_end + "\tIntergeneic insertion"
                    out_file.write('\t'.join(data) + '\n')
                    #print(barcode)
                prev_end = curr_end

print("Done!  Data are stored in: " + args.output)
