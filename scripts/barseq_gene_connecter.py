#!/usr/bin/python

import csv
import argparse, textwrap
import time

Tstart = time.process_time()

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
group = parser.add_mutually_exclusive_group()
group.add_argument("-T", "--tnseq_diff", required = False,
           action = 'store_true', help = "perform the tnseq_diff option")
group.add_argument("-R", "--ratios", required = False,
           action = 'store_true', help = "perform the ratio analysis option")

args = parser.parse_args()

out_file = open(args.output, 'w')

prev_end = 0
counter = 0


if args.ratios:
    with open(args.barseq, 'r') as f:
        readcsv_bar = csv.reader(f, delimiter='\t')
        headers = next(readcsv_bar)
        headerlen = len(headers)
        ratio_amt = headerlen - 2
        #print(headers, headers[-ratio_amt:])
        out_file.write(str(headers).strip("[]").replace(',', '\t').replace("'", "") + "\tstart" + "\tend" + "\tgene_name" "\tgene_desc\n")
        
        
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
                    gene_name  = gene_row[7]
                    gene_desc  = gene_row[8]
                    
                    
                    if curr_start == "begin" and curr_end == "end":  #skip header
                        continue
                    #print("current position: " + str(pos))
                    #print(curr_start, curr_end)

                    if int(pos) > int(curr_start) and int(pos) < int(curr_end): #not sure why but sometimes values were not ints
                        #print(gene_desc)
                        #print(curr_start, curr_end, prev_end)
                        data = barcode, pos, ratio1, ratio2, ratio3, ratio4, curr_start, curr_end, gene_name, gene_desc
                        out_file.write('\t'.join(data) + '\n')
                        #print("Debug" + barcode + " "+ str(curr_start) + " "+ str(curr_end) + " "+ str(prev_end))
                        break
                    
                    if int(prev_end) < int(pos) < int(curr_start):
                        #print("Intergenic insertion")
                        #print(curr_start, curr_end, prev_end)
                        data = barcode, pos, ratio1, ratio2, ratio3, ratio4, curr_start, curr_end + "\t\tIntergenic insertion"
                        out_file.write('\t'.join(data) + '\n')
                        #print(barcode)
                    prev_end = curr_end

    
elif args.tnseq_diff:
    with open(args.barseq, 'r') as f:
        readcsv_bar = csv.reader(f, delimiter='\t')
        headers = next(readcsv_bar)
        headerlen = len(headers)
        ratio_amt = headerlen - 2
        #print(headers, headers[-ratio_amt:])
        out_file.write("GeneID" + "\tLength\t" + str(headers).strip("[]").replace(',', '\t').replace("'", "") + "\n")


        for row in readcsv_bar:
            pos, count1, count2, count3, count4, count5, count6, count7, count8, count9, count10, count11, count12, count13, count14, count15, count16, count17, count18, count19, count20, count21, count22, count23, count24, count25 = row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9], row[10], row[11], row[12], row[13], row[14], row[15], row[16], row[17], row[18], row[19], row[20], row[21], row[22], row[23], row[24], row[25]

            with open(args.gene_info, 'r') as gene:
                readcsv_gene = csv.reader(gene, delimiter='\t')
                #next(readcsv_bar) #skip header line
                for gene_row in readcsv_gene:
                    gene_id    = gene_row[0]
                    start      = gene_row[4]
                    end        = gene_row[5]
                    gene_name  = gene_row[7]
                    gene_desc  = gene_row[8]
                    
                    
                    if start == "begin" and end == "end":  #skip header
                        continue
                    #print("current position: " + str(pos))
                    #print(curr_start, curr_end)

                    if int(pos) > int(start) and int(pos) < int(end): #not sure why but sometimes values were not ints
                        #print(gene_desc)
                        #print(curr_start, curr_end, prev_end)
                        gene_length = int(end) - int(start)
                        data = gene_id, str(gene_length), pos, count1, count2, count3, count4, count5, count6, count7, count8, count9, count10, count11, count12, count13, count14, count15, count16, count17, count18, count19, count20, count21, count22, count23, count24, count25
                        out_file.write('\t'.join(data) + '\n')
                        #print("Debug" + barcode + " "+ str(curr_start) + " "+ str(curr_end) + " "+ str(prev_end))
                    
runT = time.process_time() - Tstart
print("Done! Runtime = " + str(runT) + "\nData are stored in: " + args.output)
