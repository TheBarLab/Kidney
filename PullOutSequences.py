
# Genome assembly hg19/ GRCh37.p13
#
from Bio import Entrez, SeqIO
import random
import pandas as pd
import math


def get_seq(ids,startloc,endloc,strand ):
    Entrez.email = "A.N.Other@example.com"     # Always tell NCBI who you are
    handle = Entrez.efetch(db="nucleotide",
                           id=ids, # id="CM000663.1",#id="307603377",
                           rettype="fasta",
                           strand=strand, # 1= plus strand 2= minus strand
                           seq_start=startloc,# 226077128,
                           seq_stop=endloc)# 226078128)
    record = SeqIO.read(handle, "fasta")
    handle.close()
    print(record.seq)
    return record.seq

# create GenBank IDs
startfrom=663 #1...22,x,y
GenBankID=["na"]
for i in range(0,24):
    GenBankID.append("CM000"+str(startfrom+i)+".1")

print(GenBankID)

filedir='/home/naorsagy/Desktop/R_Methyl/DNN_files_and_scripts/'
fileout='seqs_for_cg_sites.csv'

# pull out DNA sequences and save to file
if True:
    # read csv with probe sites list and genomic coordinates
    corrdfile="/home/naorsagy/Desktop/R_Methyl/DNN_files_and_scripts/files/HM450_hg19.csv"
    df = pd.read_csv(corrdfile)

    # header: CpG_chrm	CpG_beg	CpG_end	probe_strand	probeID	probeType	orientation	probeCpGcnt	context35	probeBeg	probeEnd	length	gene	gene_HGNC	chrm_A
    # probe_strand : '+' or '-'
    # CpG_beg: 53468111
    # probeID: 'cg2343242'
    # CpG_chrm: chromosome number 1....21 ,x, y
    #print(df)

    sites=[]
    coordinates=[]
    chromosome_num=[]
    IDS=[]
    Sequences=[]
    strands=[]
    findx=10
    skip_to_index=450001
    #indx=1
    for indx,row in df.iterrows():
        #indx=indx+1
        #print(row)
        # skip indx
        if indx<skip_to_index:
            continue
        chromosome =row['chrom_num']
        if isinstance(chromosome, str): # because on rare occasions it won't be a string
            if 'M' in chromosome: # if mitochondrial
                print('skpping mitochondrial chromosome')
                continue
        if '+' in row['probe_strand']:
            strand=1
        else:
            strand=2
        startloc=int(row['CpG_beg'])-500
        endloc=int(row['CpG_beg'])+500
        try:
            print("site #"+str(indx+1)+"  ; "+row['probeID']+" : getting sequence for chromosome "+str(chromosome)+" strand direction: "+row['probe_strand']+ " ; around " +str(row['CpG_beg'])+ " DB ID:  "+ GenBankID[int(chromosome)] + "; from loc "+str(startloc)+" to "+str(endloc))
            cur_seq = get_seq(GenBankID[int(chromosome)], startloc, endloc, strand)
        except:
            print("error chromosome " + str(chromosome) )
        print(cur_seq)
        #order is probe,coordinates,CHROMOSOME,IDS,SEQ
        sites.append(row["probeID"])
        strands.append(row["probe_strand"])
        coordinates.append(row["CpG_beg"])
        chromosome_num.append(chromosome)
        IDS.append(GenBankID[int(chromosome)])
        Sequences.append(cur_seq)

        if (indx%50000 == 0 and indx!=0) or indx==df.index[-1]:
            print('writing file '+filedir+str(findx)+'_'+fileout)
            datad = {'sites': sites, 'cooridnates' : coordinates, 'strand':strands,
                    'chromosome_num' : chromosome_num, 'IDS': IDS,
                    'Sequence:': Sequences}

            data=pd.DataFrame(datad)
            print(data)
            #export data to csv
            data.to_csv(filedir+str(findx)+'_'+fileout)
            sites = []
            coordinates = []
            chromosome_num = []
            IDS = []
            Sequences = []
            strands = []
            findx=findx+1
