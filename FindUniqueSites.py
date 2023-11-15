import sys
import pandas as pd
import numpy as np
from statistics import mean
import csv
import re
from Bio import Entrez, SeqIO
from os.path import exists

# read huge cg data once
def read_huge_cg_data():
    cg_loc_file = "/home/naorsagy/Desktop/R_Methyl/HM450.hg19_min_data.xlsx"
    print("reading huge cg locations file....")
    cg_data = pd.read_excel(cg_loc_file)

    return cg_data

# pull out coordinates for cg site
def get_coordinates_for_cg_site(cg_site, cg_data):
   for indx,item in enumerate(cg_data["probeID"]):
       if item==cg_site:
           return cg_data["CpG_chrm"][indx], cg_data["CpG_beg"][indx]

   return 0,0 #  is not found

def indexes(iterable, obj):
    return (index for index, elem in enumerate(iterable) if elem == obj)

def rest_indexes(iterable, obj):
    return (index for index, elem in enumerate(iterable) if elem != obj)
# find unique sites and filter for 2000 NaN - all tissues

# sort by tissue  - input is the main file
if False:
    tissues_list = ["whole_blood", "brain", "adipose", "colon", "skeletal_muscle",
                    "fallopian_tube", "pancreas", "stomach", "prostate", "thyroid",
                    "esophagus", "placenta", "adipose_-_visceral", "skin", "nasal_epithelium",
                    "testis", "bone_marrow", "cartilage", "bone", "adrenal_gland",
                    "tongue", "cord_blood", "bladder", "nasopharynx", "spleen",
                    "saliva", "liver", "breast", "lung", "kidney"]
    #infile = '/home/naorsagy/Desktop/R_Methyl/unique_cg_by_site_filtered_2000Vals/20220913_unique_cg_all_tissues_2kvals.csv'
    infile= '/home/naorsagy/Desktop/R_Methyl/unique_cg_by_site_filtered_more2000Vals/20221123_unique_cg_all_tissues_2kvals_0_1_thr.csv'
    outfolder='/home/naorsagy/Desktop/R_Methyl/unique_cg_by_site_filtered_more2000Vals/'

    for tissue in tissues_list:
        print("finding " +tissue+".....")
        outfile=outfolder+'20221126_unique_cg_'+tissue+'.csv'

        fo = open(outfile, 'w')
        un_count=0
        head_flag=0
        foundcount=0
        with open(infile) as fp:  # scan methylation file
            while True:
                cur_line = fp.readline()
                if not cur_line:
                    break
                #print(cur_line)
                cur_line=cur_line.replace('[','')
                cur_line = cur_line.replace('(', '')
                cur_line=cur_line.replace(']', '')
                cur_line=cur_line.replace("'", '')
                line_splt = cur_line.split(",")
                #print(line_splt[0])
                #print(line_splt[0])
                if head_flag==0:
                    head_flag=1
                    fo.write(cur_line)
                if line_splt[0]==tissue:
                    un_count=un_count+1
                    fo.write(cur_line)
        if un_count==0:
            fo.write("EMPTY")
            print("NONE")
        else:
            print("found "+str(un_count)+" unique sites for "+tissue)
        fo.close()

#find unique sites for all tissues - create main file
if False:
    difference_threshold=0.1

    infile = '/home/naorsagy/Desktop/R_Methyl/tissue_methylation_v1.txt'
    outfile= '/home/naorsagy/Desktop/R_Methyl/unique_cg_by_site_filtered_more2000Vals/20221123_unique_cg_all_tissues_2kvals_0_1_thr.csv'
    outfol = "/home/naorsagy/Desktop/R_Methyl/unique_cg_by_site_filtered_more2000Vals/"
    tissues_list = ["whole_blood", "brain", "adipose", "colon", "skeletal_muscle",
                    "fallopian_tube", "pancreas", "stomach", "prostate", "thyroid",
                    "esophagus", "placenta", "adipose_-_visceral", "skin", "nasal_epithelium",
                    "testis", "bone_marrow", "cartilage", "bone", "adrenal_gland",
                    "tongue", "cord_blood", "bladder", "nasopharynx", "spleen",
                    "saliva", "liver", "breast", "lung", "kidney"]

    breakcount=0
    total_skipped=0
    breaknow=0
    fo = open(outfile, 'w')
    with open(infile) as fp: # scan methylation file
            while True:
                cur_line = fp.readline()
                if not cur_line:
                    break
                breakcount=breakcount+1
                #if breakcount==260000:
                #    break
                if breakcount==2: # if we're in the tissue's list line
                    tissues_line=cur_line.split("\t")[1:]
                    header_line=cur_line.split("\t")
                    header_line.insert(0,'tissue')
                    print(tissues_line)
                    fo.write(f"{header_line}\n")
                    continue
                elif breakcount<2:
                     continue

                cur_line=cur_line.split("\t")

                if (breakcount % 250) == 0:
                    print(cur_line[0] +": analyzing...." +str(breakcount))

                #if breakcount<248630 and breakcount>248620:
                #    print("***********************************"+cur_line[0])
                #if cur_line[0]!="cg16667631":
                    #print(cur_line[0])
                #    continue
                #else:
                #    print("*********FOUND**************************" + cur_line[0])
                #    breaknow=1
                cur_line = [float('nan') if element == 'NA' else element for element in cur_line]
                cur_line = [float('nan') if element == 'NA\n' else element for element in cur_line]
                cur_line_nums=cur_line[1:] # omit cg site name



                cur_line_nums = [float(s) if isinstance(s, str) else s for s in cur_line_nums]
                #print(cur_line_nums)
                count_nans=np.sum(np.isnan(cur_line_nums))
                #print(len(cur_line_nums))
                if count_nans>3323: # skip sites with less than 2000 values
                    print("skipping "+cur_line[0]+ ", nan number:"+str(count_nans))
                    #print(cur_line)
                    total_skipped=total_skipped+1
                    #print(cur_line_nums)
                    continue

                #find tissues locs in list
                for sel_tissue in tissues_list:
                    #cur_line_help=[]
                    cur_tissue_indices=list(indexes(tissues_line,sel_tissue))
                    cur_tissue_rest=list(rest_indexes(tissues_line,sel_tissue))
                    res_list = [cur_line_nums[i] for i in cur_tissue_indices]
                    rest_list= [cur_line_nums[i] for i in cur_tissue_rest]
                    tissue_mean=np.nanmean(res_list)  # mean for selected tissue
                    rest_mean=np.nanmean(rest_list)
                    botfive = np.nanquantile(rest_list,0.05)  # quantile(rest_of_tissues,0.05)
                    topfive = np.nanquantile(rest_list, 0.95)


                    if (tissue_mean < (botfive-difference_threshold)) or (tissue_mean > (topfive+difference_threshold)):
                        print(cur_line[0] +","+str(breakcount)+": unique site - " + sel_tissue + " mean/bot/top: " + str(tissue_mean) + "/" + str(
                            botfive) + "/" + str(topfive))
                       # all_new=cur_line
                       # print(cur_line)
                       # all_new.insert(0,sel_tissue)
                       # print(cur_line)
                        fo.write(f"{sel_tissue,cur_line}\n")

                        #if breaknow == 1:
                        #    print(cur_line[0] + "," + str(breakcount) +
                        #      ": unique site - " + sel_tissue + " mean/bot/top: " + str(tissue_mean) + "/" + str(botfive) + "/" + str(topfive))



                        # write to file

                #if breaknow == 1:
                    #print(cur_line)
                    #print("out")
                    #print(cur_line_help)
                 #   break
                    #if ((kidney_mean < (botfive - difference_threshold)) | | (
                    #        kidney_mean > (topfive + difference_threshold)))  # if unique
                    #print(res_list)

    print(" total methylaytion sites skipped:" +str(total_skipped))
    fo.close()

# create tissue and cg sites table - one file with no methylation values
if False:
    infol = "/home/naorsagy/Desktop/R_Methyl/unique_cg_by_site_filtered_more2000Vals/"
    fileout= infol+"20221126_unique_cg_site_by_tissue.xlsx"
    tissues_list = ["whole_blood", "brain", "adipose", "colon", "skeletal_muscle",
                    "fallopian_tube", "pancreas", "stomach", "prostate", "thyroid",
                    "esophagus", "placenta", "adipose_-_visceral", "skin", "nasal_epithelium",
                    "testis", "bone_marrow", "cartilage", "bone", "adrenal_gland",
                    "tongue", "cord_blood", "bladder", "nasopharynx", "spleen",
                    "saliva", "liver", "breast", "lung", "kidney",
                    "sperm"]

    fo = open(fileout, 'w')
    for tissue in tissues_list:
        print("\n***********\n"+tissue+"\n*********")
        filein = infol + '20221126_unique_cg_'+tissue+'.csv'
        #'20221126_unique_cg_'+tissue+'.csv'
        if (not exists(filein)): # if file doesn't exist
            print("no input file for " + tissue)
            continue
        cur_list=[]
        cur_list.append(tissue)
        with open(filein) as fp:
            while True:
                cur_line = fp.readline()
                if not cur_line:
                    print(tissue+" reached eof")
                    break
                cur_line=cur_line.split(" ")

                if len(cur_line[0])>24:
                    #print(cur_line[0])
                    cur_line_drl=cur_line[0].split(",")
                    #cg_name=cur_line_drl[0]
                    #print(cg_name)
                    cg_name=cur_line_drl[1]
                else:
                    cg_name=cur_line[0]
                #print(cg_name)
                cur_list.append(cg_name)
        my_string = ','.join(cur_list)
        print(my_string)
        fo.write(my_string+"\n")

    fo.close()

if False:
    infol="/home/naorsagy/Desktop/R_Methyl/unique_cg_by_tissue/"
    outfol="/home/naorsagy/Desktop/R_Methyl/unique_cg_by_site_filtered_2000NaN/"
    tissues_list=["whole_blood",        "brain",              "adipose",            "colon",              "skeletal_muscle",
                   "fallopian_tube",     "pancreas",           "stomach",            "prostate",           "thyroid",
                   "esophagus",          "placenta" ,          "adipose_-_visceral", "skin",               "nasal_epithelium",
                   "testis",             "bone_marrow",        "cartilage",          "bone",               "adrenal_gland",
                   "tongue" ,            "cord_blood"  ,       "bladder",            "nasopharynx",        "spleen",
                   "saliva"  ,           "liver"  ,            "breast",             "lung",               "kidney",
                   "sperm" ]

    #cg_data=read_huge_cg_data()


    for tissue in tissues_list:

        filein=infol+"20220907_"+tissue+"_unique_cg_sites.csv"
        fileout=outfol+"20220907_"+tissue+"_unique_cg_sites_fil.csv"

        # check if file exists
        if (not exists(filein)):
            print("no input file for "+tissue)
            continue

        fo=open(fileout, 'w')
        print("opening "+filein)

        with open(filein) as fp:
            while True:
                cur_line = fp.readline()
                if not cur_line:
                    break
                #print(cur_line)
                # find how many NaNs in line
                inc_cou=cur_line.count("NaN")
                #print(inc_cou)
                if inc_cou<2000:
                    fo.write(cur_line)
                #else:
                #    print("skipping line")

        fo.close()
        #for tissue in tissues_list:


#find unique sites for all tissues - create main file
if False:
    difference_threshold=0.1

    infile = '/home/naorsagy/Desktop/R_Methyl/tissue_methylation_v1.txt'
    outfile= '/home/naorsagy/Desktop/R_Methyl/unique_cg_by_site_filtered_more2000Vals/20230116_unique_cg_all_tissues_2kvals_0_1_thr_with_stats.csv'
    outfol = "/home/naorsagy/Desktop/R_Methyl/unique_cg_by_site_filtered_more2000Vals/"
    tissues_list = ["whole_blood", "brain", "adipose", "colon", "skeletal_muscle",
                    "fallopian_tube", "pancreas", "stomach", "prostate", "thyroid",
                    "esophagus", "placenta", "adipose_-_visceral", "skin", "nasal_epithelium",
                    "testis", "bone_marrow", "cartilage", "bone", "adrenal_gland",
                    "tongue", "cord_blood", "bladder", "nasopharynx", "spleen",
                    "saliva", "liver", "breast", "lung", "kidney"]

    breakcount=0
    total_skipped=0
    breaknow=0
    #fo = open(outfile, 'w')

    keepsite=[]
    keeptissue=[]
    keeptissuemean=[]
    keeprest5perquant=[]
    keeprest95perquant=[]
    keeploworhigh=[]

    with open(infile) as fp: # scan methylation file
            while True:
                cur_line = fp.readline()
                if not cur_line:
                    break
                breakcount=breakcount+1
                #if breakcount==260000:
                #    break
                if breakcount==2: # if we're in the tissue's list line
                    tissues_line=cur_line.split("\t")[1:]
                    header_line=cur_line.split("\t")
                    header_line.insert(0,'tissue')
                    print(tissues_line)
                    #fo.write(f"{header_line}\n")
                    continue
                elif breakcount<2:
                     continue

                cur_line=cur_line.split("\t")

                if (breakcount % 250) == 0:
                    print(cur_line[0] +": analyzing...." +str(breakcount))

                #if breakcount>10000:
                #    break
                #if breakcount<248630 and breakcount>248620:
                #    print("***********************************"+cur_line[0])
                #if cur_line[0]!="cg16667631":
                    #print(cur_line[0])
                #    continue
                #else:
                #    print("*********FOUND**************************" + cur_line[0])
                #    breaknow=1
                cur_line = [float('nan') if element == 'NA' else element for element in cur_line]
                cur_line = [float('nan') if element == 'NA\n' else element for element in cur_line]
                cur_line_nums=cur_line[1:] # omit cg site name



                cur_line_nums = [float(s) if isinstance(s, str) else s for s in cur_line_nums]
                #print(cur_line_nums)
                count_nans=np.sum(np.isnan(cur_line_nums))
                #print(len(cur_line_nums))
                if count_nans>3323: # skip sites with less than 2000 values
                    print("skipping "+cur_line[0]+ ", nan number:"+str(count_nans))
                    #print(cur_line)
                    total_skipped=total_skipped+1
                    #print(cur_line_nums)
                    continue

                #find tissues locs in list
                for sel_tissue in tissues_list:
                    #cur_line_help=[]
                    cur_tissue_indices=list(indexes(tissues_line,sel_tissue))
                    cur_tissue_rest=list(rest_indexes(tissues_line,sel_tissue))
                    res_list = [cur_line_nums[i] for i in cur_tissue_indices]
                    rest_list= [cur_line_nums[i] for i in cur_tissue_rest]
                    tissue_mean=np.nanmean(res_list)  # mean for selected tissue
                    rest_mean=np.nanmean(rest_list)
                    botfive = np.nanquantile(rest_list,0.05)  # quantile(rest_of_tissues,0.05)
                    topfive = np.nanquantile(rest_list, 0.95)


                    if (tissue_mean < (botfive-difference_threshold)) or (tissue_mean > (topfive+difference_threshold)):
                        print(cur_line[0] +","+str(breakcount)+": unique site - " + sel_tissue + " mean/bot/top: " + str(tissue_mean) + "/" + str(
                            botfive) + "/" + str(topfive))
                       # all_new=cur_line
                       # print(cur_line)
                       # all_new.insert(0,sel_tissue)
                       # print(cur_line)
                        #fo.write(f"{sel_tissue,cur_line}\n")
                        keepsite.append(cur_line[0])
                        keeptissue.append(sel_tissue)
                        keeptissuemean.append(tissue_mean)
                        keeprest5perquant.append(botfive)
                        keeprest95perquant.append(topfive)
                        if (tissue_mean < (botfive-difference_threshold)):
                            keeploworhigh.append("low")
                        else:
                            keeploworhigh.append("high")

                        #if breaknow == 1:
                        #    print(cur_line[0] + "," + str(breakcount) +
                        #      ": unique site - " + sel_tissue + " mean/bot/top: " + str(tissue_mean) + "/" + str(botfive) + "/" + str(topfive))



                        # write to file

                #if breaknow == 1:
                    #print(cur_line)
                    #print("out")
                    #print(cur_line_help)
                 #   break
                    #if ((kidney_mean < (botfive - difference_threshold)) | | (
                    #        kidney_mean > (topfive + difference_threshold)))  # if unique
                    #print(res_list)

    print(" total methylaytion sites skipped:" +str(total_skipped))
    #fo.close()
    #output_dict=dict()
    output_dict={'site':keepsite,'unique tissue':keeptissue,'tissue mean':keeptissuemean,'0.05 quantile, rest of tissues':keeprest5perquant
        , '0.95 quantile, rest of tissues': keeprest95perquant,'unique low/high':keeploworhigh}
    print(output_dict)
    data = pd.DataFrame(output_dict)
    # export data to csv
    data.to_csv(outfile)

    # find unique sites for all tissues - create main file
    if True:
        difference_threshold = 0.1

        infile = '/home/naorsagy/Desktop/R_Methyl/tissue_methylation_v1.txt'
        outfile = '/home/naorsagy/Desktop/R_Methyl/unique_cg_by_site_filtered_more2000Vals/20230116_unique_cg_all_tissues_2kvals_0_1_thr_with_stats.csv'
        outfol = "/home/naorsagy/Desktop/R_Methyl/unique_cg_by_site_filtered_more2000Vals/"
        tissues_list = ["whole_blood", "brain", "adipose", "colon", "skeletal_muscle",
                        "fallopian_tube", "pancreas", "stomach", "prostate", "thyroid",
                        "esophagus", "placenta", "adipose_-_visceral", "skin", "nasal_epithelium",
                        "testis", "bone_marrow", "cartilage", "bone", "adrenal_gland",
                        "tongue", "cord_blood", "bladder", "nasopharynx", "spleen",
                        "saliva", "liver", "breast", "lung", "kidney"]

        breakcount = 0
        total_skipped = 0
        breaknow = 0
        # fo = open(outfile, 'w')

        keepsite = []
        keeptissue = []
        keeptissuemean = []
        keeprest5perquant = []
        keeprest95perquant = []
        keeploworhigh = []

        with open(infile) as fp:  # scan methylation file
            while True:
                cur_line = fp.readline()
                if not cur_line:
                    break
                breakcount = breakcount + 1
                # if breakcount==260000:
                #    break
                if breakcount == 2:  # if we're in the tissue's list line
                    tissues_line = cur_line.split("\t")[1:]
                    header_line = cur_line.split("\t")
                    header_line.insert(0, 'tissue')
                    print(tissues_line)
                    # fo.write(f"{header_line}\n")
                    continue
                elif breakcount < 2:
                    continue

                cur_line = cur_line.split("\t")

                if (breakcount % 250) == 0:
                    print(cur_line[0] + ": analyzing...." + str(breakcount))

                # if breakcount>10000:
                #    break
                # if breakcount<248630 and breakcount>248620:
                #    print("***********************************"+cur_line[0])
                # if cur_line[0]!="cg16667631":
                # print(cur_line[0])
                #    continue
                # else:
                #    print("*********FOUND**************************" + cur_line[0])
                #    breaknow=1
                cur_line = [float('nan') if element == 'NA' else element for element in cur_line]
                cur_line = [float('nan') if element == 'NA\n' else element for element in cur_line]
                cur_line_nums = cur_line[1:]  # omit cg site name

                cur_line_nums = [float(s) if isinstance(s, str) else s for s in cur_line_nums]
                # print(cur_line_nums)
                count_nans = np.sum(np.isnan(cur_line_nums))
                # print(len(cur_line_nums))
                if count_nans > 3323:  # skip sites with less than 2000 values
                    print("skipping " + cur_line[0] + ", nan number:" + str(count_nans))
                    # print(cur_line)
                    total_skipped = total_skipped + 1
                    # print(cur_line_nums)
                    continue

                # find tissues locs in list
                for sel_tissue in tissues_list:
                    # cur_line_help=[]
                    cur_tissue_indices = list(indexes(tissues_line, sel_tissue))
                    cur_tissue_rest = list(rest_indexes(tissues_line, sel_tissue))
                    res_list = [cur_line_nums[i] for i in cur_tissue_indices]
                    rest_list = [cur_line_nums[i] for i in cur_tissue_rest]
                    tissue_mean = np.nanmean(res_list)  # mean for selected tissue
                    rest_mean = np.nanmean(rest_list)
                    botfive = np.nanquantile(rest_list, 0.05)  # quantile(rest_of_tissues,0.05)
                    topfive = np.nanquantile(rest_list, 0.95)

                    if (tissue_mean < (botfive - difference_threshold)) or (
                            tissue_mean > (topfive + difference_threshold)):
                        print(cur_line[0] + "," + str(
                            breakcount) + ": unique site - " + sel_tissue + " mean/bot/top: " + str(
                            tissue_mean) + "/" + str(
                            botfive) + "/" + str(topfive))
                        # all_new=cur_line
                        # print(cur_line)
                        # all_new.insert(0,sel_tissue)
                        # print(cur_line)
                        # fo.write(f"{sel_tissue,cur_line}\n")
                        keepsite.append(cur_line[0])
                        keeptissue.append(sel_tissue)
                        keeptissuemean.append(tissue_mean)
                        keeprest5perquant.append(botfive)
                        keeprest95perquant.append(topfive)
                        if (tissue_mean < (botfive - difference_threshold)):
                            keeploworhigh.append("low")
                        else:
                            keeploworhigh.append("high")

                        # if breaknow == 1:
                        #    print(cur_line[0] + "," + str(breakcount) +
                        #      ": unique site - " + sel_tissue + " mean/bot/top: " + str(tissue_mean) + "/" + str(botfive) + "/" + str(topfive))

                        # write to file

                # if breaknow == 1:
                # print(cur_line)
                # print("out")
                # print(cur_line_help)
                #   break
                # if ((kidney_mean < (botfive - difference_threshold)) | | (
                #        kidney_mean > (topfive + difference_threshold)))  # if unique
                # print(res_list)

        print(" total methylaytion sites skipped:" + str(total_skipped))
        # fo.close()
        # output_dict=dict()
        output_dict = {'site': keepsite, 'unique tissue': keeptissue, 'tissue mean': keeptissuemean,
                       '0.05 quantile, rest of tissues': keeprest5perquant
            , '0.95 quantile, rest of tissues': keeprest95perquant, 'unique low/high': keeploworhigh}
        print(output_dict)
        data = pd.DataFrame(output_dict)
        # export data to csv
        data.to_csv(outfile)

#find unique sites for buccal epithelium
if False:
    difference_threshold=0.1

    infile = '/home/naorsagy/Desktop/R_Methyl/tissue_methylation_v1.txt'
    buccalfile='/home/naorsagy/Desktop/R_Methyl/BuccalEpithelium/GSE94876_meth_data.csv'
    outfile= '/home/naorsagy/Desktop/R_Methyl/unique_cg_by_site_filtered_more2000Vals/20230130_unique_cg_buccalepi_0_1_with_stats.csv'
    outfol = "/home/naorsagy/Desktop/R_Methyl/unique_cg_by_site_filtered_more2000Vals/"
    tissues_list = ["whole_blood", "brain", "adipose", "colon", "skeletal_muscle",
                    "fallopian_tube", "pancreas", "stomach", "prostate", "thyroid",
                    "esophagus", "placenta", "adipose_-_visceral", "skin", "nasal_epithelium",
                    "testis", "bone_marrow", "cartilage", "bone", "adrenal_gland",
                    "tongue", "cord_blood", "bladder", "nasopharynx", "spleen",
                    "saliva", "liver", "breast", "lung", "kidney"]

    breakcount=0
    total_skipped=0
    breaknow=0
    #fo = open(outfile, 'w')

    keepsite=[]
    keeptissue=[]
    keeptissuemean=[]
    keeprest5perquant=[]
    keeprest95perquant=[]
    keeploworhigh=[]

    # read meth file
    print('read buccal methylation data...')
    buccal_data=pd.read_csv(buccalfile)


    with open(infile) as fp: # scan methylation file
            while True:
                cur_line = fp.readline()
                if not cur_line:
                    break
                breakcount=breakcount+1
                #if breakcount==260000:
                #    break
                if breakcount==2: # if we're in the tissue's list line
                    tissues_line=cur_line.split("\t")[1:]
                    header_line=cur_line.split("\t")
                    header_line.insert(0,'tissue')
                    print(tissues_line)
                    #fo.write(f"{header_line}\n")
                    continue
                elif breakcount<2:
                     continue

                cur_line=cur_line.split("\t")

                if (breakcount % 250) == 0:
                    print(cur_line[0] +": analyzing...." +str(breakcount))

                #if breakcount>10000:
                #    break
                #if breakcount<248630 and breakcount>248620:
                #    print("***********************************"+cur_line[0])
                #if cur_line[0]!="cg16667631":
                    #print(cur_line[0])
                #    continue
                #else:
                #    print("*********FOUND**************************" + cur_line[0])
                #    breaknow=1
                cur_line = [float('nan') if element == 'NA' else element for element in cur_line]
                cur_line = [float('nan') if element == 'NA\n' else element for element in cur_line]
                cur_line_nums=cur_line[1:] # omit cg site name

                cur_line_nums = [float(s) if isinstance(s, str) else s for s in cur_line_nums]
                #print(cur_line_nums)
                count_nans=np.sum(np.isnan(cur_line_nums))
                #print(len(cur_line_nums))
                if count_nans>3323: # skip sites with less than 2000 values
                    print("skipping "+cur_line[0]+ ", nan number:"+str(count_nans))
                    #print(cur_line)
                    total_skipped=total_skipped+1
                    #print(cur_line_nums)
                    continue

                #find tissues locs in list
                res_list = buccal_data.loc[buccal_data['site'] == cur_line[0]].values
                res_list=res_list[0][1:] # drop site name
                #print(res_list)
                tissue_mean = np.nanmean(res_list)  # mean for selected tissue
                if np.isnan(tissue_mean):
                    print("could not calculate mean for buccal site "+cur_line[0])
                rest_list= cur_line_nums
                rest_mean=np.nanmean(rest_list)
                botfive = np.nanquantile(rest_list,0.05)  # quantile(rest_of_tissues,0.05)
                topfive = np.nanquantile(rest_list, 0.95)


                if (tissue_mean < (botfive-difference_threshold)) or (tissue_mean > (topfive+difference_threshold)):
                    print(cur_line[0] +","+str(breakcount)+": unique site for buccal epi, mean/bot/top: " + str(tissue_mean) + "/" + str(
                        botfive) + "/" + str(topfive))
                       # all_new=cur_line
                       # print(cur_line)
                       # all_new.insert(0,sel_tissue)
                       # print(cur_line)
                        #fo.write(f"{sel_tissue,cur_line}\n")
                    keepsite.append(cur_line[0])
                    keeptissue.append('Buccal')
                    keeptissuemean.append(tissue_mean)
                    keeprest5perquant.append(botfive)
                    keeprest95perquant.append(topfive)
                    if (tissue_mean < (botfive-difference_threshold)):
                        keeploworhigh.append("low")
                    else:
                        keeploworhigh.append("high")



    print(" total methylaytion sites skipped:" +str(total_skipped))
    #fo.close()
    #output_dict=dict()
    output_dict={'site':keepsite,'unique tissue':keeptissue,'tissue mean':keeptissuemean,'0.05 quantile, rest of tissues':keeprest5perquant
        , '0.95 quantile, rest of tissues': keeprest95perquant,'unique low/high':keeploworhigh}
    print(output_dict)
    data = pd.DataFrame(output_dict)
    # export data to csv
    data.to_csv(outfile)


def unique(list1):
    # initialize a null list
    unique_list = []

    # traverse for all elements
    for x in list1:
        # check if exists in unique_list or not
        if x not in unique_list:
            unique_list.append(x)
    # print list
    return unique_list

# blood analysis
if True:
    difference_threshold=0.1

    infile = '/run/user/1000/gvfs/sftp:host=132.66.247.107/Shared_BarLab/NGDC-CNCB/blood_methylation_v1/blood_methylation_v1_ns.txt'
    outfile= '/home/naorsagy/Desktop/R_Methyl/blood_unique_cg_by_site_filtered_more2000Vals/20230418_blood_unique_cg_2kvals_0_1_thr_with_stats.csv'
    tissuesfile='/run/user/1000/gvfs/sftp:host=132.66.247.107/Shared_BarLab/NGDC-CNCB/blood_methylation_v1/tissues_list.txt'
    outfol = "/home/naorsagy/Desktop/R_Methyl/blood_unique_cg_by_site_filtered_more2000Vals/"

    file = open(tissuesfile, "r")
    tissues_list = list(csv.reader(file, delimiter="\n"))
    file.close()
    tissues_list=unique(tissues_list)
    tissues_list = [item for sublist in tissues_list for item in sublist]
    #tissues_list = [elem.replace("'", "") for elem in tissues_list]

    print(tissues_list)

    breakcount=0
    total_skipped=0
    breaknow=0
    #fo = open(outfile, 'w')

    keepsite=[]
    keeptissue=[]
    keeptissuemean=[]
    keeprest5perquant=[]
    keeprest95perquant=[]
    keeploworhigh=[]
    num_of_unique=1
    with open(infile) as fp: # scan methylation file
            while True:
                cur_line = fp.readline()
                if not cur_line:
                    break
                breakcount=breakcount+1
                # skip to line
                #if breakcount!=2 and breakcount<300000:
                #    if (breakcount % 250) == 0:
                #        print(cur_line[0] + ": analyzing...." + str(breakcount))
                #    continue
                #if breakcount==260000:
                #    break
                if breakcount==2: # if we're in the tissue's list line
                    tissues_line=cur_line.split("\t")[1:]
                    header_line=cur_line.split("\t")
                    header_line.insert(0,'tissue')
                    print(tissues_line)
                    #fo.write(f"{header_line}\n")
                    continue
                elif breakcount<2:
                     continue

                cur_line=cur_line.split("\t")

                if (breakcount % 250) == 0:
                    print(cur_line[0] +": analyzing...." +str(breakcount))


                #if breakcount<248630 and breakcount>248620:
                #    print("***********************************"+cur_line[0])
                #if cur_line[0]!="cg16667631":
                    #print(cur_line[0])
                #    continue
                #else:
                #    print("*********FOUND**************************" + cur_line[0])
                #    breaknow=1
                cur_line = [float('nan') if element == 'NA' else element for element in cur_line]
                cur_line = [float('nan') if element == 'NA\n' else element for element in cur_line]
                cur_line_nums=cur_line[1:] # omit cg site name



                cur_line_nums = [float(s) if isinstance(s, str) else s for s in cur_line_nums]
                #print(cur_line_nums)
                count_nans=np.sum(np.isnan(cur_line_nums))
                #print(len(cur_line_nums))
                if count_nans>3323: # skip sites with less than 2000 values
                    print("skipping "+cur_line[0]+ ", nan number:"+str(count_nans))
                    #print(cur_line)
                    total_skipped=total_skipped+1
                    #print(cur_line_nums)
                    continue

                #find tissues locs in list
                try:
                    for sel_tissue in tissues_list:
                        #cur_line_help=[]
                        cur_tissue_indices=list(indexes(tissues_line,sel_tissue))
                        cur_tissue_rest=list(rest_indexes(tissues_line,sel_tissue))
                        res_list = [cur_line_nums[i] for i in cur_tissue_indices]
                        rest_list= [cur_line_nums[i] for i in cur_tissue_rest]
                        tissue_mean=np.nanmean(res_list)  # mean for selected tissue
                        rest_mean=np.nanmean(rest_list)
                        botfive = np.nanquantile(rest_list,0.05)  # quantile(rest_of_tissues,0.05)
                        topfive = np.nanquantile(rest_list, 0.95)


                        if (tissue_mean < (botfive-difference_threshold)) or (tissue_mean > (topfive+difference_threshold)):
                            print(cur_line[0] +","+str(num_of_unique)+"/"+str(breakcount)+": unique site - " + sel_tissue + " mean/bot/top: " + str(tissue_mean) + "/" + str(
                                botfive) + "/" + str(topfive))
                            num_of_unique=num_of_unique+1
                            keepsite.append(cur_line[0])
                            keeptissue.append(sel_tissue)
                            keeptissuemean.append(tissue_mean)
                            keeprest5perquant.append(botfive)
                            keeprest95perquant.append(topfive)
                            if (tissue_mean < (botfive-difference_threshold)):
                                keeploworhigh.append("low")
                            else:
                                keeploworhigh.append("high")
                except:
                    print("failed anayling for site "+cur_line[0]+" which is on line "+str(breakcount))





    print(" total methylaytion sites skipped:" +str(total_skipped))
    #fo.close()
    #output_dict=dict()
    output_dict={'site':keepsite,'unique tissue':keeptissue,'tissue mean':keeptissuemean,'0.05 quantile, rest of tissues':keeprest5perquant
        , '0.95 quantile, rest of tissues': keeprest95perquant,'unique low/high':keeploworhigh}
    print(output_dict)
    data = pd.DataFrame(output_dict)
    # export data to csv
    data.to_csv(outfile)

