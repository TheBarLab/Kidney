import pandas as pd
import matplotlib.pyplot as plt
import numpy.ma as ma
import sys

meth_data_file='/home/naorsagy/Desktop/20230503_corr_analysis_for_first_paper/GSE50874_series_matrix_with_IF_eGFR.csv'
#meth_data=pd.read_csv('/home/naorsagy/Desktop/20230503_corr_analysis_for_first_paper/GSE50874_series_matrix_with_IF_eGFR.csv')
kidney_unique_sites=pd.read_csv('/home/naorsagy/Desktop/20230503_corr_analysis_for_first_paper/kidney_unique_sites_0_2.csv')
corrs_file=pd.read_csv('/home/naorsagy/Desktop/20230503_corr_analysis_for_first_paper/20230510_GSE50874_corrs.csv')



unique_sites_if_corr=[]
unique_sites_egfr_corr=[]
non_unique_sites_if_corr=[]
non_unique_sites_egfr_corr=[]
unique_counter=0
non_unique_counter=0
for index, row in corrs_file.iterrows():

    #if row['site'] in kidney_unique_sites['site']: # if unique
    if (kidney_unique_sites['site'].eq(row['site'])).any():
        unique_sites_if_corr.append(row["IF corr"])
        unique_sites_egfr_corr.append(row["eGFR corr"])
        unique_counter+=1
        print(str(index)+": ("+str(unique_counter)+"/"+str(non_unique_counter)+") -> "+row['site'] + ' is unique ')
    else:
        non_unique_sites_if_corr.append(row["IF corr"])
        non_unique_sites_egfr_corr.append(row["eGFR corr"])
        non_unique_counter+=1
        print(str(index)+": ("+str(unique_counter)+"/"+str(non_unique_counter)+") -> "+row['site'] + ' is not unique ')

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(1, 1, 1)
ax.set_xlabel('eGFR r', fontsize=15)
ax.set_ylabel('IF r', fontsize=15)
ax.set_title('IF r vs. eGFR r', fontsize=20)
colors={'low':'royalblue', 'high':'maroon'}


ax.scatter(non_unique_sites_egfr_corr
                               , non_unique_sites_if_corr
                               , c='royalblue'
                               , s=20, marker='o',alpha=0.15,linewidths=0.0)
ax.scatter(unique_sites_egfr_corr
                               , unique_sites_if_corr
                               , c='maroon'
                               , s=20, marker='o',alpha=0.6,linewidths=0.0)
plt.show()


sys.exit()

sites=[]
eGFR_corr_save=[]
IF_corr_save=[]
with open(meth_data_file) as f:
    subjectid = f.readline()
    subjectid=subjectid.split(',')
    print(subjectid)

    eGFR = f.readline()
    eGFR=eGFR.split(',')
    eGFR=eGFR[1:]
    eGFR = [i.replace("\n", "") for i in eGFR]
    eGFR = [float("NAN") if i=="NA" else float(i) for i in eGFR]
    print(eGFR)

    IF = f.readline()
    IF = IF.split(',')
    IF = IF[1:]
    IF = [i.replace("\n","") for i in IF]
    IF = [float("NAN") if i=="NA" else float(i) for i in IF]
    print(IF)

    count=0
    while True:
        count += 1
        line = f.readline()

        if not line:  # or count>5:
            break

        line = line.split(',')
        site=line[0]
        line=line[1:]
        line = [i.replace("\n", "") for i in line]
        #print(line)
        line = [float("NAN") if (i == "NA" or i =="") else float(i) for i in line]
        a = ma.masked_invalid(line)
        b = ma.masked_invalid(IF)
        msk = (~a.mask & ~b.mask)
        c = ma.masked_invalid(eGFR)
        msk2 = (~a.mask & ~c.mask)
        if_corr=str(ma.corrcoef(a[msk], b[msk])[1, 0])
        egfr_corr=str(ma.corrcoef(a[msk2], c[msk2])[1, 0])
        print(str(count)+","+site+" -> IF corr: "+str(if_corr) + " ; GFR corr:"+str(egfr_corr))
        sites.append(site)
        eGFR_corr_save.append(egfr_corr)
        IF_corr_save.append(if_corr)


dict={'site':sites,'IF corr':IF_corr_save,'eGFR corr':eGFR_corr_save}
df=pd.DataFrame(dict)
df.to_csv('/home/naorsagy/Desktop/20230503_corr_analysis_for_first_paper/20230510_GSE50874_corrs.csv')

