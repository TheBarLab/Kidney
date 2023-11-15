import pandas as pd
import sys
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import seaborn as sns

if_corr_data=pd.read_csv('/home/naorsagy/Desktop/20230503_corr_analysis_for_first_paper/Interstitium_Fibrosis_Correlations.csv')
meta_data=pd.read_csv('/home/naorsagy/Desktop/20230503_corr_analysis_for_first_paper/chinese_v1_metadata.csv')
kidney_unique_sites=pd.read_csv('/home/naorsagy/Desktop/20230503_corr_analysis_for_first_paper/kidney_unique_sites_0_2.csv')
kidney_meth_vals_unique=pd.read_csv('/home/naorsagy/Desktop/20230503_corr_analysis_for_first_paper/20220914_unique_cg_kidney.csv')
#print(meta_data)


# first line: subject id
# second line: tissue
v1_file='/home/naorsagy/Desktop/R_Methyl/tissue_methylation_v1.txt'

with open(v1_file) as f:
    subjectid = f.readline()
    subjectid=subjectid.split('\t')
    print(subjectid)

    tissue = f.readline()
    tissue=tissue.split('\t')
    print(tissue)


# find where the kidney tissues are
indices = [i for i, x in enumerate(tissue) if x == "kidney"]   # where all the kidneys are
subjectids=[subjectid[i] for i in indices]                   # check what we got

print(meta_data)

# save ages
save_ages=[]
for indx,subject in enumerate(subjectids):
    curindx=meta_data[meta_data['sample_id']== subject].index[0]
    save_ages.append(meta_data.iloc[curindx, meta_data.columns.get_loc("age")])
    print(str(indx) + ":" + str(curindx) + " : "+str(subject)+" -> "+str(meta_data.iloc[curindx, meta_data.columns.get_loc("age")]))

# scan unique sites
save_corrs=[]
for index, row in kidney_meth_vals_unique.iterrows():
    meth_vals_for_site = [row[i] for i in indices]
    res = [float(x) for x in meth_vals_for_site]
    a = ma.masked_invalid(res)
    b = ma.masked_invalid(save_ages)
    msk = (~a.mask & ~b.mask)
    print(ma.corrcoef(a[msk], b[msk])[1,0])
    save_corrs.append(ma.corrcoef(a[msk], b[msk])[1,0])

print(res)
print(save_ages)


print(if_corr_data)
print(kidney_unique_sites)

save_if_corr_for_unique_sites=[]
# screen unique sites IF
for indx,site in enumerate(kidney_unique_sites['site']):
    curindx=if_corr_data[if_corr_data['site'] == site].index[0]
    print(str(indx)+":"+str(curindx)+" : "+ str(if_corr_data.iloc[curindx, if_corr_data.columns.get_loc("IF correlation full")]))
    save_if_corr_for_unique_sites.append(if_corr_data.iloc[curindx, if_corr_data.columns.get_loc("IF correlation full")])

print(save_if_corr_for_unique_sites)
# dictionary of lists
dict = {'if': save_if_corr_for_unique_sites}
df = pd.DataFrame(dict)
# saving the dataframe
df.to_csv('GFG.csv')

data_to_plot=[save_if_corr_for_unique_sites,save_corrs]
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(1, 1, 1)
ax.set_xlabel('', fontsize=15)
ax.set_ylabel('Correlation', fontsize=15)
ax.set_title('Correlation distribution for IF and Age', fontsize=20)
colors={'low':'royalblue', 'high':'maroon'}
bp = ax.violinplot(data_to_plot,showmeans=False, showmedians=False,
        showextrema=False,widths=0.5,points=200,bw_method='silverman')
xtick_loc = [1,2]
ax.set_xticks(xtick_loc)
labels=["IF","Age"]
ax.set_xticklabels(labels)
#plt.ylim(-1,1)
plt.show()


ax=sns.violinplot(data_to_plot,alpha=0.5,inner=None)
for violin, alpha in zip(ax.collections[::2], [0.6,0.6]):
    violin.set_alpha(0.3)

#ax.set_xlabel('', fontsize=15)
ax.set_ylabel('Correlation', fontsize=20)
xtick_loc = [0,1]
ax.set_xticks(xtick_loc)
labels=["IF","Age"]
ax.set_xticklabels(labels)
plt.xticks(fontsize=15)
plt.show()