import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pickle
import pandas as pd
import numpy as np

#import sys
#sys.path.insert(0, './utility/')

def save_obj( data, name ):
    with open(name,'wb') as f:
        pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL);


def load_obj( name ):
    with open(name,'rb') as f:
        return pickle.load(f);



#ceres 23q4 dependency
'''f = open("data/CRISPRGeneDependency.csv", "r")

#extract gene name and remove HGNC
gene = f.readline()
gene = np.array(gene.replace('\n','').replace('.','_').split(','));
gene = gene[1:];
entrez = [];
i = 0;
while i <len(gene):
    start = 0;
    start = gene[i].find(' ',0);
    entrez.append(gene[i][start+2:len(gene[i])-1])
    gene[i] = gene[i][0:start];
    i = i + 1;
entrez = pd.DataFrame(entrez,index = gene, columns = ['entrez'])
entrez.to_pickle('data/entrez.pkl')
print(entrez)

#extract cell line and values
cell = f.read();
f.close();
cell = cell.split('\n');
cell = cell[0:len(cell)-1];
values = [];
cell_name = [];
for i in range(0,len(cell)):
    temp = cell[i].split(',');
    cell_name.append(temp[0]);
    values.append(temp[1:]);

#make gene and cell into pickle
dep = pd.DataFrame( values, index = cell_name, columns = gene);
dep = dep.replace({'NA':np.nan})
dep = dep.replace({'':np.nan})
dep = dep.astype('float')
print(dep)
dep.to_pickle('data/ceres_dep.pkl')
#'''

#extracting ceres gene effect
'''f = open("data/CRISPRGeneEffect.csv", "r")

#extract gene name and remove HGNC
gene = f.readline()
gene = np.array(gene.replace('\n','').replace('.','_').split(','));
gene = gene[1:];
entrez = [];
i = 0;
while i <len(gene):
    start = 0;
    start = gene[i].find(' ',0);
    entrez.append(gene[i][start+2:len(gene[i])-1])
    gene[i] = gene[i][0:start];
    i = i + 1;
entrez = pd.DataFrame(entrez,index = gene, columns = ['entrez'])
#entrez.to_pickle('data/entrez.pkl')
print(entrez)

#extract cell line and values
cell = f.read();
f.close();
cell = cell.split('\n');
cell = cell[0:len(cell)-1];
values = [];
cell_name = [];
for i in range(0,len(cell)):
    temp = cell[i].split(',');
    cell_name.append(temp[0]);
    values.append(temp[1:]);

#make gene and cell into pickle
dep = pd.DataFrame( values, index = cell_name, columns = gene);
dep = dep.replace({'NA':np.nan})
dep = dep.replace({'':np.nan})
dep = dep.astype('float')
print(dep)
dep.to_pickle('data/ceres_effect.pkl')
#'''

'''dep = pd.read_pickle('data/ceres_effect.pkl')
print(dep.shape)
#'''


#cnv
#f = open("data/cnv_22Q4", "r")
'''f = open("data/23Q4/OmicsCNGene.csv", "r")
#extract gene name and remove HGNC
gene = f.readline()
gene = np.array(gene.replace('\n','').replace('.','_').split(','));
gene = gene[1:];
i = 0;
while i <len(gene):
    start = 0;
    start = gene[i].find(' ',0);
    gene[i] = gene[i][0:start];
    i = i + 1;
print(gene)

#extract cell line and values
cell = f.read();
f.close();
cell = cell.split('\n');
cell = cell[0:len(cell)-1];
values = [];
cell_name = [];
for i in range(0,len(cell)):
    temp = cell[i].split(',');
    cell_name.append(temp[0]);
    values.append(temp[1:]);

#make gene and cell into pickle
exp = pd.DataFrame( values, index = cell_name, columns = gene);
exp = exp.replace({'NA':np.nan})
exp = exp.replace({'':np.nan})
exp = exp.astype('float')
print(exp)
exp.to_pickle('data/cnv_original.pkl')
#'''

#mut
'''f = open("data_2023/OmicsSomaticMutations.csv", "r")
#extract gene name and remove HGNC
gene = f.readline()
gene = np.array(gene.replace('\n','').split(','));
print(len(gene))
import re
from csv import reader
#extract cell line and values
cell = f.read();
f.close();
cell = cell.split('\n');
cell = cell[0:len(cell)-1];
values = [];
for i in range(0,len(cell)):
    #temp = cell[i].split(',');
    #temp = list(filter(None, re.split(r',|"(.*?)"', cell[i])))
    temp = reader([cell[i]])
    temp = next(temp)
    #temp = re.split(r',(?=")', cell[i])
    values.append(temp);
    if len(temp) != 64:
        print(len(temp))
        print(temp)
        #print(cell[i])
#make gene and cell into pickle
exp = pd.DataFrame( values, columns = gene);
exp = exp.replace({'':np.nan})
print(exp)
exp.to_pickle('data_2023/mut_original.pkl')
#'''

#Extracting sample information
'''sample_info = pd.read_csv('data/23Q4/Model.csv', index_col = 0)
print(sample_info)
sample_info.to_pickle('data/sample.pkl')
#'''

#Extracting Breast cancer info
'''sample = pd.read_pickle('data/sample.pkl')
dep = pd.read_pickle('data/ceres_dep.pkl')
sample = sample.reindex( index = dep.index)
sample = sample.loc[ sample.loc[:,'OncotreeLineage'] == 'Breast',:]
sample = sample.loc[:, ['StrippedCellLineName','LegacyMolecularSubtype','LegacySubSubtype']]
sample.loc['ACH-001065','LegacySubSubtype'] = 'HER2_pos'
sample.loc['ACH-001419','LegacySubSubtype'] = 'HER2_pos'
sample.loc['ACH-002399','LegacySubSubtype'] = 'HER2_pos'
sample.loc['ACH-001820','LegacySubSubtype'] = 'TNBC'
sample = sample.loc[ ~sample.loc[:,'LegacySubSubtype'].isna(),:]
sample.loc[ sample.loc[:,'LegacySubSubtype'].str.contains('HER2pos'),'LegacySubSubtype'] = 'HER2_pos'
sample.loc[ sample.loc[:,'LegacySubSubtype'] == 'ERneg_HER2neg','LegacySubSubtype'] = 'TNBC'
sample.loc[ sample.loc[:,'LegacySubSubtype'] == 'ERpos_HER2neg','LegacySubSubtype'] = 'ER_pos'
print(sample)
#'''

#Extracting dgidb database
'''dgidb = pd.read_csv('data/23Q4/interactions.tsv', sep = '\t')
print(dgidb)
dgidb.to_pickle('data/dgidb_interactions.pkl')
#'''
#query whether the gene is druggable genome
'''g = pd.read_csv('data/23Q4/genes_DGIdb.tsv', sep = '\t')
print(len(g.loc[:,'gene_claim_name'].unique()))
g = g.loc[:,'gene_claim_name'].unique()
druggable = pd.DataFrame( np.nan, index = g, columns = ['drug'])

import requests
headers = {
    # Already added when you pass json= but not when you pass data=
    # 'Content-Type': 'application/json',
}
print(druggable)
for i in druggable.index:
    tmp_s = '{ genes(names: [\"'+str(i) +'\"]) {nodes {name geneCategoriesWithSources {name} } } }'
    json_data = {
        'query': tmp_s,
    }
    response = requests.post('https://dgidb.org/api/graphql', headers=headers, json=json_data)
    if 'DRUGGABLE GENOME' in response.text:
        druggable.loc[i,'drug'] = 1
    else:
        druggable.loc[i,'drug'] = 0
# Note: json_data will not be serialized by requests
# exactly as it was in the original request.
#data = '{ "query": "{ genes(first:20000) {nodes {name geneCategoriesWithSources {name} } } }" }'
#response = requests.post('https://dgidb.org/api/graphql', headers=headers, data=data)
druggable.to_pickle('data/druggable_genome.pkl')
print(druggable)
print(druggable.sum())
#'''
#extending to all categories
'''g = pd.read_csv('data/23Q4/genes_DGIdb.tsv', sep = '\t')
print(len(g.loc[:,'gene_claim_name'].unique()))
g = g.loc[:,'gene_claim_name'].unique()
druggable = pd.DataFrame( np.nan, index = g, columns = ['drug'])

import requests
headers = {
    # Already added when you pass json= but not when you pass data=
    # 'Content-Type': 'application/json',
}
print(druggable)
j = 0
for i in druggable.index:
    print(j)
    tmp_s = '{ genes(names: [\"'+str(i) +'\"]) {nodes {name geneCategoriesWithSources {name} } } }'
    json_data = {
        'query': tmp_s,
    }
    j = j + 1
    response = requests.post('https://dgidb.org/api/graphql', headers=headers, json=json_data)
    if 'Categories' in response.text:
        druggable.loc[i,'drug'] = 1
    else:
        druggable.loc[i,'drug'] = 0
# Note: json_data will not be serialized by requests
# exactly as it was in the original request.
#data = '{ "query": "{ genes(first:20000) {nodes {name geneCategoriesWithSources {name} } } }" }'
#response = requests.post('https://dgidb.org/api/graphql', headers=headers, data=data)
druggable.to_pickle('data/druggable_genome_correct.pkl')
print(druggable)
print(druggable.sum())
#'''
#fix naming problem in the system
'''druggable = pd.read_pickle('data/druggable_genome_correct.pkl')
druggable = druggable.loc[ ~druggable.index.isna(),:]
drug_old = druggable.loc[ ~druggable.index.isin(reg),:] 
newg = []
for i in reg:
    tmp = i
    tmp = tmp.replace('.00','')
    tmp = tmp.replace('0','')
    tmp = tmp[0:3]+tmp[4]
    newg.append(tmp)
print(newg)
druggable = pd.DataFrame( np.nan, index = newg, columns = ['drug'])

import requests
headers = {
    # Already added when you pass json= but not when you pass data=
    # 'Content-Type': 'application/json',
}
print(druggable)
for i in druggable.index:
    tmp_s = '{ genes(names: [\"'+str(i) +'\"]) {nodes {name geneCategoriesWithSources {name} } } }'
    json_data = {
        'query': tmp_s,
    }
    response = requests.post('https://dgidb.org/api/graphql', headers=headers, json=json_data)
    #if 'DRUGGABLE GENOME' in response.text:
    if 'Categories' in response.text:
        druggable.loc[i,'drug'] = 1
    else:
        druggable.loc[i,'drug'] = 0
druggable = pd.concat([ druggable, drug_old], axis = 0)
print(druggable)
print(druggable.sum())
druggable.to_pickle('data/druggable_genome_update_correct.pkl')
#'''

#TCGA expression data
'''tcga = pd.read_csv('data/TCGA_BRCA_clinical.txt', sep = '\t', index_col = 0)
print(tcga)
print(tcga.loc[:,['ER_status','HER2_status','TNBC']])
tcga.to_pickle('data/TCGA_sample.pkl')
exp = pd.read_csv('data/TCGA_BRCA_expression_status.txt', sep = '\t', index_col = 0)
exp.to_pickle('data/exp_TCGA.pkl')
#'''

'''exp = pd.read_pickle('data/exp_TCGA.pkl')
tcga = pd.read_pickle('data/TCGA_sample.pkl')
print(exp)
s = exp.columns.str.split('-')
tmp = []
for i in s:
    tmp.append( i[0]+'-'+i[1]+'-'+i[2])
expori = exp.copy()
exp.columns = tmp
tcga = tcga.loc[:,['ER_status','HER2_status','TNBC']].dropna(how = 'any')
HER2 = tcga.loc[tcga.loc[:,'HER2_status'] == 'Positive',:]
TNBC = tcga.loc[tcga.loc[:,'TNBC'] == 'Yes',:]
ER = tcga.loc[~((tcga.index.isin(HER2.index)) | (tcga.index.isin(TNBC.index))),:]
print(HER2)
print(TNBC)
print(ER)
ER = expori.loc[:,exp.columns.isin(ER.index)]
HER2 = expori.loc[:,exp.columns.isin(HER2.index)]
TNBC = expori.loc[:,exp.columns.isin(TNBC.index)]
print(HER2)
print(TNBC)
print(ER)
tmp = [ER, HER2, TNBC]
#save_obj(tmp, 'tcga_exp_type') 
print(tmp)
print(tcga)
exp = exp.loc[:, exp.columns.isin(tcga.index)]
print(exp)
print(exp.columns[exp.columns.duplicated()])
import matplotlib.backends.backend_pdf
pdf = matplotlib.backends.backend_pdf.PdfPages('figure/exp_histogram.pdf')
for i in range(3):
    fig = plt.figure()
    plt.hist(tmp[i].sum(axis=1)/len(tmp[i].columns), bins = 20)
    pdf.savefig(fig)
    plt.close()
pdf.close()
#'''

#tcga subtype
'''tmp = load_obj('data/tcga_exp_type') 
print(tmp)
import matplotlib.backends.backend_pdf
pdf = matplotlib.backends.backend_pdf.PdfPages('figure/exp_histogram.pdf')
name = ['ER+','HER2+','TNBC']
for i in range(3):
    fig = plt.figure()
    plt.hist(tmp[i].sum(axis=1)/len(tmp[i].columns), bins = 20)
    plt.xlabel('expression percentage of TCGA pts')
    plt.ylabel('number of genes')
    s = name[i]+ ' n=' + str(len(tmp[i].columns))
    plt.title(s)
    pdf.savefig(fig)
    plt.close()
pdf.close()
#'''

#subtype assignment table
'''sample = pd.read_pickle('data/sample.pkl')
dep = pd.read_pickle('data/ceres_dep.pkl')
sample = sample.reindex( index = dep.index)
sample = sample.loc[ sample.loc[:,'OncotreeLineage'] == 'Breast',:]
sample = sample.loc[:, ['StrippedCellLineName','LegacySubSubtype']]
sample.columns = ['cell_line_name','DepMap_subtype']
sample.loc['ACH-001065','Clinical_Subtype'] = 'HER2_pos'
sample.loc['ACH-001419','Clinical_Subtype'] = 'HER2_pos'
sample.loc['ACH-002179','Clinical_Subtype'] = 'HER2_pos'
sample.loc['ACH-002399','Clinical_Subtype'] = 'HER2_pos'
sample.loc['ACH-001820','Clinical_Subtype'] = 'TNBC'
print(sample)
sample.loc[:,'DepMap_subtype'].fillna('', inplace = True)
sample.loc[ sample.loc[:,'DepMap_subtype'].str.contains('HER2pos'),'Clinical_Subtype'] = 'HER2_pos'
sample.loc[ sample.loc[:,'DepMap_subtype'] == 'ERneg_HER2neg','Clinical_Subtype'] = 'TNBC'
sample.loc[ sample.loc[:,'DepMap_subtype'] == 'ERpos_HER2neg','Clinical_Subtype'] = 'ER_pos'
print(sample)
sample.to_csv('data/cell_line_subtype.csv')
#'''

###############################
#extracting prioritizing genes#
###############################
#extracting priotorizing genes according to subtype, using prob thres at 0.5
'''sample = pd.read_pickle('data/sample.pkl')
dep = pd.read_pickle('data/ceres_dep.pkl')
sample = sample.reindex( index = dep.index)
sample = sample.loc[ sample.loc[:,'OncotreeLineage'] == 'Breast',:]
sample = sample.loc[:, ['StrippedCellLineName','LegacyMolecularSubtype','LegacySubSubtype']]
sample.loc['ACH-001065','LegacySubSubtype'] = 'HER2_pos'
sample.loc['ACH-001419','LegacySubSubtype'] = 'HER2_pos'
sample.loc['ACH-002179','LegacySubSubtype'] = 'HER2_pos'
sample.loc['ACH-002399','LegacySubSubtype'] = 'HER2_pos'
sample.loc['ACH-001820','LegacySubSubtype'] = 'TNBC'
sample = sample.loc[ ~sample.loc[:,'LegacySubSubtype'].isna(),:]
sample.loc[ sample.loc[:,'LegacySubSubtype'].str.contains('HER2pos'),'LegacySubSubtype'] = 'HER2_pos'
sample.loc[ sample.loc[:,'LegacySubSubtype'] == 'ERneg_HER2neg','LegacySubSubtype'] = 'TNBC'
sample.loc[ sample.loc[:,'LegacySubSubtype'] == 'ERpos_HER2neg','LegacySubSubtype'] = 'ER_pos'
#print(sample)
print(sample.loc[:,'LegacySubSubtype'].value_counts())
eff = pd.read_pickle('data/ceres_dep.pkl')
eff = eff.reindex(index = sample.index)
eff = eff.loc[:,eff.isna().sum() < len(eff.index) * 0.2]
dep = pd.read_pickle('data/ceres_effect.pkl')
dep = dep.reindex(index = sample.index)
dep = dep.reindex(columns = eff.columns)
er = eff.loc[ sample.loc[:,'LegacySubSubtype'] == 'ER_pos',:]
her = eff.loc[ sample.loc[:,'LegacySubSubtype'] == 'HER2_pos',:]
tnbc = eff.loc[ sample.loc[:,'LegacySubSubtype'] == 'TNBC',:]
er = er.loc[:,(er >= 0.5).sum() >= 1]
her = her.loc[:,(her >= 0.5).sum() >= 1]
tnbc = tnbc.loc[:,(tnbc >= 0.5).sum() >= 1]
druggable = pd.read_pickle('data/druggable_genome_update_correct.pkl')
druggable = druggable.loc[druggable.loc[:,'drug'] == 1,:]
er = er.loc[:,er.columns.isin(druggable.index)]
her = her.loc[:,her.columns.isin(druggable.index)]
tnbc = tnbc.loc[:,tnbc.columns.isin(druggable.index)]

tcga = load_obj('data/tcga_exp_type')
thres = 0.95
er = er.loc[:,er.columns.isin(tcga[0].index[ tcga[0].sum(axis=1) >= len(tcga[0].columns)*thres ])]
her = her.loc[:,her.columns.isin(tcga[1].index[ tcga[1].sum(axis=1) >= len(tcga[1].columns)*thres ])]
tnbc = tnbc.loc[:,tnbc.columns.isin(tcga[2].index[ tcga[2].sum(axis=1) >= len(tcga[2].columns)*thres ])]


er_gene = pd.DataFrame( np.nan, index = er.columns, columns = ['perc','median_prob','assignment'])
her_gene = pd.DataFrame( np.nan, index = her.columns, columns = ['perc','median_prob','assignment'])
tnbc_gene = pd.DataFrame( np.nan, index = tnbc.columns, columns = ['perc','median_prob','assignment'])
ceres_control = pd.read_csv('data/23Q4/common_ess_2021.csv', sep = ',',  index_col = 0)
achilles_control = pd.read_csv('data/23Q4/AchillesCommonEssentialControls.csv', sep = ' ', index_col = 0)
er_gene.loc[:,'perc'] = (er >= 0.5).sum()/len(er.index)
her_gene.loc[:,'perc'] = (her >= 0.5).sum()/len(her.index)
tnbc_gene.loc[:,'perc'] = (tnbc >= 0.5).sum()/len(tnbc.index)
#er_gene.loc[:,'median_eff'] =  dep.reindex(index = er.index, columns = er_gene.index).median()
#her_gene.loc[:,'median_eff'] =  dep.reindex(index = her.index, columns = her_gene.index).median()
#tnbc_gene.loc[:,'median_eff'] =  dep.reindex(index = tnbc.index, columns = tnbc_gene.index).median()
er_gene.loc[:,'median_prob'] =  eff.reindex(index = er.index, columns = er_gene.index).median()
her_gene.loc[:,'median_prob'] =  eff.reindex(index = her.index, columns = her_gene.index).median()
tnbc_gene.loc[:,'median_prob'] =  eff.reindex(index = tnbc.index, columns = tnbc_gene.index).median()
#er_gene.loc[:,'median_prob'] =  eff.reindex(index = er.index, columns = er_gene.index).mean()
#her_gene.loc[:,'median_prob'] =  eff.reindex(index = her.index, columns = her_gene.index).mean()
#tnbc_gene.loc[:,'median_prob'] =  eff.reindex(index = tnbc.index, columns = tnbc_gene.index).mean()
er_gene.loc[ er_gene.loc[:,'perc']*8 >= 2, 'assignment'] = 'priority'
her_gene.loc[ her_gene.loc[:,'perc'] >= 0.1, 'assignment'] = 'priority'
tnbc_gene.loc[ tnbc_gene.loc[:,'perc'] >= 0.1, 'assignment'] = 'priority'
#er_gene.loc[ er_gene.loc[:,'perc']*8 < 2, 'assignment'] = 'lower 10%'
#her_gene.loc[ her_gene.loc[:,'perc'] < 0.1, 'assignment'] = 'lower 10%'
#tnbc_gene.loc[ tnbc_gene.loc[:,'perc'] < 0.1, 'assignment'] = 'lower 10%'
er_gene.loc[ er_gene.loc[:,'perc']*8 < 2, 'assignment'] = 'not priority'
her_gene.loc[ her_gene.loc[:,'perc'] < 0.1, 'assignment'] = 'not priority'
tnbc_gene.loc[ tnbc_gene.loc[:,'perc'] < 0.1, 'assignment'] = 'not priority'
er_gene.loc[ er_gene.loc[:,'median_prob'] < 0.25, 'assignment'] = 'not priority'
her_gene.loc[ her_gene.loc[:,'median_prob'] < 0.25, 'assignment'] = 'not priority'
tnbc_gene.loc[ tnbc_gene.loc[:,'median_prob'] < 0.25, 'assignment'] = 'not priority'
er_gene.loc[ er_gene.index.isin(achilles_control.index) | er_gene.index.isin(ceres_control.index)  , 'assignment'] = 'CommonEssCore'
her_gene.loc[ her_gene.index.isin(achilles_control.index) | her_gene.index.isin(ceres_control.index)  , 'assignment'] = 'CommonEssCore'
tnbc_gene.loc[ tnbc_gene.index.isin(achilles_control.index) | tnbc_gene.index.isin(ceres_control.index)  , 'assignment'] = 'CommonEssCore'
ceres_control = pd.read_csv('data/23Q4/CRISPRInferredCommonEssentials.csv', sep = ' ',  index_col = 0)
er_gene.loc[  er_gene.index.isin(ceres_control.index)  , 'assignment'] = 'CommonEssCore'
her_gene.loc[  her_gene.index.isin(ceres_control.index)  , 'assignment'] = 'CommonEssCore'
tnbc_gene.loc[  tnbc_gene.index.isin(ceres_control.index)  , 'assignment'] = 'CommonEssCore'
print(er_gene)
print(her_gene)
print(tnbc_gene)
er_gene.to_csv('R_data/er_gene_lowthres.csv')
her_gene.to_csv('R_data/her_gene_lowthres.csv')
tnbc_gene.to_csv('R_data/tnbc_gene_lowthres.csv')
er_gene.loc[:,'gene'] = er_gene.index
her_gene.loc[:,'gene'] = her_gene.index
tnbc_gene.loc[:,'gene'] = tnbc_gene.index
er_gene.loc[:,'molecular'] = 'ER+'
her_gene.loc[:,'molecular'] = 'HER2+'
tnbc_gene.loc[:,'molecular'] = 'TNBC'
print(er_gene.loc[:,'assignment'].value_counts())
print(her_gene.loc[:,'assignment'].value_counts())
print(tnbc_gene.loc[:,'assignment'].value_counts())
er_gene = pd.concat([er_gene, her_gene, tnbc_gene], axis = 0, ignore_index =  True)
er_gene.to_csv('R_data/all_gene_prob_lowthres.csv')
#'''

#generating priority genes
'''
gene = pd.read_csv('R_data/all_gene_prob_lowthres.csv', index_col=0 )
gene = gene.loc[ gene.loc[:,'assignment'] == 'priority',:]
print(gene)
gene.loc[:,'anno'] = ''
gene.loc[:,'sel'] = 'above'
gene.loc[ gene.loc[:,'median_prob'] < 0.5, 'sel'] = 'below'
print(gene.loc[gene.loc[:,'sel'] == 'above',:])
for i in ['ER+','HER2+','TNBC']:
    tmp = gene.loc[ (gene.loc[:,'molecular'] == i) & (gene.loc[:,'sel'] == 'above'),:]
    gene.loc[tmp.loc[:,'median_prob'].nlargest(10).index,'anno'] = gene.loc[tmp.loc[:,'median_prob'].nlargest(10).index,'gene']
print(gene)
gene.to_csv('R_data/all_gene_prob_thres_lowthres.csv')
#'''

'''
gene = pd.read_csv('R_data/all_gene_prob_thres.csv', index_col = 0)
print(gene)
for i in ['ER+','HER2+','TNBC']:
    tmp = gene.loc[ (gene.loc[:,'molecular'] == i),:]
    print(tmp.loc[:,'sel'].value_counts())
gene.rename( columns = {'sel':'Top priority'}
#'''

#unclear if will use
#supplementary tables subtype
'''
gene = pd.read_csv('R_data/all_gene_prob_thres_drug_lowthres.csv', index_col = 0)
print(gene)
gene = gene.rename( columns = {'perc':'fraction of cell lines essential','median_prob':'median dependency score','sel':'top priority'})
output = []
for i in ['ER+','HER2+','TNBC']:
    tmp = gene.loc[gene.loc[:,'molecular'] == i,:]
    tmp.index = tmp.loc[:,'gene']
    print(tmp.loc[:,'top priority'].value_counts())
    tmp.loc[tmp.loc[:,'top priority'] == 'below','top priority'] = ''
    tmp.loc[tmp.loc[:,'top priority'] == 'above','top priority'] = 'yes'
    #tmp = tmp.drop(columns = ['gene','anno','assignment','Oncogene','TSG'])
    tmp = tmp.drop(columns = ['gene','anno','assignment'])
    print(tmp)
    output.append(tmp)
    #tmp.to_csv('R_data/table_'+i+'.csv')
output = pd.concat(output, axis = 0)
output = output.sort_values(by = ['molecular','median dependency score'], ascending = [True, False])
print(output)
output.to_csv('R_data/table_all_combined.csv')
#'''

##########################
#analyzing priority genes#
##########################
#MSK onco info
'''gene = pd.read_csv('R_data/all_gene_prob_thres.csv', index_col = 0)
onco = pd.read_csv('data/23Q4/cancerGeneList.tsv', sep = '\t')
onco_gene = onco.loc[onco.loc[:,'Is Oncogene'] == 'Yes',:]
gene.loc[:,'Oncogene'] = 'No'
gene.loc[gene.loc[:,'gene'].isin(onco_gene.loc[:,'Hugo Symbol']),'Oncogene'] = 'Yes'
onco_gene = onco.loc[onco.loc[:,'Is Tumor Suppressor Gene'] == 'Yes',:]
gene.loc[:,'TSG'] = 'No'
gene.loc[gene.loc[:,'gene'].isin(onco_gene.loc[:,'Hugo Symbol']),'TSG'] = 'Yes'
print(gene)
for i in ['ER+','HER2+','TNBC']:
    print(i)
    print( (gene.loc[ gene.loc[:,'molecular'] == i ,'Oncogene'] == 'Yes').sum())
    print( (gene.loc[ gene.loc[:,'molecular'] == i ,'TSG'] == 'Yes').sum())
gene.to_pickle('data/breast_priority_gene_oncoinfo.pkl')
#'''

#open target query
#gene = pd.read_pickle('data/breast_priority_gene_oncoinfo.pkl')
'''gene = pd.read_csv('R_data/all_gene_prob_thres_lowthres.csv')
ingene = gene.loc[:,'gene'].unique()
ingene = list(ingene)
from gprofiler import GProfiler
gp = GProfiler( return_dataframe = True)
result = gp.convert( organism = 'hsapiens', query = ingene, target_namespace = 'ENSG')
print(result)
#gene.loc[:,'open_target'] = np.nan 
gene.loc[:,'clinical_phase'] = np.nan 
gene.loc[:,'drug_name'] = '' 

import requests
import json
breast_neoplasm = 'EFO_0003869'
breast_carcinoma = 'EFO_0000305'
breast_disease = 'EFO_0009483'
breast_cancer = 'MONDO_0007254'
cancer = 'MONDO_0004992'
gene_id = "ENSG00000169083"
gene_id = "ENSG00000141510"
disease = "MONDO_0007254"
query_string_cnt =  """ query ChemblQuery($ensemblId: String!, $efoId: String!, $size: Int!) {
      disease(efoId: $efoId) {
        id
        chembl: evidences(
          ensemblIds: [$ensemblId]
          enableIndirect: true
          datasourceIds: ["chembl"]
          size: $size
        ) { count } }  }
    """
query_string_all =  """ query ChemblQuery($ensemblId: String!, $efoId: String!, $size: Int!) {
      disease(efoId: $efoId) {
        id
        chembl: evidences(
          ensemblIds: [$ensemblId]
          enableIndirect: true
          datasourceIds: ["chembl"]
          size: $size
        ) { count rows { clinicalPhase clinicalStatus } } }  }
    """
query_string_1 = """
        query ChemblQuery($ensemblId: String!, $efoId: String!, $size: Int!) {
      disease(efoId: $efoId) {
        id
        chembl: evidences(
          ensemblIds: [$ensemblId]
          enableIndirect: true
          datasourceIds: ["chembl"]
          size: $size
        ) {
          count
          rows { disease { id name } target { id approvedSymbol} drug { id name drugType mechanismsOfAction { rows { mechanismOfAction targets { id approvedSymbol}}}} clinicalPhase clinicalStatus urls{url}}}}}
    """
base_url = "https://api.platform.opentargets.org/api/v4/graphql"
headers = {
    # Already added when you pass json= but not when you pass data=
    # 'Content-Type': 'application/json',
}
for i in result.index:
    print(i)
    gene_id = result.loc[i,'converted']
    gene_name = result.loc[i,'incoming']
    variables = {"ensemblId": gene_id, "efoId": disease, "size": 1}
    r = requests.post(base_url, headers=headers, json={"query": query_string_cnt, "variables": variables})
    api_response = json.loads(r.text)
    send_cnt = api_response['data']['disease']['chembl']['count']
    if send_cnt > 0:
        variables = {"ensemblId": gene_id, "efoId": disease, "size": send_cnt}
        #r = requests.post(base_url, headers=headers, json={"query": query_string_all, "variables": variables})
        r = requests.post(base_url, headers=headers, json={"query": query_string_1, "variables": variables})
        api_response = json.loads(r.text)
        max_phase = 0
        drug_name = ''
        #print(api_response)
        for j in api_response['data']['disease']['chembl']['rows']:
            if j['clinicalPhase'] > max_phase:
                max_phase = j['clinicalPhase']
                #if j['drug'] != None:
            if j['clinicalPhase'] >= 2:
                if gene_name != 'TXNRD1':
                    drug_name = j['drug']['name'] + '; ' + drug_name
            #if max_phase >= 2:
            #    break
        gene.loc[gene.loc[:,'gene'] == gene_name,'clinical_phase'] = max_phase
        gene.loc[gene.loc[:,'gene'] == gene_name,'drug_name'] = drug_name
        #gene.loc[gene.loc[:,'gene'] == gene_name,'drug_name'] = drug_name
print(gene)
print(gene.loc[ gene.loc[:,'clinical_phase'] >= 2, :])
gene.to_pickle('data/breast_priority_open_target_prob_rerun_lowthres.pkl')
# Note: json_data will not be serialized by requests
# exactly as it was in the original request.
#data = '{ "query": "{ genes(first:20000) {nodes {name geneCategoriesWithSources {name} } } }" }'
#response = requests.post('https://dgidb.org/api/graphql', headers=headers, data=data)
#'''

#open target query to csv
'''
gene = pd.read_pickle('data/breast_priority_open_target.pkl')
gene.loc[ (gene.loc[:,'sel'] == 'above') & (gene.loc[:,'clinical_phase'] >= 2), 'sel'] = 'exist_drug'
gene.to_csv('R_data/breast_priority_open_target.csv')
#'''

#annotate results
'''
gene = pd.read_csv('R_data/all_gene_prob_thres_lowthres.csv', index_col = 0)
#gene = pd.read_pickle('data/breast_priority_gene_oncoinfo.pkl')
gene2 = pd.read_pickle('data/breast_priority_open_target_prob_rerun_lowthres.pkl')
#print(gene)
#print(gene2)
print(gene2.loc[ gene2.loc[:,'clinical_phase'] >= 2, :])
gene.loc[gene.loc[:,'gene'].isin(gene2.loc[gene2.loc[:,'clinical_phase']>=2,'gene']) & (gene.loc[:,'sel'] == 'above'), 'sel'] = 'exist_drug'
#gene.loc[gene.loc[:,'gene'].isin(gene2.loc[gene2.loc[:,'clinical_phase']>=2,'gene']) & (gene.loc[:,'sel'] == 'above'), 'drug_name'] = 'exist_drug'
gene.loc[:,'drug_name'] = gene2.loc[:,'drug_name'].values.tolist()
gene.loc[ gene.loc[:,'sel'] != 'exist_drug', 'drug_name'] = ''
print(gene.loc[(gene.loc[:,'sel'] == 'exist_drug'),:])
gene.to_csv('R_data/all_gene_prob_thres_drug_lowthres.csv')
#'''

#functional enrichment
'''gene = pd.read_pickle('data/breast_priority_open_target.pkl')
from gprofiler import GProfiler
gp = GProfiler( return_dataframe = True)
for i in ['ER+','HER2+','TNBC']:
    ingene = gene.loc[gene.loc[:,'molecular'] == i,'gene'].tolist()
    print(gene.loc[gene.loc[:,'molecular'] == i,:])
    #result = gp.profile( organism = 'hsapiens', query = ingene, sources = ['GO:MF'], user_threshold = 0.0001, no_evidences = False, no_iea = True)
    result = gp.profile( organism = 'hsapiens', query = ingene, sources = ['GO:MF'],  no_evidences = False, no_iea = True)
    print(result.loc[:,['name','p_value','intersections']])
#'''

#heatmap annotation
#gene = pd.read_csv('R_data/breast_priority_open_target.csv')
gene = pd.read_csv('R_data/all_gene_prob_thres_drug_lowthres.csv')
gene = gene.loc[gene.loc[:,'sel'] == 'above',:]
gene.to_csv('R_data/priority_list_lowthres.csv')
dep = pd.read_pickle('data/ceres_dep.pkl')
sample = pd.read_pickle('data/sample.pkl')
sample = sample.reindex( index = dep.index)
sample = sample.loc[ sample.loc[:,'OncotreeLineage'] == 'Breast',:]
sample = sample.loc[:, ['StrippedCellLineName','LegacyMolecularSubtype','LegacySubSubtype']]
sample.loc['ACH-001065','LegacySubSubtype'] = 'HER2+'
sample.loc['ACH-001419','LegacySubSubtype'] = 'HER2+'
sample.loc['ACH-002179','LegacySubSubtype'] = 'HER2+'
sample.loc['ACH-002399','LegacySubSubtype'] = 'HER2+'
sample.loc['ACH-001820','LegacySubSubtype'] = 'TNBC'
sample.loc[ sample.loc[:,'LegacySubSubtype'].str.contains('HER2pos'),'LegacySubSubtype'] = 'HER2+'
sample.loc[ sample.loc[:,'LegacySubSubtype'] == 'ERneg_HER2neg','LegacySubSubtype'] = 'TNBC'
sample.loc[ sample.loc[:,'LegacySubSubtype'] == 'ERpos_HER2neg','LegacySubSubtype'] = 'ER+'
#sample.loc['ACH-001396','LegacySubSubtype'] = 'ER+'
ingene = list(gene.loc[:,'gene'].unique())
dep = dep.reindex(columns = ingene)
dep = dep.reindex(index = sample.index, columns = ingene)
dep = dep.T
#dep = (dep - dep.mean()) / dep.std()
dep = dep.T
print(sample.columns)

#dep.loc[:,'type'] = sample.loc[:,'LegacySubSubtype']
dep.index = sample.loc[:,'StrippedCellLineName']
sample = sample.loc[:,['StrippedCellLineName','LegacySubSubtype']]
sample.columns = ['name','molecular_type']
anno = pd.DataFrame( False, index = dep.columns, columns = ['ER+','HER2+','TNBC'])
for i in anno.columns:
    anno.loc[gene.loc[gene.loc[:,'molecular'] == i,'gene'],i] = True 
#dep = dep.sort_values(by = 'type')
sample.index.name = 'DepMapID'
anno.columns = ['ER','HER2','TNBC']
print(dep)
print(anno)
print(sample)
sample.to_csv('R_data/priority_heatmap_anno_row.csv')
anno.to_csv('R_data/priority_heatmap_anno_col.csv')
dep.to_csv('R_data/priority_heatmap_data.csv')
print( anno.loc[(anno.loc[:,'ER'] == True) & (anno.loc[:,'HER2'] == True) & (anno.loc[:,'TNBC'] == True),:])
#dep.to_csv('R_data/priority_dep_breast_cell_line.csv')
#dep.to_csv('R_data/priority_dep_breast_cell_line.csv')
#'''
###################
#mutation analysis#
###################
#sample breast
'''sample = pd.read_pickle('data/sample.pkl')
dep = pd.read_pickle('data/ceres_dep.pkl')
sample = sample.reindex( index = dep.index)
sample = sample.loc[ sample.loc[:,'OncotreeLineage'] == 'Breast',:]
sample = sample.loc[:, ['StrippedCellLineName','LegacyMolecularSubtype','LegacySubSubtype']]
sample.loc['ACH-001065','LegacySubSubtype'] = 'HER2+'
sample.loc['ACH-001419','LegacySubSubtype'] = 'HER2+'
sample.loc['ACH-002179','LegacySubSubtype'] = 'HER2+'
sample.loc['ACH-002399','LegacySubSubtype'] = 'HER2+'
sample.loc['ACH-001820','LegacySubSubtype'] = 'TNBC'
sample = sample.loc[ ~sample.loc[:,'LegacySubSubtype'].isna(),:]
sample.loc[ sample.loc[:,'LegacySubSubtype'].str.contains('HER2pos'),'LegacySubSubtype'] = 'HER2+'
sample.loc[ sample.loc[:,'LegacySubSubtype'] == 'ERneg_HER2neg','LegacySubSubtype'] = 'TNBC'
sample.loc[ sample.loc[:,'LegacySubSubtype'] == 'ERpos_HER2neg','LegacySubSubtype'] = 'ER+'
print(sample)
sample.to_pickle('data/sample_breast.pkl')
#'''

#calculate CNV spearman correlation
'''mut_type = 'cnv'
cnv = pd.read_pickle('data/cnv_breast_ori.pkl')
sample = pd.read_pickle('data/sample.pkl')
dep = pd.read_pickle('data/ceres_dep.pkl')
sample = sample.reindex( index = dep.index)
sample = sample.loc[ sample.loc[:,'OncotreeLineage'] == 'Breast',:]
sample = sample.loc[:, ['StrippedCellLineName','LegacyMolecularSubtype','LegacySubSubtype']]
sample.loc['ACH-001065','LegacySubSubtype'] = 'HER2+'
sample.loc['ACH-001419','LegacySubSubtype'] = 'HER2+'
sample.loc['ACH-002179','LegacySubSubtype'] = 'HER2+'
sample.loc['ACH-002399','LegacySubSubtype'] = 'HER2+'
sample.loc['ACH-001820','LegacySubSubtype'] = 'TNBC'
sample = sample.loc[ ~sample.loc[:,'LegacySubSubtype'].isna(),:]
sample.loc[ sample.loc[:,'LegacySubSubtype'].str.contains('HER2pos'),'LegacySubSubtype'] = 'HER2+'
sample.loc[ sample.loc[:,'LegacySubSubtype'] == 'ERneg_HER2neg','LegacySubSubtype'] = 'TNBC'
sample.loc[ sample.loc[:,'LegacySubSubtype'] == 'ERpos_HER2neg','LegacySubSubtype'] = 'ER+'
gene = pd.read_csv('R_data/all_gene_prob_thres.csv')
#gene = gene.loc[gene.loc[:,'sel'] != 'below',:]
#gene = gene.loc[~gene.loc[:,'anno'].isna(),:]
#print(gene)
dep = dep.reindex(index = sample.index)
cnv = cnv.reindex(index = sample.index)
eff = dep.copy()
from scipy import stats
for i in ['ER+','HER2+','TNBC']:
    print(i)
    dep = eff.loc[ sample.loc[:,'LegacySubSubtype'] == i,:]
    #dep = dep.loc[:, dep.columns.isin( gene.loc[ (gene.loc[:,'sel'] == 'above') & (gene.loc[:,'molecular'] == i),'gene'])]
    dep = dep.loc[:, dep.columns.isin( gene.loc[ (gene.loc[:,'molecular'] == i),'gene'])]
    mut_sub = cnv.loc[ sample.loc[:,'LegacySubSubtype'] == i,:]
    p_val = pd.DataFrame( np.nan, index = dep.columns, columns = mut_sub.columns)
    cohen_d = pd.DataFrame( np.nan, index = dep.columns, columns = mut_sub.columns)
    for j in p_val.index:
        for k in p_val.columns:
            cohen_d.loc[j,k],p_val.loc[j,k] =  stats.spearmanr(dep.loc[:,j],mut_sub.loc[:,k],nan_policy = 'omit')
    p_val = p_val.dropna(how = 'all', axis = 1)
    cohen_d = cohen_d.dropna(how = 'all', axis = 1)
    p_val.to_pickle('data/'+mut_type+'_'+i+'_spearman_p_val.pkl')
    print(p_val)
    cohen_d.to_pickle('data/'+mut_type+'_'+i+'_spearman_r.pkl')
    print(cohen_d)
#'''

#mutational data extraction
'''mut = pd.read_csv('data/23Q4/mut_maf.txt', sep= '\t')
print(mut)
print(mut.columns)
mut.to_pickle('data/mut_maf_oncokb.pkl')
#'''
'''mut = pd.read_pickle('data/mut_maf_oncokb.pkl')
cnv = pd.read_pickle('data/cnv_breast.pkl')
mut = mut.loc[ mut.loc[:,'DepMap_ID'].isin(cnv.index),:]
print(mut)
print(mut.loc[:,'Variant_annotation'].value_counts())

mut.to_pickle('data/mut_breast.pkl')
#print(mut.iloc[0,:])
#'''
'''mut = pd.read_pickle('data/mut_breast.pkl')
mut = mut.loc[ mut.loc[:,'Variant_annotation'] != 'silent',:]
print(mut)
mut = mut.loc[:,['Hugo_Symbol','DepMap_ID']]
mut.loc[:,'mut'] = 1
mut = mut.drop_duplicates()
mut = mut.pivot( columns = 'Hugo_Symbol', index = 'DepMap_ID', values = 'mut')
mut = mut.fillna(0)
print(mut)
mut.to_pickle('data/mut_onehot.pkl')
#'''

#CNV depmap descrete value extraction
'''
cnv = pd.read_pickle('download/ccle_cnv_curate.pkl')
cnv = cnv.T
sample = pd.read_pickle('data/sample.pkl')
dep = pd.read_pickle('data/ceres_dep.pkl')
sample = sample.reindex( index = dep.index)
sample = sample.loc[ sample.loc[:,'OncotreeLineage'] == 'Breast',:]
sample = sample.loc[:, ['StrippedCellLineName','LegacyMolecularSubtype','LegacySubSubtype']]
sample.loc['ACH-001065','LegacySubSubtype'] = 'HER2_pos'
sample.loc['ACH-001419','LegacySubSubtype'] = 'HER2_pos'
sample.loc['ACH-002179','LegacySubSubtype'] = 'HER2_pos'
sample.loc['ACH-002399','LegacySubSubtype'] = 'HER2_pos'
sample.loc['ACH-001820','LegacySubSubtype'] = 'TNBC'
sample = sample.loc[ ~sample.loc[:,'LegacySubSubtype'].isna(),:]
sample.loc[ sample.loc[:,'LegacySubSubtype'].str.contains('HER2pos'),'LegacySubSubtype'] = 'HER2_pos'
sample.loc[ sample.loc[:,'LegacySubSubtype'] == 'ERneg_HER2neg','LegacySubSubtype'] = 'TNBC'
sample.loc[ sample.loc[:,'LegacySubSubtype'] == 'ERpos_HER2neg','LegacySubSubtype'] = 'ER_pos'
cnv = cnv.loc[ cnv.index.isin(sample.loc[:,'StrippedCellLineName']),:]
print(cnv)
print(cnv.loc[:,'ERBB2'])
sample.loc[:,'depmap'] = sample.index
sample.index = sample.loc[:,'StrippedCellLineName']
sample = sample.reindex(index = cnv.index)
cnv.index = sample.loc[:,'depmap']
cnv =cnv.astype('float')
print(cnv)
cnv.to_pickle('data/ccle_cnv.pkl')
#'''

#Mann whitney U test and cohen d for association
#according to subtype
'''mut_type = 'mut'
mut = pd.read_pickle('data/mut_onehot.pkl')
sample = pd.read_pickle('data/sample.pkl')
dep = pd.read_pickle('data/ceres_dep.pkl')
sample = sample.reindex( index = dep.index)
sample = sample.loc[ sample.loc[:,'OncotreeLineage'] == 'Breast',:]
sample = sample.loc[:, ['StrippedCellLineName','LegacyMolecularSubtype','LegacySubSubtype']]
sample.loc['ACH-001065','LegacySubSubtype'] = 'HER2+'
sample.loc['ACH-001419','LegacySubSubtype'] = 'HER2+'
sample.loc['ACH-002179','LegacySubSubtype'] = 'HER2+'
sample.loc['ACH-002399','LegacySubSubtype'] = 'HER2+'
sample.loc['ACH-001820','LegacySubSubtype'] = 'TNBC'
sample = sample.loc[ ~sample.loc[:,'LegacySubSubtype'].isna(),:]
sample.loc[ sample.loc[:,'LegacySubSubtype'].str.contains('HER2pos'),'LegacySubSubtype'] = 'HER2+'
sample.loc[ sample.loc[:,'LegacySubSubtype'] == 'ERneg_HER2neg','LegacySubSubtype'] = 'TNBC'
sample.loc[ sample.loc[:,'LegacySubSubtype'] == 'ERpos_HER2neg','LegacySubSubtype'] = 'ER+'
gene = pd.read_csv('R_data/all_gene_prob_thres_lowthres.csv')
#gene.loc[ (gene.loc[:,'sel'] == 'above') & (gene.loc[:,'clinical_phase'] >= 2), 'sel'] = 'exist_drug'
print(gene)
#gene = gene.loc[gene.loc[:,'sel'] != 'below',:]
#gene = gene.loc[~gene.loc[:,'anno'].isna(),:]
#print(gene)
dep = dep.reindex(index = sample.index)
eff = dep.copy()
from scipy import stats
mut = pd.read_pickle('data/mut_onehot.pkl')
mut_mean = mut.mean(axis = 0)
mut = mut.loc[:,mut_mean>0.1]
print(mut)
for i in ['ER+','HER2+','TNBC']:
    print(i)
    dep = eff.loc[ sample.loc[:,'LegacySubSubtype'] == i,:]
    #dep = dep.loc[:, dep.columns.isin( gene.loc[ (gene.loc[:,'sel'] == 'above') & (gene.loc[:,'molecular'] == i),'gene'])]
    dep = dep.loc[:, dep.columns.isin( gene.loc[ (gene.loc[:,'molecular'] == i),'gene'])]
    mut_sub = mut.loc[ sample.loc[:,'LegacySubSubtype'] == i,:]
    cut_off = 2
    if False:
        if len(mut_sub) * 0.1 < 2:
            cut_off = 2
        else:
            cut_off = len(mut_sub)*0.1
    mut_sub = mut_sub.loc[:, mut_sub.sum() >= cut_off]
    p_val = pd.DataFrame( np.nan, index = dep.columns, columns = mut_sub.columns)
    cohen_d = pd.DataFrame( np.nan, index = dep.columns, columns = mut_sub.columns)
    for j in p_val.index:
        for k in p_val.columns:
            a = dep.loc[ mut_sub.loc[:,k] == 1,j]
            b = dep.loc[ mut_sub.loc[:,k] == 0,j]
            #c = dep.loc[ mut_sub.loc[:,k] == -1,j]
            if len(b) == 0 or len(a) == 0:
                continue
            t, p = stats.mannwhitneyu( a, b)
            p_val.loc[j,k] = p
            s = np.sqrt(((len(a) - 1)*(a.std())**2 + (len(b)-1)*(b.std())**2) / (len(a)+len(b)-2))
            d = (a.mean() - b.mean())/s
            cohen_d.loc[j,k] = d
    p_val = p_val.dropna(how = 'all', axis = 1)
    cohen_d = cohen_d.dropna(how = 'all', axis = 1)
    p_val.to_pickle('data/'+mut_type+'_'+i+'_ttest_p_val_prob_lowthres.pkl')
    print(p_val)
    cohen_d.to_pickle('data/'+mut_type+'_'+i+'_ttest_cohen_d_prob_lowthres.pkl')
    print(cohen_d)
#'''

#CNV
'''mut_type = 'cnv'
#mut = pd.read_pickle('data/23Q4/oncokb_'+mut_type+'.pkl')
mut = pd.read_pickle('data/cnv_breast_derive.pkl')
#mut = pd.read_pickle('data/ccle_cnv.pkl')
print(mut)
sample = pd.read_pickle('data/sample.pkl')
dep = pd.read_pickle('data/ceres_dep.pkl')
sample = sample.reindex( index = dep.index)
mut = mut.reindex(columns = dep.columns)
sample = sample.loc[ sample.loc[:,'OncotreeLineage'] == 'Breast',:]
sample = sample.loc[:, ['StrippedCellLineName','LegacyMolecularSubtype','LegacySubSubtype']]
sample.loc['ACH-001065','LegacySubSubtype'] = 'HER2+'
sample.loc['ACH-001419','LegacySubSubtype'] = 'HER2+'
sample.loc['ACH-002179','LegacySubSubtype'] = 'HER2+'
sample.loc['ACH-002399','LegacySubSubtype'] = 'HER2+'
sample.loc['ACH-001820','LegacySubSubtype'] = 'TNBC'
sample = sample.loc[ ~sample.loc[:,'LegacySubSubtype'].isna(),:]
sample.loc[ sample.loc[:,'LegacySubSubtype'].str.contains('HER2pos'),'LegacySubSubtype'] = 'HER2+'
sample.loc[ sample.loc[:,'LegacySubSubtype'] == 'ERneg_HER2neg','LegacySubSubtype'] = 'TNBC'
sample.loc[ sample.loc[:,'LegacySubSubtype'] == 'ERpos_HER2neg','LegacySubSubtype'] = 'ER+'
sample = sample.reindex(index = mut.index)
gene = pd.read_csv('R_data/all_gene_prob_thres_lowthres.csv')
#gene.loc[ (gene.loc[:,'sel'] == 'above') & (gene.loc[:,'clinical_phase'] >= 2), 'sel'] = 'exist_drug'
print(gene)
#gene = gene.loc[gene.loc[:,'sel'] != 'below',:]
#gene = gene.loc[~gene.loc[:,'anno'].isna(),:]
#print(gene)
dep = dep.reindex(index = sample.index)
eff = dep.copy()
from scipy import stats
#mut = pd.read_pickle('data/mut_onehot.pkl')
mut_mean = mut.copy()
mut_mean[mut_mean!=0] = 1
mut_mean = mut.mean(axis = 0)
mut = mut.loc[:,mut_mean>0.1]
print(mut)
#mut.to_pickle('data/cnv_breast_sel.pkl')
if True:
    for i in ['ER+','HER2+','TNBC']:
        print(i)
        dep = eff.loc[ sample.loc[:,'LegacySubSubtype'] == i,:]
        #dep = dep.loc[:, dep.columns.isin( gene.loc[ (gene.loc[:,'sel'] == 'above') & (gene.loc[:,'molecular'] == i),'gene'])]
        dep = dep.loc[:, dep.columns.isin( gene.loc[ (gene.loc[:,'molecular'] == i),'gene'])]
        mut_sub = mut.loc[ sample.loc[:,'LegacySubSubtype'] == i,:]
        if False:
            cut_off = 2
            if False:
                if len(mut_sub) * 0.1 < 2:
                    cut_off = 2
                else:
                    cut_off = len(mut_sub)*0.1
        #mut_sub = mut_sub.loc[:, mut_sub.sum() >= cut_off]
        print(mut_sub)
        p_val = pd.DataFrame( np.nan, index = dep.columns, columns = mut_sub.columns)
        cohen_d = pd.DataFrame( np.nan, index = dep.columns, columns = mut_sub.columns)
        for j in p_val.index:
            for k in p_val.columns:
                a = dep.loc[ mut_sub.loc[:,k] == 1,j]
                b = dep.loc[ mut_sub.loc[:,k] == 0,j]
                #c = dep.loc[ mut_sub.loc[:,k] == -1,j]
                if len(b) == 0 or len(a) == 0:
                    continue
                if len(a)/(len(a)+len(b)) <= 0.1:
                    continue
                if False:
                    if len(a) == 0:
                        t, p = stats.kruskal( b, c)
                    elif len(b) == 0:
                        t, p = stats.kruskal( a, c)
                    elif len(c) == 0:
                        t, p = stats.kruskal( a, b)
                    else:
                        t, p = stats.kruskal( a, b, c)
                t, p = stats.mannwhitneyu( a, b)
                p_val.loc[j,k] = p
                s = np.sqrt(((len(a) - 1)*(a.std())**2 + (len(b)-1)*(b.std())**2) / (len(a)+len(b)-2))
                d = (a.mean() - b.mean())/s
                cohen_d.loc[j,k] = d
        p_val = p_val.dropna(how = 'all', axis = 1)
        cohen_d = cohen_d.dropna(how = 'all', axis = 1)
        p_val.to_pickle('data/'+mut_type+'_'+i+'_ttest_p_val_prob_derive_amp_lowthres.pkl')
        print(p_val)
        cohen_d.to_pickle('data/'+mut_type+'_'+i+'_ttest_cohen_d_prob_derive_amp_lowthres.pkl')
        print(cohen_d)
if False:
    for i in ['ER+','HER2+','TNBC']:
        print(i)
        dep = eff.loc[ sample.loc[:,'LegacySubSubtype'] == i,:]
        #dep = dep.loc[:, dep.columns.isin( gene.loc[ (gene.loc[:,'sel'] == 'above') & (gene.loc[:,'molecular'] == i),'gene'])]
        dep = dep.loc[:, dep.columns.isin( gene.loc[ (gene.loc[:,'molecular'] == i),'gene'])]
        mut_sub = mut.loc[ sample.loc[:,'LegacySubSubtype'] == i,:]
        if False:
            cut_off = 2
            if False:
                if len(mut_sub) * 0.1 < 2:
                    cut_off = 2
                else:
                    cut_off = len(mut_sub)*0.1
        #mut_sub = mut_sub.loc[:, mut_sub.sum() >= cut_off]
        print(mut_sub)
        p_val = pd.DataFrame( np.nan, index = dep.columns, columns = mut_sub.columns)
        cohen_d = pd.DataFrame( np.nan, index = dep.columns, columns = mut_sub.columns)
        for j in p_val.index:
            for k in p_val.columns:
                a = dep.loc[ mut_sub.loc[:,k] == -1,j]
                b = dep.loc[ mut_sub.loc[:,k] == 0,j]
                #c = dep.loc[ mut_sub.loc[:,k] == -1,j]
                if len(b) == 0 or len(a) == 0:
                    continue
                if len(a)+len(b)-2 <= 0:
                    continue
                if len(a)/(len(a)+len(b)) <= 0.1:
                    continue
                if False:
                    if len(a) == 0:
                        t, p = stats.kruskal( b, c)
                    elif len(b) == 0:
                        t, p = stats.kruskal( a, c)
                    elif len(c) == 0:
                        t, p = stats.kruskal( a, b)
                    else:
                        t, p = stats.kruskal( a, b, c)
                t, p = stats.mannwhitneyu( a, b)
                p_val.loc[j,k] = p
                s = np.sqrt(((len(a) - 1)*(a.std())**2 + (len(b)-1)*(b.std())**2) / (len(a)+len(b)-2))
                d = (a.mean() - b.mean())/s
                cohen_d.loc[j,k] = d
        p_val = p_val.dropna(how = 'all', axis = 1)
        cohen_d = cohen_d.dropna(how = 'all', axis = 1)
        p_val.to_pickle('data/'+mut_type+'_'+i+'_ttest_p_val_prob_derive_del_lowthres.pkl')
        print(p_val)
        cohen_d.to_pickle('data/'+mut_type+'_'+i+'_ttest_cohen_d_prob_derive_del_lowthres.pkl')
        print(cohen_d)
#'''

#mutation association result.
'''import statsmodels.stats.multitest as stat
gene = pd.read_csv('R_data/all_gene_prob_thres_lowthres.csv')
output = []
for i in ['ER+','HER2+','TNBC']:
    print(i)
    tmp = gene.loc[ gene.loc[:,'molecular'] == i,:]
    tmp = tmp.loc[ tmp.loc[:,'anno'] != 'below',:]
    #tmp = tmp.loc[ ~tmp.loc[:,'anno'].isna(),:]
    #print(tmp)

    #p_val = pd.read_pickle('data/cnv_'+i+'_spearman_p_val.pkl')
    #cohen_d = pd.read_pickle('data/cnv_'+i+'_spearman_r.pkl')
    p_val = pd.read_pickle('data/mut_'+i+'_ttest_p_val_prob_lowthres.pkl')
    print(p_val)
    cohen_d = pd.read_pickle('data/mut_'+i+'_ttest_cohen_d_prob_lowthres.pkl')
    p_val = p_val.loc[ p_val.index.isin( tmp.loc[:,'gene']),:]
    cohen_d = cohen_d.reindex(index = p_val.index, columns = p_val.columns)
    #print(p_val.min(axis = 1))
    for j in p_val.index:
        tmp = p_val.loc[j,:].T
        reject, correct, _ , _ = stat.multipletests(tmp, method = 'fdr_bh', alpha = 0.05)
        #result.loc[:,'fdr'] = correct
        #result.loc[:,'sel'] = 'not_sig'
        p_val.loc[j,:] = correct
    p_val = p_val.melt(ignore_index=False)
    p_val = p_val.reset_index()
    cohen_d = cohen_d.melt(ignore_index=False)
    cohen_d = cohen_d.reset_index()
    p_val.loc[:,'cohen_d'] = cohen_d.loc[:,'value']
    p_val.loc[:,'molecular'] = i
    p_val.loc[:,'assign'] = 'not_sig'
    p_val.loc[(p_val.loc[:,'value'] < 0.1) & (p_val.loc[:,'cohen_d'].abs() > 1),'assign'] = 'sig'
    p_val.loc[:,'anno'] = p_val.loc[:,'index'] + '_' + p_val.loc[:,'Hugo_Symbol']
    p_val.loc[:,'fdr'] = -np.log10(p_val.loc[:,'value'])
    #p_val.loc[ ~p_val.index.isin(p_val.loc[p_val.loc[:,'assign'] == 'sig','fdr'].abs().nlargest(5).index),'anno'] = ''
    p_val.loc[ ~p_val.index.isin(p_val.loc[ (p_val.loc[:,'assign'] == 'sig') & (p_val.loc[:,'cohen_d'] > 1), ['cohen_d','fdr']].abs().nlargest(5, ['cohen_d', 'fdr']).index),'anno'] = ''
    p_val.loc[ p_val.loc[:,'assign'] != 'sig','anno'] = ''
    #p_val = p_val.loc[ (p_val.loc[:,'value'] < 0.1),:]
    #p_val = p_val.loc[ (p_val.loc[:,'value'] < 0.1) & (p_val.loc[:,'cohen_d'].abs() > 1),:]
    print(p_val.loc[:,'assign'].value_counts())
    print(p_val.loc[p_val.loc[:,'anno'] != '',:])
    output.append(p_val)
    #print(cohen_d)
output = pd.concat(output, axis = 0)
#print(output)
output.to_csv('R_data/mutation_association_lowthres.txt', sep = '\t')
#'''
'''
output = pd.read_csv('R_data/mutation_association_lowthres.txt', sep = '\t', index_col = 0)
output = output.rename( columns = { 'index':'g1_dependency', 'variable':'g2_mutation', 'value': 'FDR'})
output = output.loc[ output.loc[:,'assign'] == 'sig',:]
output = output.drop(columns = ['assign','anno'])
output = output.sort_values( by = ['molecular','cohen_d','FDR'], ascending = [True,False,True])
print(output)
output.to_csv('R_data/mutation_association_lowthres_table.txt', sep = '\t')
#'''

#CNV result
'''
import matplotlib.backends.backend_pdf
import statsmodels.stats.multitest as stat
import seaborn as sns
gene = pd.read_csv('R_data/all_gene_prob_thres.csv')
sample = pd.read_pickle('data/sample_breast.pkl')
#mut = pd.read_pickle('data/ccle_cnv.pkl')
mut = pd.read_pickle('data/cnv_breast_derive.pkl')
sample = sample.reindex(mut.index)
pdf = matplotlib.backends.backend_pdf.PdfPages('figure/test_box.pdf')
eff = pd.read_pickle('data/ceres_dep.pkl')
output = []
for i in ['ER+','HER2+','TNBC']:
    print(i)
    tmp = gene.loc[ gene.loc[:,'molecular'] == i,:]
    tmp = tmp.loc[ tmp.loc[:,'anno'] != 'below',:]
    #tmp = tmp.loc[ ~tmp.loc[:,'anno'].isna(),:]
    #print(tmp)

    mut_mean = mut.copy()
    mut_mean = mut_mean.loc[ sample.loc[:,'LegacySubSubtype'] == i,:]
    efftmp = eff.reindex( index = mut_mean.index)
    mut_mean[mut_mean!=0] = 1
    p_val_r = pd.read_pickle('data/cnv_'+i+'_spearman_p_val.pkl')
    cohen_d = pd.read_pickle('data/cnv_'+i+'_spearman_r.pkl')
    #p_val = pd.read_pickle('data/cnv_'+i+'_ttest_p_val_prob_ccle.pkl')
    p_val = pd.read_pickle('data/cnv_'+i+'_ttest_p_val_prob_derive.pkl')
    print(p_val)
    mut_mean = mut_mean.reindex(columns = p_val.columns)
    mut_mean = mut_mean.mean(axis = 0)
    mut_mean = mut_mean[mut_mean>0.25]
    print(len(mut_mean))
    p_val = p_val.loc[:,p_val.columns.isin(mut_mean.index)]
    p_val = p_val.loc[ p_val.index.isin( tmp.loc[:,'gene']),:]
    p_val_r = p_val_r.reindex( index = p_val.index, columns = p_val.columns)
    p_val_r = p_val_r.dropna( how = 'all', axis = 1)
    #print(p_val)
    #print(p_val.min(axis = 1))
    for j in p_val.index:
        tmp = p_val.loc[j,:].T
        reject, correct, _ , _ = stat.multipletests(tmp, method = 'fdr_bh', alpha = 0.05)
        p_val.loc[j,:] = correct
        tmp = p_val_r.loc[j,:].T
        reject, correct, _ , _ = stat.multipletests(tmp, method = 'fdr_bh', alpha = 0.05)
        p_val_r.loc[j,:] = correct
    cohen_d = cohen_d.reindex( index = p_val.index, columns = p_val.columns)
    p_val_r = p_val_r.reindex( index = p_val.index, columns = p_val.columns)
    #print(cohen_d)
    p_val = p_val.melt(ignore_index=False)
    p_val = p_val.reset_index()
    cohen_d = cohen_d.melt(ignore_index=False)
    cohen_d = cohen_d.reset_index()
    p_val_r = p_val_r.melt(ignore_index=False)
    p_val_r = p_val_r.reset_index()
    p_val.loc[:,'cohen_d'] = cohen_d.loc[:,'value']
    p_val.loc[:,'r_p'] = p_val_r.loc[:,'value']
    #p_val = p_val.loc[ (p_val.loc[:,'value'] < 0.1),:]
    #p_val = p_val.loc[ (p_val.loc[:,'value'] < 0.1) & (p_val.loc[:,'cohen_d'].abs() > 0.5),:]
    #p_val = p_val.loc[ (p_val.loc[:,'value'] < 0.1) & (p_val.loc[:,'cohen_d'].abs() > 0.5) & (p_val.loc[:,'r_p'] < 0.1),:]
    p_val.loc[:,'molecular'] = i
    p_val.loc[:,'assign'] = 'not_sig'
    p_val.loc[(p_val.loc[:,'value'] < 0.1) & (p_val.loc[:,'cohen_d'].abs() > 0.5),'assign'] = 'sig'
    p_val.loc[:,'anno'] = p_val.loc[:,'index'] + '_' + p_val.loc[:,'variable']
    p_val.loc[:,'fdr'] = -np.log10(p_val.loc[:,'value'])
    #p_val.loc[ ~p_val.index.isin(p_val.loc[ (p_val.loc[:,'assign'] == 'sig') & (p_val.loc[:,'cohen_d'] > 0.5), ['fdr','cohen_d']].abs().nlargest(5, ['fdr', 'cohen_d']).index),'anno'] = ''
    p_val.loc[ ~p_val.index.isin(p_val.loc[ (p_val.loc[:,'assign'] == 'sig'), ['fdr']].abs().nlargest(5, ['fdr']).index),'anno'] = ''
    p_val.loc[ p_val.loc[:,'assign'] != 'sig','anno'] = ''
    print(p_val)
    #print(cohen_d)
    if len(p_val) > 0:
        fig = plt.figure()
        #fig, ax = plt.subplots()
        tmp = pd.concat([ efftmp.loc[:,p_val.loc[ p_val.index[0],'index']], mut.loc[:,p_val.loc[p_val.index[0],'variable']]], axis = 1)
        tmp.columns = ['g1','g2']
        tmp.loc[tmp.loc[:,'g2'] == 0,'g2'] = 0
        sns.boxplot(tmp,x = 'g2',y = 'g1')
        #plt.hist(tmp[i].sum(axis=1)/len(tmp[i].columns), bins = 20)
        plt.title(p_val.loc[p_val.index[0],'index']+'_'+p_val.loc[p_val.index[0],'variable']+'_'+i)
        pdf.savefig(fig)
        plt.close()
    output.append(p_val)
    #print(cohen_d)
output = pd.concat(output, axis = 0)
print(output)
output.to_csv('R_data/cnv_association.txt', sep = '\t')
pdf.close()
#'''

#CNV result
'''import matplotlib.backends.backend_pdf
import statsmodels.stats.multitest as stat
import seaborn as sns
gene = pd.read_csv('R_data/all_gene_prob_thres_lowthres.csv')
sample = pd.read_pickle('data/sample_breast.pkl')
#mut = pd.read_pickle('data/ccle_cnv.pkl')
mut = pd.read_pickle('data/cnv_breast_derive.pkl')
sample = sample.reindex(mut.index)
eff = pd.read_pickle('data/ceres_dep.pkl')
output = []
selgene = [ ['TFAP2C','BHLHE23_amp', 'GATA5_amp', 'GMEB2_amp', 'TAF4_amp', 'ZNF512B_amp'], ['NDUFAF4','ADCK5_amp','CYC1_amp','CYP11B1_amp','CYP11B2_amp'], ['RNF115_amp','PIAS3_amp','HJV_amp','ITGA10_amp','UBE2Q1']]
selgene_i = 0
for i in ['ER+','HER2+','TNBC']:
    print(i)
    tmp = gene.loc[ gene.loc[:,'molecular'] == i,:]
    tmp = tmp.loc[ tmp.loc[:,'anno'] != 'below',:]
    #tmp = tmp.loc[ tmp.loc[:,'sel'] != 'below',:]
    #tmp = tmp.loc[ ~tmp.loc[:,'anno'].isna(),:]
    print(tmp)
    mut_mean = mut.copy()
    mut_mean = mut_mean.loc[ sample.loc[:,'LegacySubSubtype'] == i,:]
    efftmp = eff.reindex( index = mut_mean.index)
    mut_mean[mut_mean!=0] = 1
    p_val_a = pd.read_pickle('data/cnv_'+i+'_ttest_p_val_prob_derive_amp_lowthres.pkl')
    p_val_d = pd.read_pickle('data/cnv_'+i+'_ttest_p_val_prob_derive_del_lowthres.pkl')
    cohen_da = pd.read_pickle('data/cnv_'+i+'_ttest_cohen_d_prob_derive_amp_lowthres.pkl')
    cohen_dd = pd.read_pickle('data/cnv_'+i+'_ttest_cohen_d_prob_derive_del_lowthres.pkl')
    mut_mean = mut_mean.mean(axis = 0)
    mut_mean = mut_mean[mut_mean>0.25]
    p_val_a = p_val_a.reindex(columns = cohen_da.columns)
    p_val_d = p_val_d.reindex(columns = cohen_dd.columns)
    p_val_a = p_val_a.loc[:,p_val_a.columns.isin(mut_mean.index)]
    p_val_a = p_val_a.loc[ p_val_a.index.isin( tmp.loc[:,'gene']),:]
    p_val_d = p_val_d.loc[:,p_val_d.columns.isin(mut_mean.index)]
    p_val_d = p_val_d.loc[ p_val_d.index.isin( tmp.loc[:,'gene']),:]
    cohen_da = cohen_da.reindex( index = p_val_a.index, columns = p_val_a.columns)
    cohen_dd = cohen_dd.reindex( index = p_val_d.index, columns = p_val_d.columns)
    colname = [i+'_amp' for i in p_val_a.columns]
    p_val_a.columns = colname
    cohen_da.columns = colname
    colname = [i+'_amp' for i in p_val_d.columns]
    p_val_d.columns = colname
    cohen_dd.columns = colname
    p_val = pd.concat([p_val_a, p_val_d], axis = 1)
    cohen_d = pd.concat([cohen_da, cohen_dd], axis = 1)
    for j in p_val.index:
        tmp = p_val.loc[j,:].T
        reject, correct, _ , _ = stat.multipletests(tmp, method = 'fdr_bh', alpha = 0.05)
        p_val.loc[j,:] = correct
    p_val = p_val.melt(ignore_index=False)
    p_val = p_val.reset_index()
    cohen_d = cohen_d.melt(ignore_index=False)
    cohen_d = cohen_d.reset_index()
    p_val.loc[:,'cohen_d'] = cohen_d.loc[:,'value']
    p_val.loc[:,'molecular'] = i
    p_val.loc[:,'assign'] = 'not_sig'
    p_val.loc[(p_val.loc[:,'value'] < 0.1) & (p_val.loc[:,'cohen_d'].abs() > 1),'assign'] = 'sig'
    p_val.loc[:,'anno'] = p_val.loc[:,'index'] + '_' + p_val.loc[:,'variable']
    p_val.loc[:,'fdr'] = -np.log10(p_val.loc[:,'value'])
    p_val.loc[ ~(p_val.loc[:,'index'].isin(selgene[selgene_i]) & p_val.loc[:,'variable'].isin(selgene[selgene_i])),'anno'] = ''
    #p_val.loc[ ~p_val.index.isin(p_val.loc[ (p_val.loc[:,'assign'] == 'sig') & (p_val.loc[:,'cohen_d'] > 1), ['cohen_d','fdr']].abs().nlargest(5, ['cohen_d', 'fdr']).index),'anno'] = ''
    #p_val.loc[ ~p_val.index.isin(p_val.loc[ (p_val.loc[:,'assign'] == 'sig'), ['fdr']].abs().nlargest(5, ['fdr']).index),'anno'] = ''
    p_val.loc[ p_val.loc[:,'assign'] != 'sig','anno'] = ''
    #print(p_val.loc[:,'assign'].value_counts())
    #print(p_val.loc[p_val.loc[:,'assign'] == 'sig',:])
    print(p_val.loc[p_val.loc[:,'anno'] != '',:])
    #print(cohen_d)
    output.append(p_val)
    #print(cohen_d)
    selgene_i = selgene_i + 1
output = pd.concat(output, axis = 0)
#print(output)
output.to_csv('R_data/cnv_association_derive_lowthres.txt', sep = '\t')
output = output.loc[output.loc[:,'assign'] == 'sig',:]
#output.to_csv('R_data/cnv_association_derive_short_lowthres.txt', sep = '\t')
#'''

'''
output = pd.read_csv('R_data/cnv_association_derive_short_lowthres.txt', sep = '\t', index_col = 0)
output = output.rename( columns = { 'index':'g1_dependency', 'variable':'g2_CNV', 'value': 'FDR'})
output = output.loc[ output.loc[:,'assign'] == 'sig',:]
output = output.drop(columns = ['assign','anno'])
output = output.sort_values( by = ['molecular','cohen_d','FDR'], ascending = [True,False,True])
print(output)
output.to_csv('R_data/cnv_association_lowthres_table.txt', sep = '\t')
#'''

'''
#CNV result
import matplotlib.backends.backend_pdf
import statsmodels.stats.multitest as stat
import seaborn as sns
gene = pd.read_csv('R_data/all_gene_prob_thres.csv')
sample = pd.read_pickle('data/sample_breast.pkl')
#mut = pd.read_pickle('data/ccle_cnv.pkl')
mut = pd.read_pickle('data/cnv_breast_derive.pkl')
sample = sample.reindex(mut.index)
eff = pd.read_pickle('data/ceres_dep.pkl')
output = []
for i in ['ER+','HER2+','TNBC']:
    print(i)
    tmp = gene.loc[ gene.loc[:,'molecular'] == i,:]
    tmp = tmp.loc[ tmp.loc[:,'anno'] != 'below',:]
    #tmp = tmp.loc[ ~tmp.loc[:,'anno'].isna(),:]
    #print(tmp)

    mut_mean = mut.copy()
    mut_mean = mut_mean.loc[ sample.loc[:,'LegacySubSubtype'] == i,:]
    efftmp = eff.reindex( index = mut_mean.index)
    mut_mean[mut_mean!=0] = 1
    p_val_a = pd.read_pickle('data/cnv_'+i+'_ttest_p_val_prob_derive_amp.pkl')
    p_val_d = pd.read_pickle('data/cnv_'+i+'_ttest_p_val_prob_derive_del.pkl')
    cohen_da = pd.read_pickle('data/cnv_'+i+'_ttest_cohen_d_prob_derive_amp.pkl')
    cohen_dd = pd.read_pickle('data/cnv_'+i+'_ttest_cohen_d_prob_derive_del.pkl')
    mut_mean = mut_mean.mean(axis = 0)
    mut_mean = mut_mean[mut_mean>0.3]
    p_val_a = p_val_a.reindex(columns = cohen_da.columns)
    p_val_d = p_val_d.reindex(columns = cohen_dd.columns)
    p_val_a = p_val_a.loc[:,p_val_a.columns.isin(mut_mean.index)]
    p_val_a = p_val_a.loc[ p_val_a.index.isin( tmp.loc[:,'gene']),:]
    p_val_d = p_val_d.loc[:,p_val_d.columns.isin(mut_mean.index)]
    p_val_d = p_val_d.loc[ p_val_d.index.isin( tmp.loc[:,'gene']),:]
    cohen_da = cohen_da.reindex( index = p_val_a.index, columns = p_val_a.columns)
    cohen_dd = cohen_dd.reindex( index = p_val_d.index, columns = p_val_d.columns)
    colname = [i+'_amp' for i in p_val_a.columns]
    p_val_a.columns = colname
    cohen_da.columns = colname
    colname = [i+'_amp' for i in p_val_d.columns]
    p_val_d.columns = colname
    cohen_dd.columns = colname
    p_val = pd.concat([p_val_a, p_val_d], axis = 1)
    cohen_d = pd.concat([cohen_da, cohen_dd], axis = 1)
    p_val = p_val.melt(ignore_index=False)
    p_val = p_val.reset_index()
    cohen_d = cohen_d.melt(ignore_index=False)
    cohen_d = cohen_d.reset_index()
    p_val.loc[:,'cohen_d'] = cohen_d.loc[:,'value']
    p_val.loc[:,'molecular'] = i
    p_val.loc[:,'assign'] = 'not_sig'
    tmp = p_val.loc[:,'value']
    reject, correct, _ , _ = stat.multipletests(tmp, method = 'fdr_bh', alpha = 0.05)
    p_val.loc[:,'value'] = correct
    p_val.loc[(p_val.loc[:,'value'] < 0.2) & (p_val.loc[:,'cohen_d'].abs() > 1),'assign'] = 'sig'
    p_val.loc[:,'anno'] = p_val.loc[:,'index'] + '_' + p_val.loc[:,'variable']
    p_val.loc[:,'fdr'] = -np.log10(p_val.loc[:,'value'])
    p_val.loc[ ~p_val.index.isin(p_val.loc[ (p_val.loc[:,'assign'] == 'sig') & (p_val.loc[:,'cohen_d'] > 1), ['cohen_d','fdr']].abs().nlargest(5, ['fdr', 'cohen_d']).index),'anno'] = ''
    #p_val.loc[ ~p_val.index.isin(p_val.loc[ (p_val.loc[:,'assign'] == 'sig'), ['fdr']].abs().nlargest(5, ['fdr']).index),'anno'] = ''
    p_val.loc[ p_val.loc[:,'assign'] != 'sig','anno'] = ''
    print(p_val.loc[:,'assign'].value_counts())
    print(p_val.loc[p_val.loc[:,'assign'] == 'sig',:])
    print(p_val.loc[p_val.loc[:,'anno'] != '',:])
    #print(cohen_d)
    output.append(p_val)
    #print(cohen_d)
output = pd.concat(output, axis = 0)
print(output)
#output.to_csv('R_data/cnv_association_derive.txt', sep = '\t')
#'''
