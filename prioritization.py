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
'''#f = open("data/cnv_22Q4", "r")
f = open("data_2023/OmicsCNGene.csv", "r")
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
exp.to_pickle('data_2023/cnv_original.pkl')
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

###############################
#extracting prioritizing genes#
###############################

#extracting priotorizing genes according to subtype
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
#print(sample.loc[:,'LegacySubSubtype'].value_counts())
eff = pd.read_pickle('data/ceres_dep.pkl')
eff = eff.loc[:,eff.isna().sum() < len(eff.index) * 0.2]
eff = eff.reindex(index = sample.index)
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

er_gene = pd.DataFrame( np.nan, index = er.columns, columns = ['perc','median_eff','assignment'])
her_gene = pd.DataFrame( np.nan, index = her.columns, columns = ['perc','median_eff','assignment'])
tnbc_gene = pd.DataFrame( np.nan, index = tnbc.columns, columns = ['perc','median_eff','assignment'])
ceres_control = pd.read_csv('data/23Q4/common_ess_2021.csv', sep = ',',  index_col = 0)
achilles_control = pd.read_csv('data/23Q4/AchillesCommonEssentialControls.csv', sep = ' ', index_col = 0)
er_gene.loc[:,'perc'] = (er >= 0.5).sum()/len(er.index)
her_gene.loc[:,'perc'] = (her >= 0.5).sum()/len(her.index)
tnbc_gene.loc[:,'perc'] = (tnbc >= 0.5).sum()/len(tnbc.index)
er_gene.loc[:,'median_eff'] =  dep.reindex(index = er.index, columns = er_gene.index).median()
her_gene.loc[:,'median_eff'] =  dep.reindex(index = her.index, columns = her_gene.index).median()
tnbc_gene.loc[:,'median_eff'] =  dep.reindex(index = tnbc.index, columns = tnbc_gene.index).median()
er_gene.loc[ er_gene.loc[:,'perc']*8 >= 2, 'assignment'] = 'priority'
her_gene.loc[ her_gene.loc[:,'perc'] >= 0.1, 'assignment'] = 'priority'
tnbc_gene.loc[ tnbc_gene.loc[:,'perc'] >= 0.1, 'assignment'] = 'priority'
er_gene.loc[ er_gene.loc[:,'perc']*8 < 2, 'assignment'] = 'lower 10%'
her_gene.loc[ her_gene.loc[:,'perc'] < 0.1, 'assignment'] = 'lower 10%'
tnbc_gene.loc[ tnbc_gene.loc[:,'perc'] < 0.1, 'assignment'] = 'lower 10%'
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
er_gene.to_csv('R_data/er_gene.csv')
her_gene.to_csv('R_data/her_gene.csv')
tnbc_gene.to_csv('R_data/tnbc_gene.csv')
er_gene.loc[:,'gene'] = er_gene.index
her_gene.loc[:,'gene'] = her_gene.index
tnbc_gene.loc[:,'gene'] = tnbc_gene.index
er_gene.loc[:,'molecular'] = 'ER+'
her_gene.loc[:,'molecular'] = 'HER2+'
tnbc_gene.loc[:,'molecular'] = 'TNBC'
er_gene = pd.concat([er_gene, her_gene, tnbc_gene], axis = 0, ignore_index =  True)
print(er_gene)
er_gene.to_csv('R_data/all_gene.csv')
#'''

#extracting priotorizing genes all breast cancer cell lines
'''sample = pd.read_pickle('data/sample.pkl')
dep = pd.read_pickle('data/ceres_dep.pkl')
sample = sample.reindex( index = dep.index)
sample = sample.loc[ sample.loc[:,'OncotreeLineage'] == 'Breast',:]
eff = pd.read_pickle('data/ceres_dep.pkl')
eff = eff.loc[:,eff.isna().sum() < len(eff.index) * 0.2]
eff = eff.reindex(index = sample.index)
dep = pd.read_pickle('data/ceres_effect.pkl')
dep = dep.reindex(index = sample.index)
dep = dep.reindex(columns = eff.columns)
eff = eff.loc[:,(eff >= 0.5).sum() >= 1]
druggable = pd.read_pickle('data/druggable_genome_update_correct.pkl')
druggable = druggable.loc[druggable.loc[:,'drug'] == 1,:]
eff = eff.loc[:,eff.columns.isin(druggable.index)]
er_gene = pd.DataFrame( np.nan, index = eff.columns, columns = ['perc','median_eff','assignment'])
ceres_control = pd.read_csv('data/23Q4/common_ess_2021.csv', sep = ',',  index_col = 0)
achilles_control = pd.read_csv('data/23Q4/AchillesCommonEssentialControls.csv', sep = ' ', index_col = 0)
er_gene.loc[:,'perc'] = (eff >= 0.5).sum()/len(eff.index)
er_gene.loc[:,'median_eff'] =  dep.reindex(index = eff.index, columns = er_gene.index).median()
er_gene.loc[:,'assignment'] = ''
er_gene.loc[ er_gene.loc[:,'perc'] >= 0.1, 'assignment'] = 'priority'
er_gene.loc[ er_gene.loc[:,'perc'] < 0.1, 'assignment'] = 'lower 10%'
er_gene.loc[ er_gene.index.isin(achilles_control.index) | er_gene.index.isin(ceres_control.index)  , 'assignment'] = 'CommonEssCore'
ceres_control = pd.read_csv('data/23Q4/CRISPRInferredCommonEssentials.csv', sep = ' ',  index_col = 0)
er_gene.loc[  er_gene.index.isin(ceres_control.index)  , 'assignment'] = 'CommonEssCore'
print(er_gene)
print(er_gene.loc[:,'assignment'].value_counts())
er_gene.loc[:,'gene'] = er_gene.index
#er_gene.to_csv('R_data/all_breast_gene_pri.csv')
#'''

#generating priority genes
'''gene = pd.read_csv('R_data/all_gene.csv', index_col=0 )
gene.loc[:,'anno'] = gene.loc[:,'gene']
gene.loc[:,'sel'] = 'above'
gene.loc[ gene.loc[:,'median_eff'] > -0.65, 'anno'] = ''
gene.loc[ gene.loc[:,'median_eff'] > -0.65, 'sel'] = 'below'
gene = gene.loc[ gene.loc[:,'assignment'] == 'priority',:]
print(gene)
#gene.to_csv('R_data/all_gene_PIK3CA.csv')
#'''

#supplementary tables subtype
'''gene = pd.read_csv('R_data/all_gene_PIK3CA.csv', index_col = 0)
print(gene)
gene.loc[:,'median_eff'] = -gene.loc[:,'median_eff']
gene = gene.rename( columns = {'perc':'fraction of cell lines essential','median_eff':'median genetic effect score','sel':'HIG'})
for i in ['ER+','HER2+','TNBC']:
    tmp = gene.loc[gene.loc[:,'molecular'] == i,:]
    tmp.index = tmp.loc[:,'gene']
    tmp.loc[tmp.loc[:,'HIG'] == 'below','HIG'] = ''
    tmp.loc[tmp.loc[:,'HIG'] == 'above','HIG'] = 'yes'
    tmp = tmp.drop(columns = ['gene','anno','assignment','molecular'])
    print(tmp)
    tmp.to_csv('R_data/table_'+i+'.csv')
#'''
#supp tables all breast
'''
gene = pd.read_csv('R_data/all_breast_gene_PIK3CA.csv', index_col = 0)
print(gene)
gene.loc[:,'median_eff'] = -gene.loc[:,'median_eff']
gene = gene.rename( columns = {'perc':'fraction of cell lines essential','median_eff':'median genetic effect score','sel':'HIG'})
tmp = gene
tmp.index = tmp.loc[:,'gene']
tmp.loc[tmp.loc[:,'HIG'] == 'below','HIG'] = ''
tmp.loc[tmp.loc[:,'HIG'] == 'above','HIG'] = 'yes'
tmp = tmp.drop(columns = ['gene','anno','assignment'])
print(tmp)
tmp.to_csv('R_data/table_all_breast.csv')
#'''

'''
gene = pd.read_csv('R_data/all_breast_gene_pri.csv', index_col = 0)
gene.loc[:,'anno'] = gene.loc[:,'gene']
gene.loc[:,'sel'] = 'above'
gene.loc[ gene.loc[:,'median_eff'] > -0.65, 'anno'] = ''
gene.loc[ gene.loc[:,'median_eff'] > -0.65, 'sel'] = 'below'
gene = gene.loc[ gene.loc[:,'assignment'] == 'priority',:]
print(gene)
gene.loc[ gene.loc[:,'gene'] == 'PIK3CA','sel'] = 'exist_drug'
gene.loc[ gene.loc[:,'gene'] == 'TYMS','sel'] = 'exist_drug'
print(gene.loc[gene.loc[:,'sel'] == 'above',:])
gene.to_csv('R_data/all_breast_gene_PIK3CA.csv')
#'''

##########################
#analyzing priority genes#
##########################
#MSK onco info
'''gene = pd.read_csv('R_data/all_gene_PIK3CA.csv', index_col = 0)
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
'''
gene = pd.read_pickle('data/breast_priority_open_target.pkl')
gene.loc[ (gene.loc[:,'sel'] == 'above') & (gene.loc[:,'clinical_phase'] >= 2), 'sel'] = 'exist_drug'
gene.to_csv('R_data/breast_priority_open_target.csv')
#'''

'''query KnownDrugsQuery(
  $ensgId: String!
  $size: Int
  $cursor: String
  $freeTextQuery: String
) {
  target(ensemblId: $ensgId) {
    id
    knownDrugs(cursor: $cursor, freeTextQuery: $freeTextQuery, size: $size) {
      count
      rows {
        phase
        drug {
          id
          name
          mechanismsOfAction {
            rows {
              actionType
            }
          }
        }
      }
    }
  }
}'''
#variables = {"ensgId": "ENSG00000141736","size":2 }

'''
gene = pd.read_pickle('data/breast_priority_gene_oncoinfo.pkl')
ingene = gene.loc[:,'gene'].unique()
ingene = list(ingene)
from gprofiler import GProfiler
gp = GProfiler( return_dataframe = True)
result = gp.convert( organism = 'hsapiens', query = ingene, target_namespace = 'ENSG')
print(result)
#gene.loc[:,'open_target'] = np.nan 
gene.loc[:,'clinical_phase'] = np.nan 

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
          rows { disease { id name } target { id approvedSymbol} drug { id name drugType mechanismsOfAction { rows { mechanismOfAction targets { id approvedSymbol}}}} clinicalPhase clinicalStatus}}}}
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
        r = requests.post(base_url, headers=headers, json={"query": query_string_all, "variables": variables})
        api_response = json.loads(r.text)
        max_phase = 0
        for j in api_response['data']['disease']['chembl']['rows']:
            if j['clinicalPhase'] > max_phase:
                max_phase = j['clinicalPhase']
            if max_phase >= 2:
                break
        gene.loc[gene.loc[:,'gene'] == gene_name,'clinical_phase'] = max_phase
print(gene)
gene.to_pickle('data/breast_priority_open_target.pkl')
# Note: json_data will not be serialized by requests
# exactly as it was in the original request.
#data = '{ "query": "{ genes(first:20000) {nodes {name geneCategoriesWithSources {name} } } }" }'
#response = requests.post('https://dgidb.org/api/graphql', headers=headers, data=data)
#'''

#open target drug search for further prioritized gene
'''
gene = pd.read_csv('R_data/breast_priority_open_target.csv')
gene = gene.loc[gene.loc[:,'sel'] == 'above',:]
ingene = list(gene.loc[:,'gene'].unique())
from gprofiler import GProfiler
gp = GProfiler( return_dataframe = True)
result = gp.convert( organism = 'hsapiens', query = ingene, target_namespace = 'ENSG')
print(result)
#gene.loc[:,'open_target'] = np.nan 

import requests
import json

query_string_cnt =  """
    query KnownDrugsQuery(
  $ensgId: String!
  $cursor: String
  $freeTextQuery: String
  $size: Int = 1
) {
  target(ensemblId: $ensgId) {
    id
    knownDrugs(cursor: $cursor, freeTextQuery: $freeTextQuery, size: $size) {
      count
      
    }
  }
}
   """
query_string_all = """ query KnownDrugsQuery(
  $ensgId: String!
  $cursor: String
  $freeTextQuery: String
  $size: Int
) {
  target(ensemblId: $ensgId) {
    id
    knownDrugs(cursor: $cursor, freeTextQuery: $freeTextQuery, size: $size) {
      count
      cursor
      rows {
        phase
        status
        urls {
          name
          url
        }
        disease {
          id
          name
        }
        drug {
          id
          name
          mechanismsOfAction {
            rows {
              actionType
              targets {
                id
              }
            }
          }
        }
        drugType
        mechanismOfAction
      }
    }
  }
}"""

base_url = "https://api.platform.opentargets.org/api/v4/graphql"
headers = {
    # Already added when you pass json= but not when you pass data=
    # 'Content-Type': 'application/json',
}
output = []
for i in result.index:
    print(i)
    gene_id = result.loc[i,'converted']
    gene_name = result.loc[i,'incoming']
    variables = {"ensgId": gene_id,  "size": 1}
    r = requests.post(base_url, headers=headers, json={"query": query_string_cnt, "variables": variables})
    api_response = json.loads(r.text)
    send_cnt = api_response['data']['target']['knownDrugs']['count']
    if send_cnt > 0:
        curs = 0
        print(gene_id)
        variables = {"ensgId": gene_id, "size": send_cnt}
        r = requests.post(base_url, headers=headers, json={"query": query_string_all, "variables": variables})
        api_response = json.loads(r.text)
        print(api_response)
        index = [tmp_cnt for tmp_cnt in range(send_cnt)]
        tmp = pd.DataFrame(np.nan, index = index, columns = ['gene_id','gene_name','drug_name','drug_id','disease','disease_id','actionType','drugType','moa','phase','status','trial_name','url'])
        for j in api_response['data']['target']['knownDrugs']['rows']:
            flag = 0
            for k in j['drug']['mechanismsOfAction']['rows']:
                if k['actionType'] == 'INHIBITOR':
                    for t in k['targets']:
                        if t['id'] == gene_id:
                            flag = 1
                            break
                    if flag == 1:
                        break
                if flag == 1:
                    break
            if flag == 0:
                continue
            tmp.loc[curs,'gene_id'] = gene_id
            tmp.loc[curs,'gene_name'] = gene_name
            tmp.loc[curs,'drug_name'] = j['drug']['name']
            tmp.loc[curs,'drug_id'] = j['drug']['id']
            tmp.loc[curs,'disease'] = j['disease']['name']
            tmp.loc[curs,'disease_id'] = j['disease']['id']
            tmp.loc[curs,'drugType'] = j['drugType']
            tmp.loc[curs,'moa'] = j['mechanismOfAction']
            tmp.loc[curs,'phase'] = j['phase']
            curs = curs + 1
        tmp = tmp.dropna( how = 'all')
        output.append(tmp)
        print(tmp)
output = pd.concat(output, axis = 0, ignore_index = True)
print(output)
gene.to_pickle('data/breast_priority_open_target_drug.pkl')
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

###################
#mutation analysis#
###################
#CNV translation to MSK
'''cancer_gene = pd.read_csv('data_2023/cancerGeneList.tsv', sep = '\t', index_col = 0)
cnv = pd.read_pickle('data_2023/cnv_original.pkl')
cnv = cnv.loc[:,cnv.columns.isin(cancer_gene.index)]
print(cnv)
cnv = 2**cnv - 1
def tcnv(x):
    y = 0
    if x > 2**0.75:
        y = 2
    elif x > 2**0.4:
        y = 1
    elif x < 2 **(-1.2):
        y = -2
    elif x < 2 **(-0.6):
        y = -1
    return y
cnv = cnv.applymap(tcnv)
print(cnv)
cnv = cnv.T
cnv.index.name = 'Hugo_Symbol'
print(cnv)
cnv.to_csv('OncoKB_2023/cnv_input.txt', sep = '\t')
#'''

#mutational data extraction
'''mut = pd.read_csv('data/23Q4/mut_maf.txt', sep= '\t')
print(mut)
print(mut.columns)
mut.to_pickle('data/mut_maf_oncokb.pkl')
#'''
'''
mut = pd.read_pickle('data/mut_maf_oncokb.pkl')
print(mut.iloc[0,:])
print(mut.loc[:,'Variant_annotation'].value_counts())
#'''

'''
sample = pd.read_csv('data/23Q4/OmicsProfiles.csv', sep =',', index_col = 0)
mut = pd.read_csv('data/23Q4/oncokb_mut', sep = '\t')
print(sample)
mut.loc[:,'DepMap_ID'] = ''
for i in sample.index:
    mut.loc[mut.loc[:,'Tumor_Sample_Barcode'] == i,'DepMap_ID'] = sample.loc[i,'ModelID']
#print(mut.iloc[0,:])
#print(mut.loc[:,'ONCOGENIC'].value_counts())
print(mut)
mut.to_pickle('data/23Q4/oncokb_mut.pkl')
#'''

'''
mut = pd.read_csv('data/23Q4/oncokb_cnv', sep = '\t')
print(mut.iloc[0,:])
print(mut.loc[:,'ONCOGENIC'].value_counts())
mut.to_pickle('data/23Q4/oncokb_cnv.pkl')
#'''

#oncogenic mutation according to subtype
'''mut_type = 'mut'
mut = pd.read_pickle('data/23Q4/oncokb_'+mut_type+'.pkl')
if mut_type == 'mut':
    mut = mut.loc[:, ['Hugo_Symbol','DepMap_ID','ONCOGENIC']]
else:
    mut = mut.loc[:,['HUGO_SYMBOL','SAMPLE_ID','ONCOGENIC']]
mut.columns = ['Hugo_Symbol','DepMap_ID','ONCOGENIC']
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
mut = mut.loc[ mut.loc[:,'DepMap_ID'].isin( sample.index),:]
mut = mut.loc[ mut.loc[:,'ONCOGENIC'].isin( ['Oncogenic','Likely Oncogenic']),:]
mut.loc[:,'ONCOGENIC'] = 1
mut = pd.pivot_table(mut, values = 'ONCOGENIC', index = ['DepMap_ID'], columns = ['Hugo_Symbol'], fill_value = 0)
mut = mut.reindex(index = sample.index)
mut = mut.fillna(0)
eff = pd.read_pickle('data/ceres_effect.pkl')
eff = eff.reindex(index = sample.index)
gene = pd.read_pickle('data/breast_priority_open_target.pkl')
gene.loc[ (gene.loc[:,'sel'] == 'above') & (gene.loc[:,'clinical_phase'] >= 2), 'sel'] = 'exist_drug'
from scipy import stats
for i in ['ER+','HER2+','TNBC']:
    print(i)
    mut_sub = mut.loc[ sample.loc[:,'LegacySubSubtype'] == i,:]
    if len(mut_sub) * 0.1 < 2:
        cut_off = 2
    else:
        cut_off = len(mut_sub)*0.1
    mut_sub = mut_sub.loc[:, mut_sub.sum() >= cut_off]
    print( mut_sub.sum().sort_values())
    print( mut_sub.mean().sort_values())
#'''

#oncogenic mutation all breast
'''mut_type = 'cnv'
mut = pd.read_pickle('data/23Q4/oncokb_'+mut_type+'.pkl')
if mut_type == 'mut':
    mut = mut.loc[:, ['Hugo_Symbol','DepMap_ID','ONCOGENIC']]
else:
    mut = mut.loc[:,['HUGO_SYMBOL','SAMPLE_ID','ONCOGENIC']]
mut.columns = ['Hugo_Symbol','DepMap_ID','ONCOGENIC']
sample = pd.read_pickle('data/sample.pkl')
dep = pd.read_pickle('data/ceres_dep.pkl')
sample = sample.reindex( index = dep.index)
sample = sample.loc[ sample.loc[:,'OncotreeLineage'] == 'Breast',:]
sample = sample.loc[:, ['StrippedCellLineName','LegacyMolecularSubtype','LegacySubSubtype']]
mut = mut.loc[ mut.loc[:,'DepMap_ID'].isin( sample.index),:]
mut = mut.loc[ mut.loc[:,'ONCOGENIC'].isin( ['Oncogenic','Likely Oncogenic']),:]
mut.loc[:,'ONCOGENIC'] = 1
mut = pd.pivot_table(mut, values = 'ONCOGENIC', index = ['DepMap_ID'], columns = ['Hugo_Symbol'], fill_value = 0)
mut = mut.reindex(index = sample.index)
mut = mut.fillna(0)
eff = pd.read_pickle('data/ceres_effect.pkl')
eff = eff.reindex(index = sample.index)
from scipy import stats
mut_sub = mut.copy()
if len(mut_sub) * 0.1 < 2:
    cut_off = 2
else:
    cut_off = len(mut_sub)*0.1
mut_sub = mut_sub.loc[:, mut_sub.sum() >= cut_off]
print( mut_sub.sum().sort_values())
print( mut_sub.mean().sort_values())
#'''

#t test and cohen d for association
#all breast
'''
mut_type = 'cnv'
mut = pd.read_pickle('data/23Q4/oncokb_'+mut_type+'.pkl')
if mut_type == 'mut':
    mut = mut.loc[:, ['Hugo_Symbol','DepMap_ID','ONCOGENIC']]
else:
    mut = mut.loc[:,['HUGO_SYMBOL','SAMPLE_ID','ONCOGENIC']]
mut.columns = ['Hugo_Symbol','DepMap_ID','ONCOGENIC']
sample = pd.read_pickle('data/sample.pkl')
dep = pd.read_pickle('data/ceres_dep.pkl')
sample = sample.reindex( index = dep.index)
sample = sample.loc[ sample.loc[:,'OncotreeLineage'] == 'Breast',:]
mut = mut.loc[ mut.loc[:,'DepMap_ID'].isin( sample.index),:]
mut = mut.loc[ mut.loc[:,'ONCOGENIC'].isin( ['Oncogenic','Likely Oncogenic']),:]
mut.loc[:,'ONCOGENIC'] = 1
mut = pd.pivot_table(mut, values = 'ONCOGENIC', index = ['DepMap_ID'], columns = ['Hugo_Symbol'], fill_value = 0)
mut = mut.reindex(index = sample.index)
mut = mut.fillna(0)
eff = pd.read_pickle('data/ceres_effect.pkl')
eff = eff.reindex(index = sample.index)
gene = pd.read_csv('R_data/all_breast_gene_pri.csv', index_col = 0)
from scipy import stats
dep = eff.copy()
dep = dep.loc[:, dep.columns.isin( gene.loc[ (gene.loc[:,'assignment'] == 'priority'),'gene'])]
mut_sub = mut.copy()
if len(mut_sub) * 0.1 < 2:
    cut_off = 2
else:
    cut_off = len(mut_sub)*0.1
mut_sub = mut_sub.loc[:, mut_sub.sum() >= cut_off]
#dep = dep.reindex(columns = ['TEAD4','RETSAT'])
p_val = pd.DataFrame( np.nan, index = dep.columns, columns = mut_sub.columns)
cohen_d = pd.DataFrame( np.nan, index = dep.columns, columns = mut_sub.columns)
for j in p_val.index:
    for k in p_val.columns:
        a = dep.loc[ mut_sub.loc[:,k] == 1,j]
        b = dep.loc[ mut_sub.loc[:,k] == 0,j]
        t, p = stats.ttest_ind( a, b)
        p_val.loc[j,k] = p
        s = np.sqrt(((len(a) - 1)*(a.std())**2 + (len(b)-1)*(b.std())**2) / (len(a)+len(b)-2))
        d = (a.mean() - b.mean())/s
        cohen_d.loc[j,k] = d
print(p_val)
print(cohen_d)
i = 'breast'
p_val.to_pickle('data/'+mut_type+'_'+i+'_ttest_p_val.pkl')
cohen_d.to_pickle('data/'+mut_type+'_'+i+'_ttest_cohen_d.pkl')
#'''

#according to subtype
'''mut = pd.read_pickle('data/23Q4/oncokb_cnv.pkl')
mut = mut.loc[:,['HUGO_SYMBOL','SAMPLE_ID','ONCOGENIC']]
mut.columns = ['Hugo_Symbol','DepMap_ID','ONCOGENIC']
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
mut = mut.loc[ mut.loc[:,'DepMap_ID'].isin( sample.index),:]
mut = mut.loc[ mut.loc[:,'ONCOGENIC'].isin( ['Oncogenic','Likely Oncogenic']),:]
mut.loc[:,'ONCOGENIC'] = 1
mut = pd.pivot_table(mut, values = 'ONCOGENIC', index = ['DepMap_ID'], columns = ['Hugo_Symbol'], fill_value = 0)
mut = mut.reindex(index = sample.index)
mut = mut.fillna(0)
eff = pd.read_pickle('data/ceres_effect.pkl')
eff = eff.reindex(index = sample.index)
gene = pd.read_pickle('data/breast_priority_open_target.pkl')
gene.loc[ (gene.loc[:,'sel'] == 'above') & (gene.loc[:,'clinical_phase'] >= 2), 'sel'] = 'exist_drug'
from scipy import stats
for i in ['ER+','HER2+','TNBC']:
    print(i)
    dep = eff.loc[ sample.loc[:,'LegacySubSubtype'] == i,:]
    #dep = dep.loc[:, dep.columns.isin( gene.loc[ (gene.loc[:,'sel'] == 'above') & (gene.loc[:,'molecular'] == i),'gene'])]
    dep = dep.loc[:, dep.columns.isin( gene.loc[ (gene.loc[:,'molecular'] == i),'gene'])]
    mut_sub = mut.loc[ sample.loc[:,'LegacySubSubtype'] == i,:]
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
            t, p = stats.ttest_ind( a, b)
            p_val.loc[j,k] = p
            s = np.sqrt(((len(a) - 1)*(a.std())**2 + (len(b)-1)*(b.std())**2) / (len(a)+len(b)-2))
            d = (a.mean() - b.mean())/s
            cohen_d.loc[j,k] = d
    p_val.to_pickle('data/cnv_'+i+'_ttest_p_val.pkl')
    print(p_val)
    cohen_d.to_pickle('data/cnv_'+i+'_ttest_cohen_d.pkl')
    print(cohen_d)
#'''

'''
output = []
for j in ['mut','cnv']:
    for i in ['ER+','HER2+','TNBC','breast']:
        print(i)
        p_val = pd.read_pickle('data/'+j+'_'+i+'_ttest_p_val.pkl')
        cohen_d = pd.read_pickle('data/'+j+'_'+i+'_ttest_cohen_d.pkl')
        p_val = p_val.stack(dropna = False)
        p_val = pd.DataFrame( p_val.values, index = p_val.index, columns = ['p_val'])
        p_val = p_val.reset_index(names = ['gene_1','gene_2'])
        cohen_d = cohen_d.stack(dropna = False)
        cohen_d = pd.DataFrame( cohen_d.values, index = cohen_d.index, columns = ['cohen_d'])
        p_val.loc[:,'cohen_d'] = cohen_d.values
        p_val.loc[:,'g2_status'] = j
        p_val.loc[:,'molecular'] = i
        output.append(p_val)
        #cohen_d = cohen_d.reset_index(names = ['gene_1','gene_2'])
        print(p_val)
        #print(cohen_d)
result = pd.concat(output, axis = 0)
print(result)
result.to_pickle('data/breast_priority_association_result.pkl')
#'''

#FDR correction
'''result_all = pd.read_pickle('data/breast_priority_association_result.pkl')
import statsmodels.stats.multitest as stat
result_all = result_all.dropna(how = 'any')
#result_all = result_all.loc[ result_all.loc[:,'molecular'] == 'breast',:]
result = result.loc[ result.loc[:,'molecular'] != 'breast',:]
output = []
for i in ['mut','cnv']:
    result = result_all.loc[ result_all.loc[:,'g2_status'] == i,:]
    tmp = result.loc[:,'p_val']
    result.loc[:,'fdr'] = np.nan
    reject, correct, _ , _ = stat.multipletests(tmp, method = 'fdr_bh')
    result.loc[:,'fdr'] = correct
    result.loc[:,'sel'] = 'not_sig'
    result.loc[ (result.loc[:,'fdr'] < 0.1) & (result.loc[:,'cohen_d'].abs() >= 1), 'sel'] = 'sig'
    result.loc[ :,'name'] = ''
    idx = result.loc[:,'sel'] == 'sig'
    result.loc[ idx,'name'] = result.loc[idx,'gene_1'] + '_' + result.loc[idx,'gene_2']
    print(result.loc[ result.loc[:,'sel'] == 'sig',:])
    result.loc[:,'fdr'] = -np.log10( result.loc[:,'fdr'])
    output.append(result)
output = pd.concat( output, axis = 0)
print(output)
#output.to_csv('R_data/mut_association_breast.csv')
result.to_pickle('data/breast_priority_association_fdr_mut.pkl')
#'''


'''result_all = pd.read_pickle('data/breast_priority_association_result.pkl')
import statsmodels.stats.multitest as stat
result_all = result_all.dropna(how = 'any')
result_all = result_all.loc[ result_all.loc[:,'molecular'] != 'breast',:]
#result = result.loc[ result.loc[:,'molecular'] != 'breast',:]
output = []
for i in ['mut','cnv']:
    result = result_all.loc[ result_all.loc[:,'g2_status'] == i,:]
    tmp = result.loc[:,'p_val']
    result.loc[:,'fdr'] = np.nan
    reject, correct, _ , _ = stat.multipletests(tmp, method = 'fdr_bh')
    result.loc[:,'fdr'] = correct
    result.loc[:,'sel'] = 'not_sig'
    result.loc[ (result.loc[:,'fdr'] < 0.1) & (result.loc[:,'cohen_d'].abs() >= 1), 'sel'] = 'sig'
    result.loc[ :,'name'] = ''
    idx = result.loc[:,'sel'] == 'sig'
    result.loc[ idx,'name'] = result.loc[idx,'gene_1'] + '_' + result.loc[idx,'gene_2']
    print(result.loc[ result.loc[:,'sel'] == 'sig',:])
    result.loc[:,'fdr'] = -np.log10( result.loc[:,'fdr'])
    output.append(result)
output = pd.concat( output, axis = 0)
print(output)
output.to_csv('R_data/mut_association_breast_type_specific.csv')
#result.to_pickle('data/breast_priority_association_fdr_mut.pkl')
#'''



'''result = pd.read_pickle('data/breast_priority_association_result.pkl')
import statsmodels.stats.multitest as stat
result = result.dropna(how = 'any')
result = result.loc[ result.loc[:,'g2_status'] == 'cnv',:]
tmp = result.loc[:,'p_val'].dropna()
result.loc[:,'fdr'] = np.nan
reject, correct, _ , _ = stat.multipletests(tmp, method = 'fdr_bh')
result.loc[:,'fdr'] = correct
result.loc[:,'sel'] = 0
result.loc[ (result.loc[:,'fdr'] < 0.1) & (result.loc[:,'cohen_d'].abs() >= 1), 'sel'] = 1
result.loc[ :,'name'] = ''
idx = result.loc[:,'sel'] == 1
result.loc[ idx,'name'] = result.loc[idx,'gene_1'] + '_' + result.loc[idx,'gene_2']
print(result.loc[:,'sel'].sum())
print(result.loc[ result.loc[:,'sel'] == 1,:])
#result.to_pickle('data/breast_priority_association_fdr_mut.pkl')
#'''


#heatmap annotation
'''gene = pd.read_csv('R_data/breast_priority_open_target.csv')
gene = gene.loc[gene.loc[:,'sel'] == 'above',:]
gene.to_csv('R_data/priority_list.csv')
dep = pd.read_pickle('data/ceres_effect.pkl')
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
sample.loc['ACH-001396','LegacySubSubtype'] = 'ER+'
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
#dep.to_csv('R_data/priority_dep_breast_cell_line.csv')
#dep.to_csv('R_data/priority_dep_breast_cell_line.csv')
#'''
