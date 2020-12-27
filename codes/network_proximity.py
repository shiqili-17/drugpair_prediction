import pandas as pd
import wrappers

#inpath_edge: a sif file of PPI network
def construct_network(inpath_edge):
    network = wrappers.get_network(inpath_edge, only_lcc = True)
    return network

#inpath_dict: a directory like this: {'disease':{'prostate cancer':[disease-associated genes]},'ref_drug':{'reference':[tar1,tar2,...]},'cand_drug':{'candidate 1':[]}}
#outpath_zscore: a csv file with three columns——'disease','drug','z-score'
def DiseaseDrugScore(input_network,inpath_dict,outpath_zscore):
    outdata_zscore = []
    inpath_dict = eval(open(inpath_dict).read())
    for k,v in inpath_dict['disease'].items():
        disease_name = k
        disease_genes = [str(i) for i in v]

    drug_dict = {}
    for ref_drug,ref_genes in inpath_dict['ref_drug'].items():
        drug_dict.setdefault(ref_drug,ref_genes)
    for cand_drug,cand_genes in inpath_dict['cand_drug'].items():
        drug_dict.setdefault(cand_drug,cand_genes)
    for drug_name,drug_genes_int in drug_dict.items():
        drug_genes = [str(i) for i in drug_genes_int]
        d,z,(mean,sd) = wrappers.calculate_proximity(input_network,drug_genes,disease_genes,min_bin_size=100)
        outdata_zscore.append([disease_name,drug_name,z])
    outdf_zscore = pd.DataFrame(data=outdata_zscore,columns=['disease','drug','z-score'])
    outdf_zscore.to_csv(outpath_zscore,index=False)

#outpath_sscore: a csv file with three columns——'ref_drug','cand_drug','s_score'
def DrugDrugScore(input_network,inpath_dict,outpath_sscore):
    outdata_sscore = []
    inpath_dict = eval(open(inpath_dict).read())
    for ref_drug,ref_genes in inpath_dict['ref_drug'].items():
        for cand_drug,cand_genes in inpath_dict['cand_drug'].items():

            nodes_network = set(input_network.nodes())
            nodes_from = set(ref_genes) & nodes_network 
            nodes_to = set(cand_genes) & nodes_network
            d = wrappers.get_separation(input_network, nodes_from, nodes_to)
            outdata_sscore.append([ref_drug,cand_drug,d])
    outdf_sscore = pd.DataFrame(data=outdata_sscore,columns=['ref_drug','cand_drug','s_score'])
    outdf_sscore.to_csv(outpath_sscore,index=False)