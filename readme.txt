#Description
disease.py: screen differentially expressed genes, and then obtain disease-signature.
calcTCS.py: parse gct and gctx files (Enache, et al. 2017), and then calculate transcriptional consensus signature (TCS) for each compound using the plate-normalized Level 4 L1000 data
robustTCS.py: obtain a  robust reference-signature for each reference compound
score.py: calculate concordance ratio (CR), disease-specific discordance ratio (DR) and Orthogonality Score (OS) for drug-drug-disease combination, as defined in the study of Stathias et al. (2018).
network_proximity.py: calculate the network-based proximity of drug-drug-disease combinations (Cheng, et al. 2019).
TissueProteins.py: screen tissue-specific proteins from the gene expression data of TCGA and GTEx, and generate a sif file to create tissue-specific network

# Input arguments
inpath_inst: GSE70138_Broad_LINCS_inst_info.txt
inpath_inst: GSE70138_Broad_LINCS_inst_info.txt
level4_path: GSE92742_Broad_LINCS_Level4_ZSPCINF_mlr12k_n1319138x12328.gctx
inpath_innerid: inner_LT_final.txt
inpath_lauched: lauched_Repurposing.csv
GSE70138_Broad_LINCS_inst_info.txt, GSE70138_Broad_LINCS_inst_info.txt and GSE92742_Broad_LINCS_Level4_ZSPCINF_mlr12k_n1319138x12328.gctx are downloaded from GEO database with the accession number GSE92742.
lauched_Repurposing.csv is pre-processed after being downloaded from https://clue.io/repurposing-app.
orig1_path: a folder of RNA-Seq files downloaded from TCGA
json_path: the json file downloaded from TCGA
inpath_edge: a sif file of PPI network, PPI.sif, downloaded from the work of Cheng et al. (2019)
inpath_dict: a directory like this: {'disease':{'prostate cancer':[disease-associated genes]},'ref_drug':{'reference':[tar1,tar2,...]},'cand_drug':{'candidate 1':[]}}
inpath_tumor_tpm, inpath_normal_tpm: tumorallgenesmedian_Prostate.csv, normalallgenesmedian_Prostate.csv, files preprocessed after downloading from the UCSC Xena Portal (Vivian J, et al. 2017)
inpath_idmap_ppi: nodes_entrez_ensembl.csv, containing the annotation of genes in the entire PPI, with two columns——ENTREZID, ENSEMBL
inpath_entireppi: entire_PPI.csv which is derived from the supporting information of Cheng et al. (2019)

#References
Enache, Oana M., Lahr, David L., and Natoli, Ted E., et al. (2017). The GCTx format and cmap{Py, R, M} packages: resources for the optimized storage and integrated traversal of dense matrices of data and annotations. bioRxiv, doi:10.1101/227041.
Stathias, V., Jermakowicz, A.M., and Maloof, M.E., et al. (2018). Drug and disease signature integration identifies synergistic combinations in glioblastoma. Nature Communications 9, 5315. 
Cheng, F., Kovács, I.A., and Barabási, A.-L. (2019). Network-based prediction of drug combinations. Nature Communications 10, 1197.
Vivian J, Rao AA, Nothaft FA, et al. (2017). Toil enables reproducible, open source, big biomedical data analyses. Nature Biotechnology. 35, 314-316.