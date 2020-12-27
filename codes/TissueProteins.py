'''Screen tissue-specific proteins from the gene expression data of TCGA and GTEx'''
import pandas as pd
#cutoff: int or float type. Genes whose median expression value was greater than the cut-off value were screened out
#inpath_tumor_tpm: tumorallgenesmedian_Prostate.csv
#inpath_normal_tpm: normalallgenesmedian_Prostate.csv
#inpath_idmap_ppi: nodes_entrez_ensembl.csv, containing the annotation of genes in the entire PPI, with two columns——ENTREZID,ENSEMBL
#prcaexpr: prcaexpr_1tpm_del.csv a csv file derived from the out file of ScreenGenesTPM() function, i.e. outpath_tpm, with three columns——Entrezid (without version number), tumor_TPM_median, normal_TPM_median
#inpath_prcainppi: the output of PrcaInPPI() function, i.e. PrcaInPPI, a list of prostate cancer-associated proteins
#inpath_entireppi: entire_PPI.csv
#outpath_prcasif: a sif file with three columns, used as the input to construct a network


def ScreenGenesTPM(cutoff,inpath_tumor_tpm,inpath_normal_tpm,outpath_tpm):
    indf_tumor = pd.read_csv(inpath_tumor_tpm)
    indf_normal = pd.read_csv(inpath_normal_tpm)
    tumor_genes = indf_tumor['Gene_id'][indf_tumor['TPM_median']>cutoff].tolist()
    normal_genes = indf_normal['Gene_id'][indf_normal['TPM_median']>cutoff].tolist()
    outdf_tumor = indf_tumor[indf_tumor['Gene_id'].isin(tumor_genes+normal_genes)]
    outdf_normal = indf_normal[indf_normal['Gene_id'].isin(tumor_genes+normal_genes)]
    outdf_tumor = outdf_tumor.rename(columns = {'TPM_median':'tumor_TPM_median'})
    outdf_normal = outdf_normal.rename(columns = {'TPM_median':'normal_TPM_median'})
    outdf_expr = pd.concat([outdf_tumor,outdf_normal],axis = 1)
    outdf_expr.to_csv(outpath_tpm,index = False)

'''Count how many of the genes expressed in prostate tissue are in the PPI network
   Take the intersection of genes expressed in prostate tissue and also included in the entire PPI'''
def PrcaInPPI(inpath_idmap_ppi,prcaexpr,outpath_prcainppi):
    indf_ppi = pd.read_csv(inpath_idmap_ppi)
    emsembl_ppi = indf_ppi['ENSEMBL'].tolist()
    print('len(ppi_genes): {}'.format(len(emsembl_ppi)))
    indf_prca = pd.read_csv(prcaexpr)
    emsembl_prca = indf_prca['Gene_id'].tolist()
    print('len(prca_genes): {}'.format(len(emsembl_prca)))

    inner_ensembl = [i for i in set(emsembl_ppi)&set(emsembl_prca)]
    print('len(inner_genes): {}'.format(len(inner_ensembl)))
    
    inner_entrezid = indf_ppi['ENTREZID'][indf_ppi['ENSEMBL'].isin(inner_ensembl)].tolist()
    print('len(inner_entrezid): {}'.format(len(inner_entrezid)))
    open(outpath_prcainppi,'w').write(str(inner_entrezid))
    
'''Generate a sif file according to the proteins expressed in the specific tissue'''
def PrcaEdge(inpath_prcainppi,inpath_entireppi,outpath_prcappi,outpath_prcasif):
    inlist_prcainppi = eval(open(inpath_prcainppi).read())
    indf_entireppi = pd.read_csv(inpath_entireppi)
    outdata = []
    outcolumn = ['Protein_A_Entrez_ID','Protein_B_Entrez_ID']
    with open(outpath_prcasif,'w') as outsif:
        for ix in indf_entireppi.index:
            node_a = indf_entireppi.iat[ix,0]
            node_b = indf_entireppi.iat[ix,1]
            if node_a in inlist_prcainppi or node_b in inlist_prcainppi:
                outdata.append([node_a,node_b])
                outsif.write('\t'.join([str(node_a),'1',str(node_b)])+'\n')
            else:
                continue
        print('edges in prcappi: {}'.format(len(outdata)))
        outdf_prcappi = pd.DataFrame(data = outdata,columns = outcolumn)
