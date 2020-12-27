'''Calculate the Concordance Ratio, Discordance Ration and Orthogonality Score for each small molecule compared with a reference compound and the disease'''
import pandas as pd
import os,glob

#cutoff: int or float type, the TCS score of the reference compound needs to be greater than the cutoff
#inpath_innerids: inner_LT_final.txt
#infolder_lauchedtcs07: the output of lauchedTCSinnerGene() function in calcTCS.py, i.e. outfolder_lauchedtcs07
#infolder_reftcs: the output of RobustTCS() function in calcTCS.py, i.e. outfolder_reftcs
#detransLINCS_path: the output of transID() function in disease.py, i.e. detransLINCS_path

def S_CRscore(cutoff,infolder_lauchedtcs07,infolder_reftcs,outfolder_score):
    for ref_inpath in glob.glob(infolder_reftcs+'/BRD*'):
        ref_fname = ref_inpath.strip().split('/')[-1]
        ref_pertid = ref_fname.split('_')[0]
        ref_tcsnum = ref_fname.split('_')[1].split('.')[0]
        
        if float(ref_tcsnum) >= cutoff:
            score_outfolder = outfolder_score+'/'+ref_pertid
            if os.path.exists(score_outfolder) == True:
                print('the path already exists: {}'.format(score_outfolder))
            else:
                os.makedirs(score_outfolder)
                ref_indf = pd.read_csv(ref_inpath,index_col=0)
                ref_genes = ref_indf.index.tolist()
                ref_dic = {gene_id:ref_indf.at[gene_id,'val'] for gene_id in ref_genes}
                with open(score_outfolder+'/'+ref_pertid+'_scr.txt','w') as scr:
                    scr.write('pert_id\tS\tCR\n')
                    for drug_inpath in glob.glob(infolder_lauchedtcs07+'/BRD*'):
                        outlist = []
                        drug_pertid = drug_inpath.strip().split('/')[-1].split('_')[0]
                        drug_indf = pd.read_csv(drug_inpath,index_col=0)
                        drug_genes = drug_indf.index.tolist()
                        inner_genes = set(drug_genes)&(set(ref_genes))
                        s_score = float(len(inner_genes))/len(drug_genes)
                        same,oppo = 0,0
                        drug_dic = {j:drug_indf.at[j,'val'] for j in drug_genes}
                        for k,v in ref_dic.items():
                            if k in drug_genes:
                                if v == drug_dic[k]:
                                    same += 1
                                else:
                                    oppo += 1
                        cr_score = (same+0.00001)/(oppo+0.00001)
                        outlist.append(str(drug_pertid))
                        outlist.append(str(s_score))
                        outlist.append(str(cr_score))
                        outline = '\t'.join(outlist)+'\n'
                        scr.write(outline)

def DRscore(inpath_innerids,detransLINCS_path,infolder_lauchedtcs07,infolder_reftcs,outfolder_score):
    innerid_df = pd.read_csv(inpath_innerids,sep='\t')
    innerid_list = innerid_df['LINCS_ID'].tolist()
    disease_df = pd.read_csv(detransLINCS_path)
    disease_genes = disease_df['id'].tolist()
    disease_dict = {disease_df.iat[i,0]:disease_df.iat[i,1] for i in disease_df.index}
    
    refid_list = os.listdir(outfolder_score)
    for ref_inpath in glob.glob(infolder_reftcs+'/BRD*'):
        ref_pertid = ref_inpath.strip().split('/')[-1].split('_')[0]
        if ref_pertid in refid_list:
            ref_indf = pd.read_csv(ref_inpath)
            ref_genes = ref_indf['id'].tolist()
            ref_dict = {ref_indf.iat[i,0]:ref_indf.iat[i,1] for i in ref_indf.index}
            other_genes = (set(innerid_list)-set(ref_genes))&set(disease_genes)
            
            with open(outfolder_score+'/'+ref_pertid+'/'+ref_pertid+'_dr.txt','w') as drf:
                drf.write('pert_id\tDR\n')
                for drug_inpath in glob.glob(infolder_lauchedtcs07+'/BRD*'):
                    outlist = []
                    drug_pertid = drug_inpath.strip().split('/')[-1].split('_')[0]
                    drug_indf = pd.read_csv(drug_inpath)
                    inner_drug_df = drug_indf[drug_indf['id'].isin(other_genes)]
                    inner_drug_df.index = range(len(inner_drug_df))
                    drug_dict = {inner_drug_df.iat[i,0]:inner_drug_df.iat[i,1] for i in inner_drug_df.index}
                    same,oppo = 0.00001,0.00001
                    for k,v in drug_dict.items():
                        if v == disease_dict[k]:
                            same += 1
                        else:
                            oppo += 1
                    dr_score = oppo/same
                    outlist.append(drug_pertid)
                    outlist.append(str(dr_score))
                    outline = '\t'.join(outlist)+'\n'
                    drf.write(outline)

def ConcatScore(outfolder_score):
    for ref_pertid in os.listdir(outfolder_score):
        scr_path = glob.glob(outfolder_score+'/'+ref_pertid+'/*_scr.txt')[0]
        dr_path = glob.glob(outfolder_score+'/'+ref_pertid+'/*_dr.txt')[0]
        scr_df = pd.read_csv(scr_path,sep='\t',index_col=0)
        dr_df = pd.read_csv(dr_path,sep='\t',index_col=0)
        concat_df = pd.concat([scr_df,dr_df],axis=1)
        concat_df.to_csv(outfolder_score+'/'+ref_pertid+'/'+ref_pertid+'_score.csv')

def NormalizeOS(outfolder_score):
    import os,glob,math
    import pandas as pd
    
    for ref_pertid in os.listdir(outfolder_score):
        score_path = glob.glob(outfolder_score+'/'+ref_pertid+'/*_score.csv')[0]
        score_df = pd.read_csv(score_path,index_col=0)
        pert_ids = score_df.index.tolist()
        cr_max = float(score_df['CR'].max())
        dr_max = float(score_df['DR'].max())
        with open(outfolder_score+'/'+ref_pertid+'/'+ref_pertid+'_normalize.txt','w') as outf:
            outf.write('small_molecule\tS\tCR\tDR\tOS\n')
            for pert_id in score_df.index:
                outlist = []
                s = score_df.at[pert_id,'S']
                cr = float(score_df.at[pert_id,'CR'])/cr_max
                dr = float(score_df.at[pert_id,'DR'])/dr_max
                os = math.sqrt((1-cr)**2+dr**2)
                outlist.append(pert_id)
                outlist.append(str(s))
                outlist.append(str(cr))
                outlist.append(str(dr))
                outlist.append(str(os))
                outline = '\t'.join(outlist)+'\n'
                outf.write(outline)