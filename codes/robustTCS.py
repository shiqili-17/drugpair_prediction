import pandas as pd
import os, glob, shutil

#inpath_lauched: lauched_Repurposing.csv
#inpath_innerids: inner_LT_final.txt
#inpath_inst: GSE70138_Broad_LINCS_inst_info.txt
#infolder_txt02: the output of parseGct(lincs,level4_path) in calcTCS.py, i.e. txt_02
#infolder_trans03: the output of parseGct() function in calcTCS.py, i.e. trans_03

'''id in the Repurposing_Hub: BRD-K17443395-065-03-0,its corresponding pert-id in the LINCS is BRD-K17443395'''
def GetId(longid):
    short_ids = longid.strip().split('-')
    out_id = short_ids[0]+'-'+short_ids[1]
    return out_id

#outpath_candids: the file path of a list containing pert_ids of prostate cancer drugs
def CandRefIds(infolder_txt02,inpath_lauched,outpath_candids):
    all_trtids = os.listdir(infolder_txt02)
    lauched_indf = pd.read_csv(inpath_lauched,index_col=0)
    innerpro_ids = []
    for name in lauched_indf.index:
        indication_list = [i.strip() for i in str(lauched_indf.at[name,'Indication']).strip().split(',')]
        if 'prostate cancer' in indication_list:
            longid_list1 = [i.strip() for i in lauched_indf.at[name,'Id'].strip().split(',')]
            pertid_list1 = [GetId(longid=j) for j in longid_list1]
            deprecated_ids =  lauched_indf.at[name,'Deprecated ID']
            if str(deprecated_ids) == 'nan':
                pertid_list2 = []
            else:
                longid_list2 = [i.strip() for i in deprecated_ids.strip().split(',')]
                pertid_list2 = [GetId(longid=j) for j in longid_list2]
            id_list = pertid_list1 + pertid_list2
            for each_id in id_list:
                if each_id not in innerpro_ids and each_id in all_trtids:
                    innerpro_ids.append(each_id)
                else:
                    pass
    print('the number of all_prostate_ids:{}'.format(len(innerpro_ids)))
    open(outpath_candids+'/innerpro_ids','w').write(str(innerpro_ids))

#outfolder_copy: copy files of data treating PC-3 cells from infolder_trans03
def CopyRefPC3(inpath_idlist,infolder_trans03,outfolder_copy):
    for each_id in eval(open(inpath_idlist).read()):
        old_path = infolder_trans03+'/'+each_id+'/PC3.csv'
        new_path = outfolder_copy+'/'+each_id+'.csv'
        shutil.copyfile(old_path,new_path)

'''Calculate respective Robust TCS for each reference compound,
   genes that were concordantly over/underexpressed in at least half of the treated prostate cancer cell lines'''
def RobustTCS(inpath_innerids,infolder_trans,outfolder_reftcs):
    innerid_df = pd.read_csv(inpath_innerids,sep='\t')
    innerid_list = innerid_df['LINCS_ID'].tolist()
    
    for trans_path in glob.glob(infolder_trans+'/BRD*'):
        pert_id = trans_path.strip().split('/')[-1].split('.')[0]
        pc3_indf = pd.read_csv(trans_path,index_col=0)
        cutoff = 0.5*len(pc3_indf.columns.tolist()) #50%samples
        ids_list, values_list = [], []
        up = ((pc3_indf==1).sum(axis=1)).tolist()
        dn = ((pc3_indf==-1).sum(axis=1)).tolist()
        genes = pc3_indf.index.tolist()
        for ix,g in enumerate(genes):
            if up[ix]>cutoff:
                ids_list.append(g)
                values_list.append(1)
            elif dn[ix]>cutoff:
                ids_list.append(g)
                values_list.append(-1)
            else:
                continue
        data = {'id':ids_list,'val':values_list}
        ref_outdf = pd.DataFrame(data=data)
        inner_outdf = ref_outdf[ref_outdf['id'].isin(innerid_list)]
        tcs_num = len(inner_outdf.index.tolist())
        outfname = pert_id+'_'+str(tcs_num)+'.csv'
        inner_outdf.to_csv(outfolder_reftcs+'/'+outfname,index=False)

#outpath_tcsnum: a csv file with three columns: pert_id, pert_iname, tcs_num
def RobustInnerNum(inpath_inst,infolder_reftcs,outpath_tcsnum):
    inst_indf = pd.read_csv(inpath_inst,sep = '\t')
    id_list,name_list,num_list = [],[],[]
    for reftcs_path in glob.glob(infolder_reftcs+'/*.csv'):
        each_fname = reftcs_path.strip().split('/')[-1]
        pert_id = each_fname.split('_')[0]
        tcs_num = each_fname.split('_')[1].split('.')[0]
        pert_iname = [i for i in set(inst_indf['pert_iname'][inst_indf['pert_id']==pert_id])][0]
        
        id_list.append(pert_id)
        name_list.append(pert_iname)
        num_list.append(tcs_num)
    outdf = pd.DataFrame(data = {'pert_iname':name_list,'robust_tcs':num_list},
                        index = pd.Index(id_list,name='pert_id'))
    outdf.to_csv(outpath_tcsnum)