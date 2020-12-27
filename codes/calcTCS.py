import pandas as pd
import os, glob, shutil
from cmapPy.pandasGEXpress.parse import parse
import cmapPy.pandasGEXpress.concat as cg
import cmapPy.pandasGEXpress.write_gct as wg
from collections import Counter

#inpath_inst: GSE70138_Broad_LINCS_inst_info.txt
#level4_path: GSE92742_Broad_LINCS_Level4_ZSPCINF_mlr12k_n1319138x12328.gctx
#inpath_innerid: inner_LT_final.txt
#inpath_lauched: lauched_Repurposing.csv

'''Screen compounds meeting the criteria: 'pert_type'='trt_cp'&'pert_time'=24'''
#all_pertid： the list of the screened compounds
#all_instdf: the dataframe of the inst_info of the screened compounds
def getSmcp(inpath_inst):
    meta_inst = pd.read_csv(inpath_inst,sep = '\t')
    all_instdf = meta_inst[(meta_inst['pert_type']=='trt_cp')&(meta_inst['pert_time']==24)]
    all_pertid = set(all_instdf['pert_id'])
    print("the number of compounds：%d"%(len(all_pertid)))
    return all_instdf,all_pertid

'''Parse the gctx file: GSE92742_Broad_LINCS_Level4_ZSPCINF_mlr12k_n1319138x12328.gctx'''
#lincs = GSE70138
# gct_01: a folder containg folders for each compound with the name of its pert_id, containing the GCT file of each cell sample treated by the compound
# txt_02: text files derived from the .gct, with the sample as the column name
# trans_03:csv files derived from the .txt, converting the value to 1, if its z-score>=1; to -1, if z-score<=-1; otherwise, to 0'''
def trans(x):
    if x <= -1:
        return -1
    elif x >= 1:
        return 1
    else:
        return 0
def parseGct(lincs,level4_path):
    for _pertid in all_pertid:
        os.mkdir(lincs+'/gct_01/'+_pertid)
        os.mkdir(lincs+'/txt_02/'+_pertid)
        os.mkdir(lincs+'/trans_03/'+_pertid)
        gct_path = lincs+'/gct_01/'+_pertid
        txt_path = lincs+'/txt_02/'+_pertid
        trans_path = lincs+'/trans_03/'+_pertid
        _instdf1 = all_instdf[all_instdf['pert_id']==_pertid]
        _cellids = set(_instdf1['cell_id'])
        for _cellid in _cellids:
            _instdf2 = _instdf1[_instdf1['cell_id']==_cellid]
            _instids = set(_instdf2['inst_id'])
            gct_list = []
            for _instid in _instids:
                per_gctoo = parse(level4_path,cid = _instid)
                per_info = _instdf2[_instdf2['inst_id']==_instid]
                per_info.set_index('inst_id',inplace = True)
                per_gctoo.col_metadata_df = per_info
                gct_list.append(per_gctoo)
            gctoos = cg.hstack(gct_list)
            gctpath = gct_path+'/'+_cellid+'.gct'
            wg.write(gctoos,gctpath)
            
            with open(txt_path+'/'+_cellid+'.txt','w') as out_txt:
                for ix,line in enumerate(open(gctpath).readlines()):
                    if lincs == 'GSE92742_1':
                        if ix == 2 or ix >12:
                            out_txt.write(line)
                        else:
                            continue
                    elif lincs == 'GSE70138_2':
                        if ix == 2 or ix >13:
                            out_txt.write(line)
                        else:
                            continue
            
            indf_txt = pd.read_csv(txt_path+'/'+_cellid+'.txt',sep='\t')
            col_list = indf_txt.columns.tolist()[1:]
            for col in col_list:
                indf_txt[col] = indf_txt[col].apply(lambda x: trans(x))
            indf_txt.to_csv(trans_path+'/'+_cellid+'.csv',index = False)

'''Determine whether the number of cell lines treated by each compound >2
   Generate TCS for each compound treating multiple cell lines'''
# tcs_04: csv files derived from files in trans_03, containing only genes in the TCS.
def TCSGene(infolder_trans03,outfolder_tcs04):
    for tans3_path in glob.glob(infolder_trans03):
        _pertid = tans3_path.split('/')[-1]
        csv_ps = glob.glob(tans3_path+'/*.csv')
        if len(csv_ps) > 2:
            df3_list,agg,tcs = [],[],[]
            cutoff_cell = len(csv_ps)*0.3
            for csv_p in csv_ps:
                indf = pd.read_csv(csv_p,index_col=0)
                df3_list.append(indf)
                genes = indf.index.tolist()
                cutoff_sample = len(indf.columns.tolist())*0.2
                up = ((indf==1).sum(axis=1)).tolist()
                dn = ((indf==-1).sum(axis=1)).tolist()
                for ix,g in enumerate(genes):
                    if up[ix]>cutoff_sample or dn[ix]>cutoff_sample:
                        agg.append(g)
                    else:
                        continue
            count_dict = dict(Counter(agg))
            for k,v in count_dict.items():
                if v > cutoff_cell:
                    tcs.append(k)
                else:
                    continue

            try:
                cat_df = pd.concat(df3_list,axis=1)
            except:
                continue
            tcs_df = cat_df.loc[tcs,:]
            out_p = outfolder_tcs04+'/'+_pertid+'.csv'
            tcs_df.to_csv(out_p)
        else:
            single_cell[lincs].append(_pertid)
    open('single_cell_dict','w').write(str(single_cell))
# transtcs_05:csv files derived from files in tcs_04, with 1/-1 representing the gene concordantly over/under-expressed treated by the compound across different cell lines
def TCSinLINCS(infolder_tcs04,outfolder_transtcs05):
    for tcs4_path in glob.glob(infolder_tcs04+'/*.csv'):
        pert_id = tcs4_path.split('/')[-1].split('.')[0]
        transtcs5_path = outfolder_transtcs05+'/'+pert_id+'.csv'
        tcs4_df = pd.read_csv(tcs4_path,index_col=0)
        index = tcs4_df.index.tolist()
        count_dic = {}
        up = ((tcs4_df==1).sum(axis=1)).tolist()
        dn = ((tcs4_df==-1).sum(axis=1)).tolist()
        for i,gid in enumerate(index):
            if up[i]>=dn[i]:
                count_dic[gid]=1
            else:
                count_dic[gid]=-1
        new_index = list(count_dic.keys())
        new_val = list(count_dic.values())
        outdf = pd.DataFrame({'id':new_index,'val':new_val})
        outdf.to_csv(transtcs5_path,index=False)

#innertcga_06: csv files derived from files in transtcs_05, containing only genes in both TCGA and LINCS-TCS
def TCSinnerGene(inpath_innerid,infolder_transtcs05,outfolder_innertcga06):
    inner_id_df = pd.read_csv(inpath_innerid,sep='\t')
    inner_ids = inner_id_df['LINCS_ID'].tolist()
    for path5 in glob.glob(infolder_transtcs05+'/*.csv'):
        pert_id = path5.split('/')[-1].split('.')[0]
        indf5 = pd.read_csv(path5)
        outdf = indf5[indf5['id'].isin(inner_ids)]
        tcs_num = len(outdf.index.tolist())
        outdf.to_csv(outfolder_innertcga06+'/'+pert_id+'_'+str(tcs_num)+'.csv',index=False)

#lauchedtcs_07: Copy files from innertcga_06 of compounds labeled 'lauched' in the https://clue.io/repurposing-app
def lauchedTCSinnerGene(inpath_lauched,infolder_innertcga06,outfolder_lauchedtcs07):
    lauched_df = pd.read_csv(inpath_lauched)
    lauched_ids = []
    for j in [7,11]:
        for i in lauched_df.index:
            lauched_data = lauched_df.iat[i,j]
            if isinstance(lauched_data,str):
                _data = lauched_data.strip().split(',')
                for each in _data:
                    _idlist = each.strip().split('-')
                    _id = '-'.join([_idlist[0],_idlist[1]])
                    if _id not in lauched_ids:
                        lauched_ids.append(_id)
    for comp_path in glob.glob(infolder_innertcga06+'/*.csv'):
        comp_fname = comp_path.split('/')[-1]
        comp_id = comp_fname.split('_')[0]
        if comp_id in lauched_ids:
            new_path = outfolder_lauchedtcs07+'/'+comp_fname
            shutil.copyfile(comp_path,new_path)
    print(len(os.listdir(outfolder_lauchedtcs07)))