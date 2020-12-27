'''Genes differentially expressed between tumor and adjacent normal samples'''
import pandas as pd
import numpy as np
from scipy import stats
import os, glob, shutil, gzip, json, math, operator

#orig1_path: a folder of RNA-Seq files downloaded from TCGA
#json_path: the json file downloaded from TCGA

'''1.Unzip the RNA-seq files downloaded from TCGA (orig1_path) and copy them to a new folder(decompress_path)'''
def decompressGzip(orig1_path,decompress_path):
    os.mkdir(decompress_path)
    for orig2_path in glob.glob(orig1_path+'/*'):
        for gz_path in glob.glob(orig2_path+'/*.gz'):
            gz_fname = gz_path.split('/')[-1]
            txt_fname = gz_fname.replace('.gz','')
            txt_path = decompress_path+'/'+txt_fname
            gz_file = gzip.GzipFile(gz_path)
            open(txt_path,'wb').write(gz_file.read())
            gz_file.close()
    print(len(os.listdir(decompress_path)))

'''2.Rename the expression files in the format of 'TCGA-***' using the metadata.json file'''
def renameFile(json_path,decompress_path,rename_path):
    os.mkdir(rename_path)
    with open(json_path) as infile:
        meta_json = json.load(infile)
        old_new = {each['file_name'][:-3]:each['associated_entities'][0]['entity_submitter_id'] for each in meta_json}
        for old_path in glob.glob(decompress_path+'/*.txt'):
            old_fname = old_path.split('/')[-1]
            new_fname = old_new[old_fname]+'.txt'
            new_path = rename_path+'/'+new_fname
            shutil.copyfile(old_path,new_path)

'''3.Clsssify the files into tumor, normal and replicate categories'''
def classify(rename_path,classify_path):
    os.mkdir(classify_path),
    os.mkdir(classify_path+'/tumor')  #01
    os.mkdir(classify_path+'/normal') #11
    os.mkdir(classify_path+'/replicate')
    tumor_path = classify_path+'/tumor'
    normal_path = classify_path+'/normal'
    rep_path = classify_path+'/replicate'
    tumor_l, normal_l, replicate_l = [], [], []
    for i in os.listdir(rename_path): #TCGA-ZG-A9NI-01A-11R-A41O-07.txt
        if i[13:15] == '01':          #tumor
            for j in os.listdir(rename_path):
                if j[13:15] == '11' and j[:13] == i[:13] and j not in normal_l:
                    normal_l.append(j)
                    tumor_l.append(i)
                    tumor_oldpath = rename_path+ '/' + i
                    tumor_newpath = tumor_path + '/' + i
                    normal_oldpath = rename_path+ '/' + j
                    normal_newpath = normal_path+ '/' + j
                    shutil.copyfile(tumor_oldpath, tumor_newpath)
                    shutil.copyfile(normal_oldpath, normal_newpath)
                elif j[13:15] == '11' and j[:13] == i[:13] and j in normal_l:
                    rep_oldpath = rename_path+ '/'+ j
                    rep_newpath = rep_path + '/' + j
                    replicate_l.append(j)
                    shutil.copyfile(rep_oldpath, rep_newpath)
                else:
                    continue
    print('normal_l:%d'% len(normal_l))
    print('tumor_l:%d'% len(tumor_l))
    print('replicate_l:%d'% len(replicate_l))
    print('all_l:%d'% len(os.listdir(rename_path)))

'''4.Merge the tumor and normal file of the same patient'''
def concatTN(classify_path,concat_path):
    num = 1
    concat_df= pd.DataFrame()
    for tumor_path in glob.glob(classify_path+'/tumor/*.txt'):
        tumor_fname = tumor_path.split('/')[-1] #TCGA-J4-A83J-01A-11R-A36G-07.txt
        tsample_name = tumor_fname[:-4]
        tumor_df = pd.read_csv(tumor_path, sep='\t',header=None,index_col=None,names=['Gene_ID',tsample_name])
        partner_tumor = tumor_fname[:13] #TCGA-J4-A83J-
        if num == 1:
            concat_df = pd.concat([concat_df,tumor_df], axis=1)
        else:
            concat_df = pd.concat([concat_df, tumor_df[tsample_name]], axis=1)
        num += 1
        for normal_path in glob.glob(classify_path + '/normal/*.txt'):
            normal_fname = normal_path.split('/')[-1]
            partner_normal = normal_fname[:13]
            if partner_normal == partner_tumor:
                nsample_name = normal_fname[:-4]
                normal_df = pd.read_csv(normal_path, sep='\t',header=None,index_col=None,names=['Gene_ID',nsample_name])
                concat_df = pd.concat([concat_df, normal_df[nsample_name]], axis=1)
                num += 1
    concat_df.to_csv(concat_path,index = False)

'''5.Drop genes which are null in more than half of the samples, or absent from the LINCS fiels, and concert the gene id into the format without version number'''
def filterNull(concat_path, filternull_path):
    concat_df = pd.read_csv(concat_path)
    cutoff = 0.5*(len(concat_df.columns.tolist())-1)
    count = list((concat_df==0).sum(axis=1).items())
    keep = []
    for i in count:
        if i[1] < cutoff:
            keep.append(i[0])
    print('number of genes that has expressed values:%d'% len(keep))
    fliter_null_df = concat_df.iloc[keep, :]
    fliter_null_df.to_csv(filternull_path, index = False)

def innerLINCS(filternull_path,innerid_path,innerLINCS_path):
    innerids = pd.read_csv(innerid_path,sep='\t')['TCGA_ID'].tolist()
    orig_df = pd.read_csv(filternull_path)
    with open(innerLINCS_path,'w') as outfile:
        for ix,line in enumerate(open(filternull_path,'r').readlines()):
            if ix == 0:
                outfile.write(line)
            else:
                _list = line.strip().split(',')
                _id = _list[0].strip().split('.')[0]
                if _id in innerids:
                    del(_list[0])
                    _list.insert(0,_id.strip().split('.')[0])
                    outline = ','.join(_list)+'\n'
                    outfile.write(outline)
                else:
                    continue

'''6.Calculate FoldChange, log2FoldChange, P_value, rank, FDR'''
def compute(innerLINCS_path,compute_path1):
    indf = pd.read_csv(innerLINCS_path)
    gene_num = len(indf.index.tolist())
    sample_num = len(indf.columns.tolist())-1
    pair_num = sample_num/2
    print('gene_number: %d'% gene_num)
    print('pair_number: %d'% pair_num)
    even, odd = [], []
    for i in range(0, int(sample_num)):
        if i % 2 != 0:
            odd.append(i)
            even.append(i + 1)
    ids, fcs, logs, ps = [], [], [], []
    for idx in indf.index:
        name = indf.iat[idx, 0]
        sumfc = 0
        for j in odd:
            sumfc += (indf.iat[idx, j] + 0.000001) / (indf.iat[idx, j + 1] + 0.000001)
        fc = sumfc / pair_num
        log = math.log2(fc)
        p = stats.ttest_ind(indf.iloc[idx, odd], indf.iloc[idx, even])[1]
        ids.append(name)
        fcs.append(fc)
        logs.append(log)
        ps.append(p)
    data1 = {'Gene_name': ids, 'FC': fcs, 'log2FC': logs, 'P_value': ps}
    out1 = pd.DataFrame(data=data1)
    out1.to_csv(compute_path1, index=False)

def rankPvalue(compute_path1,rank_path2):
    indf1 = pd.read_csv(compute_path1)
    df1 = indf1.sort_values(by='P_value', ascending=True).reset_index()
    ranks = []
    for r in df1.index:
        ranks.append(r + 1)
    datar = {'rank': ranks}
    df_rank = pd.DataFrame(data=datar)
    out_r1 = pd.concat([df1, df_rank], axis=1)
    out_r2 = out_r1.drop(['index'], axis=1)
    out_r2.to_csv(rank_path2, index=False)

def computeFDRinnerLINCS_path,rank_path2,fdr_path3):
    orig_df = pd.read_csv(innerLINCS_path)
    gene_num = len(orig_df.index.tolist())
    before = pd.read_csv(rank_path2)
    fdrs = []
    for ix,line in enumerate(open(rank_path2,'r').readlines()):
        if ix != 0:
            line = line.strip().split(',')
            fdr = float(line[3]) * gene_num / int(line[4])
            fdrs.append(fdr)
            data_fdr = {'FDR': fdrs}
            df_fdr = pd.DataFrame(data=data_fdr)
            out_fdr = pd.concat([before,df_fdr],axis = 1)
    out_fdr.to_csv(fdr_path3,index = False)

'''7.Screen differentially expressed genes'''
def chooseGene(fdr_path3,choose_path):
    indf = pd.read_csv(fdr_path3)
    filter1 = indf[(indf['log2FC'] > 1) | (indf['log2FC'] < -1)]
    filter2 = filter1[filter1['FDR'] < 0.1]
    filter2.to_csv(choose_path, index=False)
    print('number of remain genes: %d'%len(filter2.index.tolist()))

'''8.Convert the expression value into 1/-1 according whether the gene is over-/down-expressed
  Convert Ensemble ID to Entrezid ID, and csv file to txt'''
def transTCGA(choose_path,detransTCGA_path):
    indf = pd.read_csv(choose_path)
    with open(detransTCGA_path, 'w') as outfile:
        outfile.write('TCGA_ID\tvalue\n')
        for i in indf.index:
            outlist = []
            name = indf.iat[i, 0]
            log = indf.iat[i, 2]
            if log > 0:
                _log = 1
            else:
                _log = -1
            outlist.append(name)
            outlist.append(str(_log))
            outline = '\t'.join(outlist) + '\n'
            outfile.write(outline)

def transID(innerid_path,detransTCGA_path,detransLINCS_path):
    id_metadf = pd.read_csv(innerid_path,sep='\t')
    before_deg = pd.read_csv(detransTCGA_path,sep='\t')
    
    lincs_list,val_list = [],[]
    for idx in before_deg.index:
        tcga_id = before_deg.iat[idx,0]
        value = before_deg.iat[idx,1]

        lincs_line = id_metadf['LINCS_ID'][id_metadf['TCGA_ID']==tcga_id]
        lincs_id = str(lincs_line).strip().split('\n')[0].split('    ')[-1]
        lincs_list.append(lincs_id)
        val_list.append(value)
        
        after_df = pd.DataFrame(data = {'val':val_list},index = pd.Index(lincs_list,name = 'id'))
        after_df.to_csv(detransLINCS_path)