import pandas as pd
import numpy as np
import sys, shutil,os

# /persist_disk/mocks_and_isolates/eu_20220817-120350_mocks_and_isolates/isolates/220816-210858/220816-210858
#python a.py 220816-210858.k2.report 0.01 220816-210858.k2.out  ../190320035227.2.1.1/190320035227.2.1.1.R1.filtered.fa ../190320035227.2.1.1/190320035227.2.1.1.R2.filtered.fa 190320035227.2.1.1

def report_parser(x, y):
    """
    Parse reports and construct a list for mtb tabel backbone
    """
    # ~ self.LOGGER.info(f"Parsing {self.report_file} file")
    l = []
    with open(x) as fin:
        for line in fin:
            if int(line.split('\t')[2]) >= 3 and int(line.split('\t')[4]) not in [0,1]:
                taxid = int(line.split('\t')[4])
                counts = int(line.split('\t')[2])
                out = f'{taxid}\t{counts}'
                l.append(out.split('\t'))
                    # ~ self.info_list.append(out.split('\t'))
    
    df = pd.DataFrame(l, columns=['taxid', 'counts'])
    df['taxid'] = df['taxid'].astype(int)
    df['counts'] = df['counts'].astype(int)
    df['perc']= (df['counts']/df['counts'].sum())*100
    # ~ print(df)
    taxid_list = df['taxid'].loc[df['perc'] <= float(y)].to_list()
    #add not hit no rank 
    taxid_list.append(0)
    taxid_list.append(1)
    # ~ print(taxid_list)
    return taxid_list

tax_list_bellow_perc_cutoff = report_parser(sys.argv[1], sys.argv[2])

def get_read_id(z):

    tax_list_bellow_perc_cutoff = report_parser(sys.argv[1], sys.argv[2])

    info_list =[]
    with open(z) as fin:
        for line in fin:
            _id = line.split('\t')[1]
            tax = line.split('\t')[2]
           
            for i in   tax_list_bellow_perc_cutoff:
                if int(i) == int(tax):
                    out =f'{tax}\t{_id}'
                    
                    info_list.append(out.split('\t'))
                     
    df = pd.DataFrame(info_list, columns=['original_taxid', 'read_id'])
    #df_outname = sys.argv[4].replace("R1.filtered.fa", 'id_equivalence.tsv')
    #df.to_csv(df_outname, sep='\t', index=False)
    return df


def get_fastq_from_id(r1_fq,r2_fq):
        
        r1_dict, r2_dict = {}, {}
        
        with open(r1_fq) as f:
            lines=f.readlines()
        headR1 = [item[:-1] for item in lines[::4]]
        readR1 = [item[:-1] for item in lines[1::4]]
        qualR1 = [item[:-1] for item in lines[3::4]]

        
        for h,r,q in zip(headR1, readR1,qualR1):
            h = h[1:].split(' ')[0]
            r = r.strip('\n')
            q = q.strip('\n')
            r1_dict[h] = [r,q]
        
        with open(r2_fq) as f:
            lines=f.readlines()
        headR2 = [item[:-1] for item in lines[::4]]
        readR2 = [item[:-1] for item in lines[1::4]]
        qualR2 = [item[:-1] for item in lines[3::4]]

        
        for h,r,q in zip(headR2, readR2,qualR2):
            h = h[1:].split(' ')[0]
            r = r.strip('\n')
            q = q.strip('\n')
            r2_dict[h] = [r,q]
       
        return r1_dict, r2_dict
            
def generate_fq_files():
    
    tax = get_read_id(sys.argv[3])
    
    out_r1_fname = sys.argv[4].replace('.R1.fq','.R1.from_k2.fq')
    out_r2_fname = sys.argv[5].replace('.R2.fq','.R2.from_k2.fq')
    
    
    out_r1 = open(out_r1_fname, 'w')
    out_r2 = open(out_r2_fname, 'w')
    

    
    r1,r2 = get_fastq_from_id(sys.argv[4],sys.argv[5])
    
    for index, row in tax.iterrows():
        
        #print(row['read_id'], r1.get(row['read_id'])[0], r1.get(row['read_id'])[1])
        
        _r1_head = '@' + row['read_id'] + '_r1'
        _r1_seq  = r1.get(row['read_id'])[0]
        _r1_qual = r1.get(row['read_id'])[1]
        
        _r1 = f'{_r1_head}\n{_r1_seq}\n+\n{_r1_qual}\n'
        
        _r2_head = '@' + row['read_id'] + '_r2'
        _r2_seq  = r2.get(row['read_id'])[0]
        _r2_qual = r2.get(row['read_id'])[1]
        
        _r2 = f'{_r2_head}\n{_r2_seq}\n+\n{_r2_qual}\n'

        
        
        out_r1.writelines(_r1)
        out_r2.writelines(_r2)
        
        




generate_fq_files()


