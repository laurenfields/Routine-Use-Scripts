# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 15:10:33 2022

@author: lawashburn
"""

import csv
import pandas as pd
import os
import numpy as np

working_directory =r"D:\DIA\20221114_DSD_quant_exports\compiled" #Directory with output folders/files from PEAKS
output_dir = r"D:\DIA\20221114_DSD_quant_exports\compiled_IDs" #Output directory
file_type2 = '20221114_' #Folder name prefix
ppm_err = 20
logP = 37.6


def get_folder_names_with_strings(str_list):
    full_list = os.listdir(working_directory)
    final_list = [nm for ps in str_list for nm in full_list if ps in nm]
    return final_list

def get_file_names_with_strings(str_list):
    full_list = os.listdir(directory_path)
    final_list = [nm for ps in str_list for nm in full_list if ps in nm]
    return final_list


file_query2 = file_type2 #search query
folder_search = (get_folder_names_with_strings([file_query2]))

id_num_storage = []
exp_name_storage_final = []

for w in folder_search:
    peptide_name_storage = []
    exp_name_storage = []
    mz_storage = []
    scan_storage = []

    pep_file_loc_pep = []
    pep_file_loc_file = []
    pep_file_scan = []
    pep_file_mz = []

    exp_name = w
    exp_name_storage_final.append(exp_name)
    directory_path = working_directory + '\\' + w
    denovo_query = 'de novo only peptides.csv'
    denovo_search_file = (get_file_names_with_strings([denovo_query]))
    denovo_search_file_path = directory_path + '\\' + denovo_search_file[0]
    denovo = pd.read_csv(denovo_search_file_path)
    denovo_filter1 = denovo[denovo['Accession'].notna()]
    denovo_filter2 = denovo_filter1[denovo_filter1["Accession"].str.contains(":")==False]
    denovo_filter3 =  denovo_filter2[denovo_filter2['ppm'] <= ppm_err]
    denovo_filter4 =  denovo_filter3[denovo_filter3['ppm'] >= -(ppm_err)]
    #denovo_filter4 =  denovo_filter3[denovo_filter3['ppm'] >= -10]
    #denovo_filter4 =  denovo_filter4[denovo_filter4['-10lgP'] >= 42]
    denovo_peptides = denovo_filter4['Peptide'].values.tolist()
    denovo_source_files = denovo_filter4['Source File'].values.tolist()
    denovo_scans = denovo_filter4['Scan'].values.tolist()
    denovo_mz = denovo_filter4['m/z'].values.tolist()
    
    for a in denovo_peptides:
        pep_file_loc_pep.append(a)
    for b in denovo_source_files:
        pep_file_loc_file.append(b)
    for c in denovo_scans:
        pep_file_scan.append(c)
    for d in denovo_mz:
        pep_file_mz.append(d)
    
    db_search_query = 'DB search psm.csv'
    db_search_file = (get_file_names_with_strings([db_search_query]))
    db_search_file_path = directory_path + '\\' + db_search_file[0]
    db_search = pd.read_csv(db_search_file_path)
    db_filter1 = db_search[db_search['Accession'].notna()]
    db_filter2 = db_filter1[db_filter1["Accession"].str.contains(":")==False]
    db_filter3 =  db_filter2[db_filter2['ppm'] <= ppm_err]
    db_filter4 =  db_filter3[db_filter3['ppm'] >= -(ppm_err)]
    db_filter4 =  db_filter4[db_filter4['-10lgP'] >= logP]
    db_peptides = db_filter4['Peptide'].values.tolist()
    db_source_files = db_filter4['Source File'].values.tolist()
    db_scans = db_filter4['Scan'].values.tolist()
    db_mz = db_filter4['m/z'].values.tolist()
    
    for c in db_peptides:
        pep_file_loc_pep.append(c)
    for d in db_source_files:
        pep_file_loc_file.append(d)
    for c in db_scans:
        pep_file_scan.append(c)
    for d in db_mz:
        pep_file_mz.append(d)
        
    peptide_origin_glossary = pd.DataFrame()
    peptide_origin_glossary['Peptide'] = pep_file_loc_pep
    peptide_origin_glossary['Source File'] = pep_file_loc_file
    peptide_origin_glossary['Scan'] = pep_file_scan
    peptide_origin_glossary['m/z'] = pep_file_mz

    file_index = peptide_origin_glossary['Source File'].values.tolist()
    file_index_no_dups = []
    for j in file_index:
        if j not in file_index_no_dups:
            file_index_no_dups.append(j)

    number_files = len(file_index_no_dups)
    prot_pep_query = 'protein-peptides.csv'
    prot_pep_search_file = (get_file_names_with_strings([prot_pep_query]))
    prot_pep_search_file_path = directory_path + '\\' + prot_pep_search_file[0]
    prot_pep = pd.read_csv(prot_pep_search_file_path)
    prot_pep_filter1 =  prot_pep[prot_pep['ppm'] <= ppm_err]
    prot_pep_filter2 =  prot_pep_filter1[prot_pep_filter1['ppm'] >= -(ppm_err)]
    prot_pep_filter3 =  prot_pep_filter2[prot_pep_filter2['Unique'] == 'Y']
    prot_pep_filter3 =  prot_pep_filter3[prot_pep_filter3['-10lgP'] >= logP]
    peptide_IDs = prot_pep_filter3['Peptide'].values.tolist()
    
    peptide_IDs_formatted = []
    
    for e in peptide_IDs:
        e_status = []
        if e[1] == '.':
            e_start = e[2:]
            e_status.append(e_start)
        else:
            e_start = e
            e_status.append(e_start)

        if e_start[-2] == '.':
            e_finish = e_start[:-2]
            peptide_IDs_formatted.append(e_finish)
        else:
            e_finish = e_start
            peptide_IDs_formatted.append(e_finish)

    for f in peptide_IDs_formatted:
        peptide_origin_glossary_filtered =  peptide_origin_glossary[peptide_origin_glossary['Peptide'] == f]
        if len(peptide_origin_glossary_filtered) > 0:
            peptide_experiments = peptide_origin_glossary_filtered['Source File'].values.tolist()
            peptide_scans = peptide_origin_glossary_filtered['Scan'].values.tolist()
            peptide_mzs= peptide_origin_glossary_filtered['m/z'].values.tolist()
            peptide_name_storage.append(f)
            exp_name_storage.append(peptide_experiments)
            scan_storage.append(peptide_scans)
            mz_storage.append(peptide_mzs)
            
        else:
            pass

    id_table_report = pd.DataFrame()
    id_table_report['Peptide'] = peptide_name_storage
    id_table_report['Source File(s)'] = exp_name_storage
    id_table_report['m/z(s)'] = mz_storage
    id_table_report['Scan(s)'] = scan_storage
    id_table_report['Source File(s)'] = id_table_report['Source File(s)'].apply(lambda x: ','.join(map(str, x)))
    id_table_report['m/z(s)'] = id_table_report['m/z(s)'].apply(lambda x: ','.join(map(str, x)))
    id_table_report['Scan(s)'] = id_table_report['Scan(s)'].apply(lambda x: ','.join(map(str, x)))

    num_IDs = len(id_table_report)
    id_num_storage.append(num_IDs)
    out_path = output_dir + '\\' + exp_name + '_ID_report.csv'
    with open(out_path,'w',newline='') as filec:
                           writerc = csv.writer(filec)
                           id_table_report.to_csv(filec,index=False) 

id_report_counts = pd.DataFrame()
id_report_counts['Experiment'] = exp_name_storage_final
id_report_counts['Number of IDs'] = id_num_storage

out_path = output_dir + '\\all_IDs_summary.csv'
with open(out_path,'w',newline='') as filec:
                       writerc = csv.writer(filec)
                       id_report_counts.to_csv(filec,index=False) 