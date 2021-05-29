# -*- coding: utf-8 -*-
"""
Created on Tue May 29 02:03:54 2021

@author: KimMacdonald
"""

#created & tested on PC using:
# pandas v1.1.3
# python v3.8.5
#and on laptop with:
# pandas v1.0.5
# python v3.8.3

    
#import packages
import glob
import os
import subprocess
import fnmatch
import pandas as pd
import numpy as np
import re


# Store current date in a variable:
from datetime import datetime
Today = datetime.today().strftime('%Y-%m-%d')   # output is like '2021-01-26'  


#Read in All_Metadata.csv: (produced by ConnectToSQLserverDatabase_pyodbc.py script in https://github.com/Kim-Macdonald/ConnectToDatabase repository)
df_metadata1 = pd.read_csv("Path/To/All_Metadata.csv")
#print(df_metadata)


#---------------remove NaNs from df:------------------
df_metadata2 = df_metadata1.replace(np.nan, '')
#print (df_metadata2) 


#-----------Remove Unwanted Values from Ct fields only (not whole Df):-----------

# Values to remove (as list): 
tags = ['negative', 'Negative', 'NEGA', 'NEG', 'neg', 'Neg', 'missing info', 'UNDET', 'UNDE', 'SIMPLEXA', 'SGENE:', 'SEEGENE', 'GENEXPERT', 'GENEX', 'GENX', 'nan', 'NaN', 'E=NEG', 'E=', 'E:', 'E-', 'E', 'N=', 'N2=', 'N2:', 'N:', 'N', 'N2', 'S=', 'S:', 'S', 'R=', 'R:', 'R', 'CT=', 'CT', 'C:', 'ORF1AB=', 'ORF1AB:', 'ORF1=', 'ORF1:', 'ORF=', 'ORF:', 'ORF1AB', 'OF1AB:', 'ORF1', 'ORF', 'ORF1', 'ORF', 'RNA', 'NA', '\rN/A', '\rN/A (repeated)', 'na', 'R=N/A', 'E=N/A', '   @FAM E gene Positive  \(Ct. ', '\)', ' @FAM  gene Positive  \(Ct. ', '   @Quasar 670 N gene Positive  \(Ct. ', '   @Quasar 67  gene Positive  \(Ct. ', '   @Cal Red 610 RdRP gene Positive  \(Ct. ', '   @Cal ed 61 dP gene Positive  \(Ct. ', '\>40', '\>4', '-', 'UDT', '\/A', '\/A \(repeated', 'eg', 'G', 'GX', 'GXPT', 'IMPLXA', 'UD', ' \(repeated', 'X']
#This will remove the E and R etc from the CIDs in sample column as well (if apply to whole df and not certain columns)
#So only apply this to the Ct columns as below

# First format df as string, so the next str.replace method will work (otherwise you'll get an error):
df_metadata2 = df_metadata2.astype(str)

# Iterate through each item in list (stored in tags variable) to remove from df one at a time, to avoid error:
for tag in tags:
  df_metadata2[['ncov_qpcr_e_sarbeco_result','ncov_qpcr_rdrp_lee_result','ncov_qpcr_n2_result','ncov_qpcr_n_sarbeco_result','ncov_qpcr_orf1_result']] = df_metadata2[['ncov_qpcr_e_sarbeco_result','ncov_qpcr_rdrp_lee_result','ncov_qpcr_n2_result','ncov_qpcr_n_sarbeco_result','ncov_qpcr_orf1_result']].apply(lambda col: col.str.replace(tag, ''))
#print(df_metadata2)  
#Works for individual value - e.g. S= AND the list
#will remove partial matches as well as exact matches 
#(so be careful - use exact match removal only for the risky stuff - as below)

# Remove these exact values only (after the above are done) or will mess up Cts (e.g. will remove '.' from 29.4 so is 294) (or remove 0 from 40, so Ct becomes 4, etc)
# And only remove from Ct columns, or will replace 0 in index column etc. 
df_metadata2[['ncov_qpcr_e_sarbeco_result','ncov_qpcr_rdrp_lee_result','ncov_qpcr_n2_result','ncov_qpcr_n_sarbeco_result','ncov_qpcr_orf1_result']] = df_metadata2[['ncov_qpcr_e_sarbeco_result','ncov_qpcr_rdrp_lee_result','ncov_qpcr_n2_result','ncov_qpcr_n_sarbeco_result','ncov_qpcr_orf1_result']].replace('.', '')
df_metadata2[['ncov_qpcr_e_sarbeco_result','ncov_qpcr_rdrp_lee_result','ncov_qpcr_n2_result','ncov_qpcr_n_sarbeco_result','ncov_qpcr_orf1_result']] = df_metadata2[['ncov_qpcr_e_sarbeco_result','ncov_qpcr_rdrp_lee_result','ncov_qpcr_n2_result','ncov_qpcr_n_sarbeco_result','ncov_qpcr_orf1_result']].replace('0', '')
df_metadata2[['ncov_qpcr_e_sarbeco_result','ncov_qpcr_rdrp_lee_result','ncov_qpcr_n2_result','ncov_qpcr_n_sarbeco_result','ncov_qpcr_orf1_result']] = df_metadata2[['ncov_qpcr_e_sarbeco_result','ncov_qpcr_rdrp_lee_result','ncov_qpcr_n2_result','ncov_qpcr_n_sarbeco_result','ncov_qpcr_orf1_result']].replace('0.0', '')
# to do on whole df, for first line above, would just use:
#df_metadata2 = df_metadata2.replace('.', '')

#print(df_metadata2)



#-----------------Combine All Cts (in 5 columns) into 1 column:-------------

# Create new column for the combined Cts: 
df_metadata2.insert(4, "Ct_combo", "")
#print (df_metadata2)

# Copy the df into a new df, or python may start getting confused as to what df it's working on (metadata2 altered too many times) and give an error: 
df_metadata4 = df_metadata2.copy()
#print (df_metadata4[['ncov_qpcr_e_sarbeco_result', 'ncov_qpcr_rdrp_lee_result']])


# If e gene field not empty, use that value, otherwise take RdRp value if it's there, and if not take n2, then n, then orf1 last:
conditions = [
    df_metadata4['ncov_qpcr_e_sarbeco_result'].ne(''),
    df_metadata4['ncov_qpcr_rdrp_lee_result'].ne(''),
    df_metadata4['ncov_qpcr_n2_result'].ne(''),
    df_metadata4['ncov_qpcr_n_sarbeco_result'].ne(''),
    df_metadata4['ncov_qpcr_orf1_result'].ne(''),
]

# Tells it to use the e gene result if it's not empty (as defined above), then the RdRp value (if rdrp field not empty), then n2, then n, then orf1 last:
choices=[df_metadata4['ncov_qpcr_e_sarbeco_result'], df_metadata4['ncov_qpcr_rdrp_lee_result'], df_metadata4['ncov_qpcr_n2_result'], df_metadata4['ncov_qpcr_n_sarbeco_result'], df_metadata4['ncov_qpcr_orf1_result']]

# Fill the Ct_combo field with values as detailed above (combined Cts):
df_metadata4['Ct_combo'] = np.select(conditions, choices, default='')

# Check that Ct_combo field is combining values from fields properly:
#print (df_metadata4['Ct_combo']) 
# Will have NaNs in whole df (if don't take care of these at start)


# #------------------Save Output:----------------- (only use if need to test that you're getting what you want each step along the way)
# df_metadata4.to_csv("Path/To/All_Metadata_cleaned.csv", index=False)

# #-----------------------------------------------



#---------------Prepare Cleaned metadata df with combined Cts:---------------------
df_CleanMetadata1 = df_metadata4.copy()

# drop extra index:
df_CleanMetadata2 = df_CleanMetadata1.drop(columns=df_CleanMetadata1.columns[0])

# # Use this to list columns in a file, when writing code above etc, so you don't have to type out
# df_colList = list(df_CleanMetadata1.columns)
# print(df_colList)

#remove NaNs from df:
df_CleanMetadata2 = df_CleanMetadata2.replace(np.nan, '')

# Define column order of metadata df:
df_CleanMetadata2 = df_CleanMetadata2[['containerid', 'second_containerid', 'seq_containerid', 'Ct_combo', 'collection_date', 'ncov_qpcr_e_sarbeco_result', 'ncov_qpcr_rdrp_lee_result', 'ncov_qpcr_n_sarbeco_result', 'ncov_qpcr_n2_result', 'ncov_qpcr_orf1_result']]
print(df_CleanMetadata2)


#-----------------Read in Run Summary w updated Lineages:---------- (result file from UpdateQCsummaryLineages_v2.py at https://github.com/Kim-Macdonald/LineageUpdater)

# With Newest Lineages:
# path = "Path/To/Runs_CombinedQCsummary_PlusNewestLineages_*.xlsx"

# for filename in glob.glob(path):
#     with open(filename, 'r') as f:
#         print(f)
#         QCsummary0a = pd.read_csv(f)


#-----------------OR Read in Run Summary without updated Lineages (and run lineage updater script after this): (result file from addVoCcalls_RunNum_v2b.py at https://github.com/Kim-Macdonald/VoCcaller)

# Without Newest Lineages: 
# Read in QCsummary excel file (Sheet1):
df_QCsummary0 = pd.read_excel('Path/To/Runs_CombinedQCsummary.xlsx', sheet_name=0)
#print(df_QCsummary0)


#----------Merge QCsummary with CleanMetadata (to combine CIDs later):-------------

# Some samples will have sample IDs that match the primary containerid field. Some will match the secondary_containerid field. Some will match the seq_containerid field. 
# So need to match fastq sample IDs to all 3 fields containerids for metadata matching (b/c can be from any one of the 3)

# Merge sample (CID) with primary CID from plover (in clean metadata)
df_QCsummary_CleanMetadata_merge1 = pd.merge(df_QCsummary0, df_CleanMetadata2, how='left', left_on='sample', right_on='containerid')
# only keep columns needed (matched by sample and CID so keep those, and Ct_combo, and coll date:)
df_QCsummary_CleanMetadata_merge1 = df_QCsummary_CleanMetadata_merge1[['sample','containerid', 'Ct_combo', 'collection_date_y']]
# Rename the containerid col to CID (and same for other 2 merges below), so can append dfs with same column names (or will add extra columns):
df_QCsummary_CleanMetadata_merge1 = df_QCsummary_CleanMetadata_merge1.rename(columns={'containerid':'CID'})
#print(df_QCsummary_CleanMetadata_merge1) #works (46 columns, inc index)

# Merge sample (CID) with secondary CID from plover (in clean metadata)
df_QCsummary_CleanMetadata_merge2 = pd.merge(df_QCsummary0, df_CleanMetadata2, how='left', left_on='sample', right_on='second_containerid')
# only keep columns needed (matched by sample and 2ndary CID so keep those, and Ct_combo, and coll date:)
df_QCsummary_CleanMetadata_merge2 = df_QCsummary_CleanMetadata_merge2[['sample', 'second_containerid', 'Ct_combo', 'collection_date_y']]
# Rename the second_containerid col to CID, so can append dfs with same column names (or will add extra columns):
df_QCsummary_CleanMetadata_merge2 = df_QCsummary_CleanMetadata_merge2.rename(columns={'second_containerid':'CID'})
#print(df_QCsummary_CleanMetadata_merge2) #works (46 columns, inc index)

# Merge sample (CID) with seq CID from plover (in clean metadata)
df_QCsummary_CleanMetadata_merge3 = pd.merge(df_QCsummary0, df_CleanMetadata2, how='left', left_on='sample', right_on='seq_containerid')
# only keep columns needed (matched by sample and seq CID so keep those, and Ct_combo, and coll date:)
df_QCsummary_CleanMetadata_merge3 = df_QCsummary_CleanMetadata_merge3[['sample', 'seq_containerid', 'Ct_combo', 'collection_date_y']]
# Rename the seq_containerid col to CID, so can append dfs with same column names (or will add extra columns):
df_QCsummary_CleanMetadata_merge3 = df_QCsummary_CleanMetadata_merge3.rename(columns={'seq_containerid':'CID'})
#print(df_QCsummary_CleanMetadata_merge3) #works (46 columns, inc index)


     #------------Merge the 3 merged dfs (above) together (append)--------

# Append the 3 dfs (add to bottom of last, as samples likely won't match in 3). Then can sort to remove duplicates:
df_QCsummary_CleanMetadata_merge4 = df_QCsummary_CleanMetadata_merge1.append(df_QCsummary_CleanMetadata_merge2, ignore_index=True)

df_QCsummary_CleanMetadata_merge5 = df_QCsummary_CleanMetadata_merge4.append(df_QCsummary_CleanMetadata_merge3, ignore_index=True)
#print(df_QCsummary_CleanMetadata_merge5)

#Drop extra index column:
df_QCsummary_CleanMetadata_merge6 = df_QCsummary_CleanMetadata_merge5.reset_index(drop=True)
#print(df_QCsummary_CleanMetadata_merge6)

# Convert sample column to str, or get error when sort:
df_QCsummary_CleanMetadata_merge6 = df_QCsummary_CleanMetadata_merge6.astype(str)
# Sort to remove duplicate sample IDs, then sort by CID match (so CID matches sort above no match - so when keep first of duplicates, it keeps the matched sample data instead of throwing it out)
df_QCsummary_CleanMetadata_merge7 = df_QCsummary_CleanMetadata_merge6.sort_values(by=['sample', 'CID'], ascending = (True, True))  
#Save this to file to check sorting is correct: (it is, matches are at top)
#df_QCsummary_CleanMetadata_merge7.to_csv("Path/To/All_Metadata_cleaned_SORTEDwithDUPS_MatchedToSampleIDs.csv", index=False)
# Remove duplicate sampleIDs in merged metadata file:
df_QCsummary_CleanMetadata_merge7 = df_QCsummary_CleanMetadata_merge7.drop_duplicates(['sample'], keep='first') 
#print(df_QCsummary_CleanMetadata_merge7)
#This took the df from 104407 rows to 34769 rows (the original # rows in the Run summary file - which is correct)


# Don't need CID column anymore - can leave as sample
# Only need: 'sample', 'Ct_combo', 'collection_date_y' fields for metadata.tsv
df_QCsummary_CleanMetadata_merge8 = df_QCsummary_CleanMetadata_merge7[['sample', 'Ct_combo', 'collection_date_y']]
#print(df_QCsummary_CleanMetadata_merge8)
# Then change column names to sample, ct, date (for proper column naming in metadata.tsv):
df_QCsummary_CleanMetadata_merge8 = df_QCsummary_CleanMetadata_merge8.rename(columns={'Ct_combo':'ct', 'collection_date_y':'date'})
#print(df_QCsummary_CleanMetadata_merge9)

#------This df has NaNs - remove before saving final output----
#remove NaNs from df:
df_QCsummary_CleanMetadata_merge8 = df_QCsummary_CleanMetadata_merge8.replace(np.nan, '')
#print(df_QCsummary_CleanMetadata_merge8)

#nan is actually a string now for some reason, so np.nan above won't work. 
#So just remove the string 'nan'
df_QCsummary_CleanMetadata_merge9 = df_QCsummary_CleanMetadata_merge8.replace('nan', '')
#print(df_QCsummary_CleanMetadata_merge9)

#------------------Save Output:-----------------
#----------Save CleanedMetadata with Matched sampleID/CID column and Cts and Coll Dates----------

df_QCsummary_CleanMetadata_merge9.to_csv("Path/To/All_Metadata_cleaned_MatchedToSampleIDs.csv", index=False)

