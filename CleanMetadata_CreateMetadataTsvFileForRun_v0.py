# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 18:48:28 2021

@author: kmacdonald
"""

#created & tested on PC using:
# pandas v1.1.3
# python v3.8.5
#and on laptop with:
# pandas v1.0.5
# python v3.8.3
#Server has pandas v1.2.3 and python v3.8.8
    
#import packages I need
import glob
import os
import subprocess
import fnmatch
import pandas as pd
import numpy as np
import re

#-----Store things in variables for later----------
# Store current date in a variable:
from datetime import datetime
Today = datetime.today().strftime('%Y-%m-%d')   # output is like '2021-01-26'  

#save the current working directory (cwd) to a variable to use in everything below. 
#For us, this would be the MiSeqRunID directory (in analysis_by_run) for each run - it changes each time we analyze a different run, so i want to pull this from where ever I am, so I don't have to enter it each time:
cwdPath = os.getcwd()

#Define variable using last part of directory/path (MiSeqRunID directory name)
#This will be used to name your files uniquely:
MiSeqRunID = os.path.basename(os.path.normpath(cwdPath))
#print(MiSeqRunID)



#---------------Read in All_Metadata.csv (produced by another script (scheduled)):--------------
#Another 2 scheduled scripts will pull the All_metadata.csv file from plover db (at 10pm) and also transfer it to server at 11pm daily. 
df_metadata1 = pd.read_csv("/path/To/analysis_by_run/All_Metadata.csv", dtype=object)
#print(df_metadata1)

#Have to declare dtypes on import above or get DtypeWarning that columns 5,6,7 have mixed types, and you need to specify dtype option on import or set low_memory = False (which does nothing)
#python will use a bunch of memory trying to poorly guess at dtypes for the columns otherwise, and slow things down.

#---------------remove NaNs from df:------------------
df_metadata2 = df_metadata1.replace(np.nan, '')
#print (df_metadata2) 


#-----------Remove Unwanted Values from Ct fields only (not whole Df):-----------

# Values to remove (as list): 
tags = ['   @FAM E gene Positive  \(Ct. ', '   @FAM E gene Positive  (Ct. ', ' @FAM  gene Positive  \(Ct. ', '   @FAM  gene Positive  (Ct. ', '   @FAM  gene Positive  \(Ct. ', '   @Quasar 670 N gene Positive  \(Ct. ', '   @Quasar 670 N gene Positive  (Ct. ', '   @Quasar 670  gene Positive  \(Ct. ', '   @Quasar 670  gene Positive  (Ct. ', '   @Quasar 67  gene Positive  \(Ct. ', '   @Quasar 67  gene Positive  (Ct. ', '   @Quasar 670  gene Positive  (Ct. ', '   @Cal Red 610 RdRP gene Positive  \(Ct. ', '   @Cal Red 610 RdRP gene Positive  (Ct. ', '   @Cal ed 61 dP gene Positive  \(Ct. ', '   @Cal ed 61 dP gene Positive  (Ct. ', '   @Cal ed 610 dP gene Positive  \(Ct. ', '   @Cal ed 610 dP gene Positive  (Ct. ', 'negative', 'Negative', 'NEGA', 'NEG', 'neg', 'Neg', 'missing info', 'Undetermined', 'UNDET', 'UNDE', 'SIMPLEXA', 'SGENE:', 'SEEGENE', 'GENEXPERT', 'GENEX', 'GENX', 'nan', 'NaN', 'E=NEG', 'E=', 'E:', 'E-', 'E', 'N=', 'N2=', 'N2:', 'N:', 'N', 'N2', 'S=', 'S:', 'S', 'ORF1AB=', 'ORF1AB:', 'ORF1=', 'ORF1:', 'ORF=', 'ORF:', 'ORF1AB', 'OF1AB:', 'ORF1', 'ORF', 'ORF1', 'ORF', 'R=', 'R:', 'R', 'CT=', 'CT', 'C:', 'RNA', 'NA', '\rN/A', '\rN/A (repeated)', '\r/A (repeated)', 'na', '\r/A', 'R=N/A', 'E=N/A', '/A', '/A (repeated)', ' (repeated)', '>40', '\>40', '\>4', '-', 'UDT', '\/A', '\/A \(repeated', 'eg', 'G', 'GX', 'GXPT', 'IMPLXA', 'UD', ' \(repeated', 'X', 'AL', 'LAT  CUV', 'WAB', '\)',  ')']
#This will remove the E and R etc from the CIDs in sample/CID columns, So only apply this to the Ct columns as below

# First format df as string, so the next str.replace method will work (otherwise you'll get an error):
df_metadata2 = df_metadata2.astype(str)

# Iterate through each item in list (stored in tags variable) to remove from df one at a time, to avoid error:
for tag in tags:
  df_metadata2[['ncov_qpcr_e_sarbeco_result','ncov_qpcr_rdrp_lee_result','ncov_qpcr_n2_result','ncov_qpcr_n_sarbeco_result','ncov_qpcr_orf1_result']] = df_metadata2[['ncov_qpcr_e_sarbeco_result','ncov_qpcr_rdrp_lee_result','ncov_qpcr_n2_result','ncov_qpcr_n_sarbeco_result','ncov_qpcr_orf1_result']].apply(lambda col: col.str.replace(tag, '', regex=False))
#print(df_metadata2)  
#WORKS for individual value (e.g. S=) AND the list
#will remove partial matches as well as exact matches 
#(so be careful - use exact match removal only for the risky stuff - as below)

# Remove these exact values only (after the above are done) or will mess up Cts (e.g. will remove '.' from 29.4 so is 294) (or remove 0 from 40, so Ct becomes 4, etc)
# And only remove from Ct columns, or will replace 0 in index column etc. 
df_metadata2[['ncov_qpcr_e_sarbeco_result','ncov_qpcr_rdrp_lee_result','ncov_qpcr_n2_result','ncov_qpcr_n_sarbeco_result','ncov_qpcr_orf1_result']] = df_metadata2[['ncov_qpcr_e_sarbeco_result','ncov_qpcr_rdrp_lee_result','ncov_qpcr_n2_result','ncov_qpcr_n_sarbeco_result','ncov_qpcr_orf1_result']].replace('.', '')
df_metadata2[['ncov_qpcr_e_sarbeco_result','ncov_qpcr_rdrp_lee_result','ncov_qpcr_n2_result','ncov_qpcr_n_sarbeco_result','ncov_qpcr_orf1_result']] = df_metadata2[['ncov_qpcr_e_sarbeco_result','ncov_qpcr_rdrp_lee_result','ncov_qpcr_n2_result','ncov_qpcr_n_sarbeco_result','ncov_qpcr_orf1_result']].replace('0', '')
df_metadata2[['ncov_qpcr_e_sarbeco_result','ncov_qpcr_rdrp_lee_result','ncov_qpcr_n2_result','ncov_qpcr_n_sarbeco_result','ncov_qpcr_orf1_result']] = df_metadata2[['ncov_qpcr_e_sarbeco_result','ncov_qpcr_rdrp_lee_result','ncov_qpcr_n2_result','ncov_qpcr_n_sarbeco_result','ncov_qpcr_orf1_result']].replace('0.0', '')
# to do on whole df, for first line above would just use:
#df_metadata2 = df_metadata2.replace('.', '')
#print(df_metadata2)



#-----------------Combine All Cts into 1 field:-------------

# Create new column for the combined Cts: 
df_metadata2.insert(4, "Ct_combo", "")
#print (df_metadata2)

# Copy the df into a new df, or python may start getting confused as to what df it's working on (metadata2 altered too many times) and give an error: 
df_metadata4 = df_metadata2.copy()
#print (df_metadata4[['ncov_qpcr_e_sarbeco_result', 'ncov_qpcr_rdrp_lee_result']])


# If e gene field not empty, use that value, otherwise take RdRp value if it's there, and if not take n2, then n, then orf1 last:
conditions1 = [
    df_metadata4['ncov_qpcr_e_sarbeco_result'].ne(''),
    df_metadata4['ncov_qpcr_rdrp_lee_result'].ne(''),
    df_metadata4['ncov_qpcr_n2_result'].ne(''),
    df_metadata4['ncov_qpcr_n_sarbeco_result'].ne(''),
    df_metadata4['ncov_qpcr_orf1_result'].ne(''),
]

# Tells it to use the e gene result if it's not empty (as defined above), then the RdRp value (if rdrp field not empty), then n2, then n, then orf1 last:
choices1 = [df_metadata4['ncov_qpcr_e_sarbeco_result'], df_metadata4['ncov_qpcr_rdrp_lee_result'], df_metadata4['ncov_qpcr_n2_result'], df_metadata4['ncov_qpcr_n_sarbeco_result'], df_metadata4['ncov_qpcr_orf1_result']]

# Fill the Ct_combo field with values as detailed above (combined Cts):
df_metadata4['Ct_combo'] = np.select(conditions1, choices1, default='')

# Check that Ct_combo field is combining values from fields properly:
#print (df_metadata4[['Ct_combo', 'ncov_qpcr_e_sarbeco_result', 'ncov_qpcr_rdrp_lee_result']]) 
# Will have NaNs in whole df (if don't take care of these at start)


# #------------------Format Cleaned metadata dfs:-----------------

df_CleanMetadata1 = df_metadata4.copy()

# Drop extra index row:
df_CleanMetadata2 = df_CleanMetadata1.drop(columns=df_CleanMetadata1.columns[0])
#print(df_CleanMetadata2)

#remove NaNs from df:
df_CleanMetadata2 = df_CleanMetadata2.replace(np.nan, '')

# # Use this to list columns in a file, when writing code below etc, so you don't have to type out column names
# df_colList = list(df_CleanMetadata1.columns)
# print(df_colList)

# Define column order of metadata df:
df_CleanMetadata2 = df_CleanMetadata2[['containerid', 'second_containerid', 'seq_containerid', 'Ct_combo', 'collection_date', 'ncov_qpcr_e_sarbeco_result', 'ncov_qpcr_rdrp_lee_result', 'ncov_qpcr_n_sarbeco_result', 'ncov_qpcr_n2_result', 'ncov_qpcr_orf1_result']]
#print(df_CleanMetadata2)


# Save All_metadata_cleaned.csv for UpdatedLineages script to merge with RunSummary for All Runs (done later by another script to backfill any missing/old metadata):
# saves in analysis_by_run directory (so is easily available for other scripts):
df_CleanMetadata2.to_csv("/path/to/analysis_by_run/All_Metadata_cleaned.csv", index=False)



#---------Create list of FastqIDs to match to cleaned metadata-------

#Run bash commands to generate the FastqList (assumes fastqs are stored in a separate directory from analysis, as below):
bashCommand1 = "ls ../../direct_fastq_symlinks_by_run/" + MiSeqRunID + "/*.fastq.gz | cut -d/ -f5 | cut -d. -f1 | cut -d_ -f1 | sort -u | sed '/Undetermined/d' >FastqList.txt"

subprocess.run(bashCommand1, shell=True, check=True)


#-----------------Read in FastqID list:----------
#Read in lists of sampleIDs:
df_FastqList1 = pd.read_table("FastqList.txt", header=None)
#print(df_FastqList1)


#Add column headers to match the other file:
#Rename Column 0 as 'sample'
df_FastqList2 = df_FastqList1.rename(columns={0: "sample"})
#print(df_FastqList2)


#-----------Create CID column (From fastqIDs) and Fill------------

df_FastqList2.insert(1, "CID1", "NA")
#print (df_FastqList2)

# Parse SampleID out of new fastq sampleIDs:
LibNum_split_QC = df_FastqList2['sample'].str.split("-")
#print(LibNum_split_QC)
# store the 2nd value between -'s as a variable (this is for samples that start with E or R (named like R1234567890-201-D-E03)):
LibNum0_QC = df_FastqList2['sample'].str.split("-").str[1]
#print(LibNum0_QC)
# store the 3rd value btwn -'s as a variable (this is for pos/neg cntrls (named like: NEG20210331-nCoVWGS-201-D)):
LibNum1_QC = df_FastqList2['sample'].str.split("-").str[2]
#print(LibNum1_QC)
# store the 1st value (before 1st -) - the sampleID, as a variable:
SampleID_QC = df_FastqList2['sample'].str.split("-").str[0]
#print(SampleID_QC)
#store the LibNum0_QC length in a variable:
SampleLength_QC = LibNum0_QC.str.len()
#print(SampleLength_QC)
#print(SampleLength_QC.loc[1]) 



#REPLACE SAMPLE string in sample column with SampleID variable (just the CID):
# Will leave controls named as-is (With 'sample' ID) (SampleLength_QC > 5 (=nCoVWGS)), since won't match to metadata
# Will extract CID part of fastqID (SampleID_QC) from non-control samples (SampleLength_QC < 5)

conditions2 = [
    (SampleLength_QC < 5),    
    (SampleLength_QC > 5)
]

choices2 = [SampleID_QC, df_FastqList2['sample']]

df_FastqList2['CID1'] = np.select(conditions2, choices2, default= df_FastqList2['sample'])
#print(df_FastqList2['CID1'])  #correct





#----------Merge FastqID list with CleanMetadata (match to 3 possible sample IDs, then combine CIDs later (with append)):-------------

#Merge sample (CID) with primary CID from plover (in clean metadata)
df_FastqList2_CleanMetadata_merge1 = pd.merge(df_FastqList2, df_CleanMetadata2, how='left', left_on='CID1', right_on='containerid')
#only keep columns needed (matched by sample and CID so keep those, and Ct_combo, and coll date:)
df_FastqList2_CleanMetadata_merge1 = df_FastqList2_CleanMetadata_merge1[['sample','containerid', 'Ct_combo', 'collection_date']]
# Rename the containerid col to CID (and same for other 2 merges below), so can append dfs with same column names (or will add extra columns):
df_FastqList2_CleanMetadata_merge1 = df_FastqList2_CleanMetadata_merge1.rename(columns={'containerid':'CID'})
#print(df_QCsummary_CleanMetadata_merge1) 

#Merge sample (CID) with secondary CID from plover (in clean metadata)
df_FastqList2_CleanMetadata_merge2 = pd.merge(df_FastqList2, df_CleanMetadata2, how='left', left_on='CID1', right_on='second_containerid')
#only keep columns needed (matched by sample and 2ndary CID so keep those, and Ct_combo, and coll date:)
df_FastqList2_CleanMetadata_merge2 = df_FastqList2_CleanMetadata_merge2[['sample', 'second_containerid', 'Ct_combo', 'collection_date']]
# Rename the second_containerid col to CID, so can append dfs with same column names (or will add extra columns):
df_FastqList2_CleanMetadata_merge2 = df_FastqList2_CleanMetadata_merge2.rename(columns={'second_containerid':'CID'})
#print(df_QCsummary_CleanMetadata_merge2) 

#Merge sample (CID) with seq CID from plover (in clean metadata)
df_FastqList2_CleanMetadata_merge3 = pd.merge(df_FastqList2, df_CleanMetadata2, how='left', left_on='CID1', right_on='seq_containerid')
#only keep columns needed (matched by sample and seq CID so keep those, and Ct_combo, and coll date:)
df_FastqList2_CleanMetadata_merge3 = df_FastqList2_CleanMetadata_merge3[['sample', 'seq_containerid', 'Ct_combo', 'collection_date']]
# Rename the seq_containerid col to CID, so can append dfs with same column names (or will add extra columns):
df_FastqList2_CleanMetadata_merge3 = df_FastqList2_CleanMetadata_merge3.rename(columns={'seq_containerid':'CID'})
#print(df_QCsummary_CleanMetadata_merge3) 


   
      #------------Merge the 3 merged dfs together (append)--------

# Append the 3 dfs (add to bottom of last, as samples likely won't match in 3). Then can sort to remove duplicates:
df_FastqList2_CleanMetadata_merge4 = df_FastqList2_CleanMetadata_merge1.append(df_FastqList2_CleanMetadata_merge2, ignore_index=True)

df_FastqList2_CleanMetadata_merge5 = df_FastqList2_CleanMetadata_merge4.append(df_FastqList2_CleanMetadata_merge3, ignore_index=True)
#print(df_QCsummary_CleanMetadata_merge5)

#Drop extra index column:
df_FastqList2_CleanMetadata_merge6 = df_FastqList2_CleanMetadata_merge5.reset_index(drop=True)
#print(df_QCsummary_CleanMetadata_merge6)

# Convert sample column to str, or get error when sort:
df_FastqList2_CleanMetadata_merge6 = df_FastqList2_CleanMetadata_merge6.astype(str)
# Sort to remove duplicate sample IDs, then sort by CID match (so CID matches sort above no match - so when keep first of duplicates, it keeps the matched sample data instead of throwing it out)
df_FastqList2_CleanMetadata_merge7 = df_FastqList2_CleanMetadata_merge6.sort_values(by=['sample', 'CID'], ascending = (True, True))  
#Save this to file to check sorting is correct: (it is, matches are at top)
#df_QCsummary_CleanMetadata_merge7.to_csv("path/to/output/All_Metadata_cleaned_SORTEDwithDUPS_MatchedToSampleIDs.csv", index=False)
# Remove duplicate sampleIDs in merged metadata file:
df_FastqList2_CleanMetadata_merge7 = df_FastqList2_CleanMetadata_merge7.drop_duplicates(['sample'], keep='first') 
#print(df_QCsummary_CleanMetadata_merge7)


# Don't need CID column anymore - can leave as sample
# Only need: 'sample', 'Ct_combo', 'collection_date_y'
df_FastqList2_CleanMetadata_merge8 = df_FastqList2_CleanMetadata_merge7[['sample', 'Ct_combo', 'collection_date']]
#print(df_QCsummary_CleanMetadata_merge8)
# Then change column names to sample, ct, date (as per ncov-tools requirements):
df_FastqList2_CleanMetadata_merge8 = df_FastqList2_CleanMetadata_merge8.rename(columns={'Ct_combo':'ct', 'collection_date':'date'})
#print(df_QCsummary_CleanMetadata_merge9)

#------This df has NaNs - remove before saving final output----
#replace NaNs with NAs for ncov-tools (in date column - need to be handled separately from other columns dealt with below, or it won't get rid of these nan's):
df_FastqList2_CleanMetadata_merge8['date'] = df_FastqList2_CleanMetadata_merge8['date'].replace(np.nan, 'NA')
#print(df_QCsummary_CleanMetadata_merge8)

#nan is actually a string now for some reason, so np.nan won't work. 
#So just replace the string 'nan' with NA for ncov-tools
df_FastqList2_CleanMetadata_merge9 = df_FastqList2_CleanMetadata_merge8.replace('nan', 'NA')
# cells are empty, so replace '' with NA for ncovtools:
df_FastqList2_CleanMetadata_merge10 = df_FastqList2_CleanMetadata_merge9.replace('', 'NA')
# Print final metadata df to screen so user can see it worked (or comment out if not desired - e.g. if automated and no user will see output)
#print(df_FastqList2_CleanMetadata_merge10)



#------------------Save Output:-----------------
#----------Save CleanedMetadata with Matched sampleID/CID column and Cts and Coll Dates----------

# Save as tsv file (save as csv but specify sep=\t to create tsv), and drop index (not needed for ncov-tools metadata.tsv)
# saves in directory you run it in, the MiseqRunID directory:
df_FastqList2_CleanMetadata_merge10.to_csv("metadata.tsv", sep='\t', index=False)


#---Clean up:---------------
#Run Bash Commands to Remove Unnecessary Files: 
bashCommand3 = "rm FastqList.txt" 

subprocess.run(bashCommand3, shell=True, check=True)

