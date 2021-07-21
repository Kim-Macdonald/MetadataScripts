# MetadataScripts
Scripts for cleaning and merging metadata with artic/ncov-tools results (e.g remove text from Ct fields, merge multiple Ct fields into 1, merge metadata with other results)


## CleanMetadata_CombineCts_MergeWithQCsummary.py

Creates a cleaned metadata.csv file for use with ncov-tools later 

1. cleans text from Ct fields, 
2. combines Cts in multiple columns into 1, 
3. formats dates - so are left with 3 columns for metadata.tsv: sample, ct, date) 
4. then merges with Result Summary file from mergeQCresults_plusMissing (https://github.com/Kim-Macdonald/mergeQCresults_plusMissing), VoCcaller (https://github.com/Kim-Macdonald/VoCcaller), LineageUpdater (https://github.com/Kim-Macdonald/LineageUpdater), to produce a file of metadata for every sample in result summary file. (this can be used to fill in metadata for all runs to-date, AFTER ncov-tools analysis has completed).


## CleanMetadata_CreateMetadataTsvFileForRun_v0.py

Creates metadata.tsv for ncov-tools for each sequencing run in the directory it's run in, by doing the following:

  1. cleans text from Ct fields, 
  2. combines Cts in multiple columns into 1
  3. matches fastq sample IDs to 3 potential ID fields in the metadata file (container id, secondary container id, sequencing id)
  4. formats dates - so are left with 3 columns for metadata.tsv: sample, ct, date. 
  5. Then merges with a list of Fastq files for the run to match metadata for only those samples in the directory being analyzed, and outputs the metadata.tsv file for ncov-tools (with sample, ct, date columns).




