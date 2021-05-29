# MetadataScripts
Scripts for cleaning and merging metadata with artic/ncov-tools results (e.g remove text from Ct fields, merge multiple Ct fields into 1, merge metadata with other results)


## CleanMetadata_CombineCts_MergeWithQCsummary.py

Creates metadata.tsv for ncov-tools (cleans text from Ct fields, combines Cts in multiple columns into 1, formats dates - so are left with 3 columns for metadata.tsv: sample, ct, date) then merges with Result Summary file from mergeQCresults_plusMissing (https://github.com/Kim-Macdonald/mergeQCresults_plusMissing), VoCcaller (https://github.com/Kim-Macdonald/VoCcaller), LineageUpdater (https://github.com/Kim-Macdonald/LineageUpdater)).

