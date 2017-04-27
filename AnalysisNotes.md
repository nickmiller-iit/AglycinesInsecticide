# Notes on analysis of soybean aphid low-dose insecticide RNA-Seq experiment.

## General overview

Analysis steps will be recorded in a Makefile for repeatability. Results and observations will be recorded in these notes. Git will be used to track changes in the Makefile and these notes. Intermediate output and data files will not generally be included in the git repo because the files are usually very large. We may add the raw data to the repo later so the repo includes everything needed to repeat the analysis.

## Initial QC

To see how our raw data looks, ran fastqc (version 0.11.5) on the raw reads.

In general, quality scores look pretty good, R2 ("reverse") reads tend to be slightly lower quality. Sample C12, rep2 was rather lower quality, but probably not fatal after quality trimming. To a lesser extend Sample T12, rep1 was also slightly lower quality. In both cases, the low quality was observed in both lanes of sequencing, suggesting the issue is with the library.

Looking at overrepresented sequences, most R1 ("forward") reads had hits to "TruSeq Adapter, Index i", where i is in {2, 3, 5, 5, 6, 7, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 19, 20, 21, 22, 25, 27}. Most R2 reads had hits to "Illumina Single End PCR primer 1". Hits ranged from about 97% match over 36 nt to 100% match over 50nt.
