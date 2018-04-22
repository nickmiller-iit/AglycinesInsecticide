# Notes on analysis of soybean aphid low-dose insecticide RNA-Seq experiment.

## General overview

Analysis steps will be recorded in a Makefile for repeatability. Results and observations will be recorded in these notes. Git will be used to track changes in the Makefile and these notes. Intermediate output and data files will not generally be included in the git repo because the files are usually very large. We may add the raw data to the repo later so the repo includes everything needed to repeat the analysis.

## Initial QC

To see how our raw data looks, ran fastqc (version 0.11.5) on the raw reads.

In general, quality scores look pretty good, R2 ("reverse") reads tend to be slightly lower quality. Sample C12, rep2 was rather lower quality, but probably not fatal after quality trimming. To a lesser extend Sample T12, rep1 was also slightly lower quality. In both cases, the low quality was observed in both lanes of sequencing, suggesting the issue is with the library.

Looking at overrepresented sequences, most R1 ("forward") reads had hits to "TruSeq Adapter, Index i", where i is in {2, 3, 5, 5, 6, 7, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 19, 20, 21, 22, 25, 27}. Most R2 reads had hits to "Illumina Single End PCR primer 1". Hits ranged from about 97% match over 36 nt to 100% match over 50nt.

## Concatenating reads

Before trimming and further QC, concatenated reads from the sample sample and read direction together. This just reduces the number of input files we have to deal with. At the same time, decompress the files because not all programs handle gzipped fastq. After decompressing and concatenating, used the perl script pe-sync-2-files.pl, written by John Garbe to check the forward and reverse reads are still synced up.

## Trimming reads and QC

Trimming was done using Trimmomtic v0.36. Parameters are in the Makefile. Post trimming QC showed that

* Retained paired reads were between 6.6 and 13 million paired reads per sample. The only exception was sample C12, rep2 with 4.2 million paired reads. This is probably to be expected, given this library gave lower quality overall.
* Unpaired reads were generally 5-10% the number of paired reads for read 1 and 0.5-1 % for read 2.
* Unpaired reads generally hade more variation in quality scores than paired reads. This probably makes sense given that unpaired reads had a pair that did not pass trimming.
* No adaptor or primer sequences were detected in any set of paired reads.
* Low levels of adaptor/primer sequences were detected in a few sets of unpaired reads
  - C12, rep 2: read 1 (1.7%)
  - C12, rep 3: read1 (0.8%)
  - C12, rep 4: read 1 (0.12%)
  - E12, rep 3: read 1 (0.15%)
  - T12, rep 1: read 1 (0.3%), read 2(0.57%)

Overall, this looks pretty respectable. We are probably OK with the remaining low levels of adaptor contamination. We are not doing de-novo assembly so, at worst, those reads will fail to map.

## Mapping reads


Reads were mapped to the genome using GSNAP v2017-04-24 with default parameters, except that -N 1 was set to enable splicing. Resulting BAM files were sorted and idexed with samtools v1.4.

GSNAP will only map paired reads or single reads, not a mixture of both. We could work around this by doing paired and single reads separatedly and then merging the output. Chose not to do this right away and just align paired reads as in most cases we will only gain a few % extra reads and the quality of the single reads is lower. 

## Updating annotations

Initial inspection of read alignments suggested thare are a number of points in the genome where appreciable numbers of reads are aligned, but there is not a reference annotation. Used StringTie v 1.3.3b to assemble transcripts from each sample based on aligned reads. Because StringTie was supplied with the reference annotation, it kept track of the assembled exons & transcripts that overlap with existing reference annotations. Stringtie was then used to merge assemblies from individual samples into a single GFF. A quick analysis with gffcompare v0.9.8 shows that the merged assembly contains 33032 exons and 11509 loci that are not in the original reference annotation.

After attempting to count reads, it became apparent that there were a number of genes identified by stringtie that overlapped with genes from the reference annotation. This causes problems because the assignment of a read to a gene becomes ambiguous when it maps to >1 gene annotation. The cautious approach is to only retain those genes identified by stingtie that do not overlapp with genes from the reference annotation. Used gffread v0.9.8d to convert the refrerence annotations from gff3 to gtf. Then used bedtools v2.26.0 to isolate features from the merged StringTie assembly that do not overlap at all with features in the reference annotation. Finally combined the two to produce a GTF file that contains orginal reference annotations plus novel genes found by StringTie.

## Read counting

Used StringTie to generate read counts. Set the output to be compatible with the Ballgown R package. This is by the same group as StringTie, and has its own approach to differential expression analysis. However (according to the docs) it can easily format data for DESeq2 and EdgeR, so it seems as good a way to go as any.

This turns out not to be acurate. Ballgown will only give expression levels as FPKM, *except* for exon-level data. Since we plan to use DESeq2, we need read counts at gene level. Fortunately featureCounts v1.5.2 can give us exactly what we want. It counts fragments at the exon level and then groups by gene id, giving us the fragment counts per gene per sample.

## Updating to new assembly

As of Feb 2018, there is an updated and significantly improved assembly, labeled version 6.0. This is the one we should use for genome-related publications. Luckily, because everything was set up using a Makefile, it should be easy to simply rerun everything described above, starting with read mapping.

## Doubts about stringtie

Early attemts to test for differential gene expression suggest something fishy with the combined GTF. We get ridiculous numbers of genes lighting up as significantly differentially expressed. Most of these appear to have come from the stringtie assebly. Ballgown aslo reveals that the stringtie genes often have many putative isoforms. The maker gene set genes do not, because we just get the maker gene model.

## Trying kalliso

Kallisto and other "pseudoalignment" transcript quantifiers are the next big thing in RNA-Seq analysis. In addition to being much faster, kallisto explicity works with transcripts (including isoforms). This is a potentially a good thing because it should get us around the longstanding problem of what to do with muliple isoforms in a reference.

### With maker gene set

Kallisto works with a fasta file of transcript sequences. We already have this for the maker gene set, so started with makefile targets to index that and quantify transcripts.

### With modified maker gene set

Three tanscripts appeared to be causing a few of the 12-hour samples to behave oddly in PCA plots. Removed the transcripts with `seqkit grep` and re-ran the analysis.

## BUSCO analysis

Initial playing around with kalliso and sleuth, plus some eyballing of the gene models in Apollo makes me suspect that the origiginal gene set is fragmented. RAN a BUSCO analysis on the OGS protein sequences.

Surprised to find the rate of fragmented BUSCOs is pretty low!

	C:90.9%[S:87.0%,D:3.9%],F:2.7%,M:6.4%,n:1658

	1507	Complete BUSCOs (C)
	1443	Complete and single-copy BUSCOs (S)
	64	Complete and duplicated BUSCOs (D)
	44	Fragmented BUSCOs (F)
	107	Missing BUSCOs (M)
	1658	Total BUSCO groups searched

