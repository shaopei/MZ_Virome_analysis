# MZ_Virome_analysis
## MZ_Virome_Functional_profiling.bash
Scripts to do the following steps:
+ Overlapping reads were joined using SeqPrep (St John, 2012), and the adaptor sequence was trimmed using Scythe (Buffalo, 2011). Then, TagCleaner (Schmieder et al., 2010) was used to detect and remove WGA2 tag sequences. Next, low-quality bases (Phred quality scores < 30) were trimmed using Sickle (Joshi NA, 2011). 
+ The joined, trimmed reads were then mapped onto Integrated Gene Catalogs (IGC), an integrated catalog of reference genes in the human gut microbiome (Li et al., 2014) by BLASTX using DIAMOND v0.7.5 (Buchfink et al., 2015) with maximum e-value cutoff 0.001, and maximum number of target sequences to report set to 25. DIAMOND outputs were organized by an in-house python script to generate a read hits abundance matrix to import to QIIME 1.9 (Quantitative Insights Into Microbial Ecology) (Caporaso et al., 2010) for further analysis.
+ The matrix was then annotated according to the KEGG annotation of each gene provided by IGC. The annotated abundance matrix was rarified (subsampling without replacement) to 2,000,000 read hits per sample. The KEGG functional profile was then generated by QIIME 1.9 using the command summarize_taxa_through_plots.py.

## MZ_Virome_Denovo_assembly_script.bash
Scripts to do the following steps:
+ The adaptor sequence was trimmed using Scythe (Buffalo, 2011). Then, TagCleaner (Schmieder et al., 2010) was used to detect and remove WGA2 tag sequences. Next, low-quality bases (Phred quality scores < 30) were trimmed using Sickle (Joshi NA, 2011).
+ pair.py python script can identify the trimmed reads that remain paired (forward and reverse).
+ The paired reads were then assembled using Integrated metagenomic assembly pipeline for short reads (InteMAP) (Lai et al., 2015) with insert size 325 bp ± 100 bp. Each sample was assembled separately. 
+ After the first run of assembly, the reads were mapped to the assembled contigs using Bowtie 2 v.2.2.8 (Langmead and Salzberg, 2012) with the following parameter: --local --maxins 800. 
+ The pairs of reads that aligned concordantly at least once were then submitted for the second run of assemble by InteMAP. 
+ Contigs larger than 500 bp from this second assembly were used for subsequent analysis. 
+ Genes were predicted from the assembled contigs that were larger than 500 bp using GeneMarkS v.4.32 (Besemer et al., 2001). 
