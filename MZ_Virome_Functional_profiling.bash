### QC
# use reads that are flagged as high quality (i.e. pass filtered)
#use fastqc to examine the quality of reads
mkdir fastqc
fastqc -o fastqc *.fastq.gz


## short-insert-size library
# Overlapping reads were joined using SeqPrep (St John, 2012). 
# the adapters has been removed from the Merged reads.
sample=H1
mkdir output_SeqPrep
SeqPrep -y G -A AGATCGGAAGAGCACACGTC -B AGATCGGAAGAGCGTCGTGT -f ${sample}_R1.fastq.gz -r ${sample}_R2.fastq.gz \
-1 ./output_SeqPrep/SeqPreped_${sample}_R1.fastq.gz -2 ./output_SeqPrep/SeqPreped_${sample}_R2.fastq.gz \
-3 ./output_SeqPrep/SeqPreped_${sample}_discarded_R1.fastq.gz -4 ./output_SeqPrep/SeqPreped_${sample}_discarded_R2.fastq.gz \
-s ./output_SeqPrep/SeqPreped_${sample}_Merged.fastq.gz -E ./output_SeqPrep/SeqPreped_${sample}_readable_alignment.txt.gz


# the adapters has been removed from the Merged reads, so there is no need to run scythe.
# TagCleaner (Schmieder et al., 2010) was used to detect and remove WGA2 tag sequences
# identify the WGA2 tag sequences using tagcleaner -predict
cd output_SeqPrep
mkdir SeqPrepMerged_tagcleaned
gzip -d SeqPreped_${sample}_Merged.fastq.gz

tagcleaner -verbose -predict -64 -log -fastq SeqPreped_${sample}_Merged.fastq

# trim the WGA2 tag sequences
tagcleaner -verbose -log -info -splitall 0 -mm5 4 -mm3 4 -minlen 10 -tag5 TGTGTTGGGTGTGTTTGGNNNNNNNNNN -tag3 NNNNNNNNNNCCAAACACACCCAACACA -fastq SeqPreped_${sample}_Merged.fastq -out SeqPrepMerged_tagcleaned/SeqPreped_tagcleaned_${sample}_Merged.fastq

# low-quality bases (Phred quality scores < 30) were trimmed using Sickle (Joshi NA, 2011)
cd SeqPrepMerged_tagcleaned
mkdir SeqPrepMerged_tagcleaned_sickled
sickle se -q 30 -t sanger -l 30 -f SeqPreped_tagcleaned_${sample}_Merged.fastq -o SeqPrepMerged_tagcleaned_sickled/5S_PF_SeqPrepMerged_tagcleaned_sickled_${sample}.fastq

# The joined, trimmed reads were then mapped onto Integrated Gene Catalogs (IGC), an integrated catalog of reference genes in the human gut microbiome (Li et al., 2014)
# by BLASTX using DIAMOND v0.7.5 (Buchfink et al., 2015) with maximum e-value cutoff 0.001, and maximum number of target sequences to report set to 25. 
cd SeqPrepMerged_tagcleaned_sickled
mkdir diamond_output
mkdir diamond.temp
# make database from the protein sequences from IGC where diamond can blastx to
/programs/diamond-0.7.5/diamond makedb --in IGC.pep -d IGC.pep -p 16
# BLASTX using DIAMOND v0.7.5 
/programs/diamond-0.7.5/diamond blastx -v -d IGC.pep -q 5S_PF_SeqPrepMerged_tagcleaned_sickled_${sample}.fastq -a diamond_output/5S_PF_SeqPrepMerged_tagcleaned_sickled_${sample}_matches -t diamond.temp 

cd diamond_output
/programs/diamond-0.7.5/diamond view -a 5S_PF_SeqPrepMerged_tagcleaned_sickled_${sample}_matches.daa -o 5S_PF_SeqPrepMerged_tagcleaned_sickled_${sample}_matches.m8

# DIAMOND outputs were organized by an in-house python script m8_parser.py to generate a read hits abundance matrix (diamond_blastx_m8_parsed_subject.txt)
python m8_parser.py





