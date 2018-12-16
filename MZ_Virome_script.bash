### QC
# use reads that are flagged as high quality (i.e. pass filtered)
#use fastqc to examine the quality of reads
mkdir fastqc
fastqc -o fastqc *.fastq.gz


## short-insert-size library
# Overlapping reads were joined using SeqPrep (St John, 2012), Reads that are merged, the adapters has been removed during this process
sample=H1
mkdir output_SeqPrep
SeqPrep -y G -A AGATCGGAAGAGCACACGTC -B AGATCGGAAGAGCGTCGTGT -f ${sample}_R1.fastq.gz -r ${sample}_R2.fastq.gz -1 ./output_SeqPrep/SeqPreped_${sample}_R1.fastq.gz -2 ./output_SeqPrep/SeqPreped_${sample}_R2.fastq.gz -3 ./output_SeqPrep/SeqPreped_${sample}_discarded_R1.fastq.gz -4 ./output_SeqPrep/SeqPreped_${sample}_discarded_R2.fastq.gz -s ./output_SeqPrep/SeqPreped_${sample}_Merged.fastq.gz -E ./output_SeqPrep/SeqPreped_${sample}_readable_alignment.txt.gz


#Scythe of SeqprepMerge
mkdir SeqPrepMerged_scythed_p0.3
scythe -t -a truseq_adapters.fasta -p 0.3 \
-o SeqPrepMerged_scythed_p0.3/scythed_SeqPrepMerged_${sample}.fastq \
-m SeqPrepMerged_scythed_p0.3/scythed_SeqPrepMerged_${sample}_matches.txt ./output_SeqPrep/SeqPreped_${sample}_Merged.fastq.gz

# TagCleaner (Schmieder et al., 2010) was used to detect and remove WGA2 tag sequences
tagcleaner -verbose -predict -64 -log -fastq ./output_SeqPrep/SeqPreped_${sample}_Merged.fastq.gz

#tagclean trim
echo SeqPreped_P2_Merged.fastq
/home/sc2457/apps/tagcleaner-standalone-0.16/tagcleaner -verbose -log -info -splitall 0 -mm5 4 -mm3 4 -minlen 10 \
-tag5 TGTGTTGGGTGTGTTTGGNNNNNNNNNN -tag3 NNNNNNNNNNCCAAACACACCCAACACA \
-fastq SeqPreped_P2_Merged.fastq -out SeqPrepMerged_tagcleaned/SeqPreped_tagcleaned_P2_Merged

#calcualte the amount of reads and splitted reads 
for f in *_Merged.fastq; do echo $f; cat $f |grep '^@' -c; cat $f |grep '^@' |grep 2$ -c ;done > SeqPrepMerged_tagcleaned_total_vs_splited_reads_count.txt

for f in *_R1.fastq; do echo $f; cat $f |grep '^@' -c; cat $f |grep '^@' |grep 2$ -c ;done > tagcleaned_total_vs_splited_reads_count_R1.txt
for f in *_R2.fastq; do echo $f; cat $f |grep '^@' -c; cat $f |grep '^@' |grep 2$ -c ;done > tagcleaned_total_vs_splited_reads_count_R2.txt



## large-insert-size library
# the adaptor sequence was trimmed using Scythe (Buffalo, 2011)
scythe -t -a truseq_adapters.fasta -p 0.3 -o PF_scythed_p0.3/scythed_${sample}_R1.fastq -m PF_scythed_p0.3/scythed_${sample}_R1_matches.txt ${sample}_R1.fastq.gz
scythe -t -a truseq_adapters.fasta -p 0.3 -o PF_scythed_p0.3/scythed_${sample}_R2.fastq -m PF_scythed_p0.3/scythed_${sample}_R2_matches.txt ${sample}_R2.fastq.gz

# TagCleaner (Schmieder et al., 2010) was used to detect and remove WGA2 tag sequences
# identify the WGA2 tag sequences using tagcleaner -predict
cd PF_scythed_p0.3
mkdir PF_scythed_p0.3_tagcleaned

tagcleaner -verbose -predict -64 -log -fastq scythed_${sample}_R1.fastq
tagcleaner -verbose -predict -64 -log -fastq scythed_${sample}_R2.fastq
# trim the WGA2 tag sequences
tagcleaner -verbose -log -info -splitall 0 -mm5 4 -mm3 4 -minlen 10 -tag5 TGTGTTGGGTGTGTTTGGNNNNNNNNNN -tag3 NNNNNNNNNNCCAAACACACCCAACACA -fastq scythed_${sample}_R1.fastq -out PF_scythed_p0.3_tagcleaned/scythed_tagcleaned_${sample}_R1
tagcleaner -verbose -log -info -splitall 0 -mm5 4 -mm3 4 -minlen 10 -tag5 TGTGTTGGGTGTGTTTGGNNNNNNNNNN -tag3 NNNNNNNNNNCCAAACACACCCAACACA -fastq scythed_${sample}_R2.fastq -out PF_scythed_p0.3_tagcleaned/scythed_tagcleaned_${sample}_R2


#low-quality bases (Phred quality scores < 30) were trimmed using Sickle (Joshi NA, 2011)
cd PF_scythed_p0.3_tagcleaned
mkdir scythed_tagcleaned_sickled
sickle se -q 30 -t sanger -l 30 -f scythed_tagcleaned_${sample}_R1.fastq -o scythed_tagcleaned_sickled/8L_PF_scythed_tagcleaned_sickled_${sample}_R1.fastq
sickle se -q 30 -t sanger -l 30 -f scythed_tagcleaned_${sample}_R2.fastq -o scythed_tagcleaned_sickled/8L_PF_scythed_tagcleaned_sickled_${sample}_R2.fastq


