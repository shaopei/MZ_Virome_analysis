### QC
# use reads that are flagged as high quality (i.e. pass filtered)
#use fastqc to examine the quality of reads
mkdir fastqc
fastqc -o fastqc *.fastq.gz

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

#use in house python script to identify pairs of reads
cd scythed_tagcleaned_sickled
mkdir scythed_tagcleaned_sickled_paired
for f in *_R1.fastq; do python pair.py $f; done

# De novo assembled using Integrated metagenomic assembly pipeline for short reads (InteMAP) (Lai et al., 2015) with insert size 325 bp ± 100 bp. 
cd scythed_tagcleaned_sickled_paired
# The first run of assembly. Each sample was assembled separately.
mkdir ${sample}
mv 8L_PF_scythed_tagcleaned_sickled_${sample}_paired_R1.fastq 8L_PF_scythed_tagcleaned_sickled_${sample}_paired_R2.fastq ${sample}/.
cd ${sample}
echo 8L_PF_scythed_tagcleaned_sickled_${sample}_paired_R1.fastq 8L_PF_scythed_tagcleaned_sickled_${sample}_paired_R2.fastq >  ${sample}_8L_files
echo -libraryname y1 -insertsize 325 100 -type sanger -technology illumina-long > ${sample}_library_info_file
python InteMAP_v1.0/runInteMAP.py ${sample}_8L_files ${sample}_library_info_file

# After the first run of assembly, the reads were mapped to the assembled contigs using Bowtie 2 v.2.2.8 (Langmead and Salzberg, 2012) with the following parameter: --local --maxins 800. 
mkdir ${sample}_8Lpp_bowtie2_local
mkdir ${sample}_unpp_bowtie2_local
echo ${sample} > ${sample}_8Lraw_mapto_8Lraw.out.fa_bowtie2.log
bowtie2-build out.fa out.fa.bowtie
bowtie2 --local --maxins 800 -p 36 -x out.fa.bowtie -1 8L_PF_scythed_tagcleaned_sickled_${sample}_paired_R1.fastq -2 8L_PF_scythed_tagcleaned_sickled_${sample}_paired_R2.fastq -S out.bowtie.sam_8L_800 --al-conc ${sample}_8Lpp_bowtie2_local/${sample}_pp_800_R%.fastq --un-conc ${sample}_unpp_bowtie2_local/${sample}_8L_800_unpp_R%.fastq >> ${sample}_8Lraw_mapto_8Lraw.out.fa_bowtie2.log 2>&1

#The pairs of reads that aligned concordantly at least once were then submitted for the second run of assemble by InteMAP. Each sample was assembled separately.
cd ${sample}_8Lpp_bowtie2_local
echo ${sample}_pp_800_R1.fastq ${sample}_pp_800_R2.fastq > ${sample}_pp_files
python InteMAP_v1.0/runInteMAP.py ${sample}_pp_files ../${sample}_library_info_file > ${sample}_8Lpp_InteMAP.log 2>&1

# soft link all contigs from the second assembly to a folder InterMAP_out.fa_pp and combine them into a file
cd ../../
mkdir InterMAP_out.fa_pp
ln -s ${sample}/${sample}_8Lpp_bowtie2_local/out.fa InterMAP_out.fa_pp/${sample}_8Lpp_out.fa
cd InterMAP_out.fa_pp
python add_prefix_and_Combine_fasta.py 
#output is Combined_InteMAP_pp_contigs.fa

#Contigs larger than 500 bp from this second assembly were used for subsequent analysis.
python filter_contig_by_size.py Combined_InteMAP_pp_contigs.fa 500 
# output is Combined_InteMAP_pp_contigs.fa_LargerThan500

#Genes were predicted from the assembled contigs that were larger than 500 bp using GeneMarkS v.4.32 (Besemer et al., 2001). 
gmsn.pl --faa --pdf --phage Combined_InteMAP_pp_contigs.fa_LargerThan500
# output is Combined_InteMAP_pp_contigs.fa_LargerThan500.faa
# substitue tab for _
sed -i $'s/\t/_/g'  Combined_InteMAP_pp_contigs.fa_LargerThan500.faa

#The amino acid sequence of the predicted genes were then used in a BLASTp search against the customized viral NR database of 621,095 sequences downloaded from http://pathology.wustl.edu/virusseeker/data/VirusDBNR_20131107_ID98.tgz (Lim et al., 2015) using DIAMOND v.0.7.5 (Buchfink et al., 2015) with maximum e-value cutoff 0.001 and maximum number of target sequences to report set to 25. 
mkdir diamond_output
mkdir diamond.temp
/programs/diamond-0.7.5/diamond makedb --in VirusDBNR_20131107_ID98.fa -d VirusDBNR_20131107_ID98 -p 16
/programs/diamond-0.7.5/diamond blastp -p 8 -v -d VirusDBNR_20131107_ID98 -q Combined_InteMAP_pp_contigs.fa_LargerThan500.faa -a diamond_output/Combined_InteMAP_pp_LargerThan500_Lim_k25_matches -t diamond.temp
#output is Combined_InteMAP_pp_LargerThan500_Lim_k25_matches.daa
cd diamond_output
/programs/diamond-0.7.5/diamond view -a Combined_InteMAP_pp_LargerThan500_Lim_k25_matches.daa -o Combined_InteMAP_pp_LargerThan500_Lim_k25_matches.m8

# Using the BLASTp results (Combined_InteMAP_pp_LargerThan500_Lim_k25_matches.m8) as input, 
# the taxonomy of each gene was assigned by the lowest-common-ancestor algorithm in MEtaGenome ANalyzer (MEGAN5) v.5.11.3 (Huson et al., 2011) 
# with the following parameters: Min Support: 1, Min Score: 40.0, Max Expected: 0.01, Top Percent: 10.0, Min-Complexity filter: 0.44. 
# The taxonomy of the contigs was then assigned using a voting system: the taxonomy of each contig was determined by the majority of taxonomic assignments of the genes on the contig.


