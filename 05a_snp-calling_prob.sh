#This script processes metagenome data to reconstruct phylogeny of symbionts and mitochondria using a probabilistic SNP-identification approach.
cp ref/Gamma*fasta ref/Delta*.fasta ref/Spiro.fasta ref/Mito.fasta .
#copy input files and reference files in the working directory

for sample in $(cut -f1 sample_info.list); do
	echo "$sample"
	cp $DATAPATH/${sample}_R?.fq.gz . 

	echo "  quality filtering of input files by trimmomatic v0.36..."
	unpigz ${sample}_R?.fq.gz
	java -jar $TOOLPATH/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 2 -phred33\
	 ${sample}_R1.fq ${sample}_R2.fq\
	 ${sample}_pre1.fastq ${sample}_1_unpaired.fastq\
	 ${sample}_pre2.fastq ${sample}_2_unpaired.fastq\
	 LEADING:10 SLIDINGWINDOW:4:2 MINLEN:50
	  #first 10-bases off, window of 4 bases to check quality >2, shorter than 50 to dump.
	rm  ${sample}_R?.fq.gz *unpaired.fastq

	echo "  splitting symbiont reads into different symbiont species (bbsplit.sh from bbmap v36.86)..."
	bbsplit.sh in1=${sample}_pre1.fastq in2=${sample}_pre2.fastq\
	 ref=Delta1a.fasta,Delta1b.fasta,Delta3.fasta,Delta4.fasta,Gamma1.fasta,Gamma3.fasta,Spiro.fasta\
	 basename=${sample}_%_#.fastq ambiguous2=toss refstats=${sample}_bbsplit.refstats threads=2 overwrite
	
	echo "  bbmap mapping clean reads to a reference mtDNA with 95% identity filter (for potential NUMTs removal)..."
	bbmap.sh in1=${sample}_pre1.fastq in2=${sample}_pre2.fastq ref=Mito.fasta\
	 out1=${sample}_mito_1.fastq out2=${sample}_mito_2.fastq\
	 minid=0.95 killbadpairs pairedonly mappedonly maxindel=1000 pairlen=1000\
	 nodisk overwrite threads=2
	rm ${sample}_pre?.fastq

	echo "  adding information-tags and filtering of mapped files for all symbionts and mitochondria refs"
	#Extract relevant values from <sample_infp.list> for sample ID, library ID, and platform unit.
	#<sample_info.list> contains the following data without header (tab-delimited)
	 # SampleID	LibraryID	PlatformUnit	COI-Haplotype	Location_ID	Location	ReadLength	SampledYear
	 #Example
	 # OalgCAVL_A11	3561.BB	MPIPZ	A	CAVL	Cavoli	150	2016
	 # OalgCAVL_A12	3561.BJ	MPIPZ	A	CAVL	Cavoli	150	2016
	 # OalgCAVL_A13	3561.BM	MPIPZ	A	CAVL	Cavoli	150	2016
	LIB_ID=$(grep -w $sample sample_info.list | cut -f 2)
	PU=$(grep -w $sample sample_info.list | cut -f 3)
	for ref in Mito Gamma1 Gamma3 Delta1a Delta1b Delta3 Delta4 Spiro; do
		bbmap.sh in1=${sample}_${ref}_1.fastq in2=${sample}_${ref}_2.fastq ref=$ref.fasta\
		outm=${sample}_$ref.sam local killbadpairs pairedonly pairlen=2000\
		rgid=$sample rgsm=$sample rglb=$LIB_ID rgpu=$PU rgpl=ILLUMINA\
		threads=2 nodisk ow covstats=${sample}_$ref.covstats
		 #'local' to softclip ugly ends of alignment.
		 #'killbadpairs' to mark unmapped if reads are discordantly mapped or in a wrong insert size.
		 #'pairedonly' to keep only paired reads.
		samtools view -bS -F 4 -@ 2 ${sample}_${ref}.sam > ${sample}_${ref}.bam
		 #converting sam to bam file for storage (samtools v1.3.1)
		 # '-b' to BAM-output, '-S' to auto detect input-format, '-F 4' to exclude any unmapped reads.
		rm ${sample}_${ref}.sam
		samtools view -h -q 20 -@ 2 ${sample}_${ref}.bam | grep -v XS:i\
		 | samtools view -buS -@ 2 - | samtools sort -@ 2 - > ${sample}_${ref}_MinQ20_sorted.bam
		 # Filter bam files to remove poorly mapped reads (non-unique mappings (XS:i) and mapping quality score < 20)
		 #'-h to include headers -bu to output uncompressed BAM, -S to auto detect input format'. BAM sorted.
		bam clipOverlap --in ${sample}_${ref}_MinQ20_sorted.bam\
		  --out ${sample}_${ref}_MinQ20_clipped.bam  --storeOrig OC --stats --unmapped
		 #clipping overlapped read pairs (bam from bamUtil v1.0.4)
		samtools view -hb -F 4 -@ 2 ${sample}_${ref}_MinQ20_clipped.bam\
		 > ${sample}_${ref}_MinQ20_clipped_mapped.bam
		 #Remove the reads that became unmapped in the clipping		 
		java -Xmx60g -jar $TOOLPATH/picard.jar MarkDuplicates I=${sample}_${ref}_MinQ20_clipped_mapped.bam\
		 O=${sample}_${ref}_MinQ20_clipped_mapped_dedup.bam M=${sample}_${ref}_MinQ20_clipped_mapped_dupstat.txt\
		 VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true
		 #Deduplication with picard v2.9.2 
		samtools sort -o ${sample}_to_${ref}_MinQ20_dedup_sorted.bam ${sample}_${ref}_MinQ20_clipped_mapped_dedup.bam -@ 2
		samtools index -b ${sample}_to_${ref}_MinQ20_dedup_sorted.bam
		 #resulting bam file sorting and indexing
		echo "    done with $ref"
	done; echo
done; ehco

for reffile in Mito Gamma1 Gamma3 Delta1a Delta1b Delta4 Spiro; do
	echo "$reffile is being processed..."
	echo "  indexing reference file and make dictionary file"
	samtools faidx $reffile.fasta
	java -Xmx60g -jar $TOOLPATH/picard.jar CreateSequenceDictionary \
	 R=${reffile}.fasta O=${reffile}.dict
	
	echo "  gatk identifying target intervals for indel realignment (GATK v3.7)"
	ls *_to_${reffile}_MinQ20_dedup_sorted.bam > Samples_To_${reffile}_dedup_sorted_bam.list
	java -Xmx40g -jar $TOOLPATH/GenomeAnalysisTK.jar \
	  -T RealignerTargetCreator \
	  -R $reffile.fasta \
	  -I Samples_To_${reffile}_dedup_sorted_bam.list \
	  -o ${reffile}_indel.intervals \
	  -rf BadMate \
	  -nt 24
	  
	echo "  gatk re-aligning around indels"
	java -Xmx40g -jar $TOOLPATH/tools/GenomeAnalysisTK.jar \
	   -T IndelRealigner \
	   -R $reffile.fasta \
	   -I Samples_To_${reffile}_dedup_sorted_bam.list \
	   -targetIntervals ${reffile}_indel.intervals \
	   --consensusDeterminationModel USE_READS \
	   --nWayOut _realigned.bam
	mkdir -p mapped/$reffile; mv *_to_${reffile}_MinQ20_dedup_sorted_realigned.ba? mapped/$reffile
	  #Realigned sorted BAM files and their index files (bai) are stored in mapped/$reffile
	rm ${sample}_to_${ref}_MinQ20_dedup_sorted.bam*
	
	echo "  ANGSD step1 for depth-count for the next SNP-identification step (v0.929)"
	ls mapped/$reffile/*realigned.bam > bam_list.txt
	n_ind=$(cat bam_list.txt | wc -l)          #Calculate the number of samples
	samtools faidx $reffile.fasta
	angsd -P 24 -b bam_list.txt -out ${n_ind}Samples_to_${reffile}_MinQ20_qc -ref $reffile.fasta\
	 -trim 0 -remove_bads 0 -uniqueOnly 0 -only_proper_pairs 1 -doQsDist 1 -doDepth 1 -doCounts 1 -maxDepth 24000
	  #Remove_bads 0 and -uniqueOnly 0 as they were already filtered from the bam files.
	  #Distribution of read coverage is stored in <${n_ind}Samples_to_${reffile}_MinQ20_qc.info>.
	
	echo "  ANGSD step2 for processing for haploids with inbreeding adjustment"
	mindepth=$(sed -e '1,/Global_depth/d' ${n_ind}Samples_to_${reffile}_MinQ20_qc.info| awk '$2>0.05 {print $1}'| head -n1)
	maxdepth=$(sed -e '1,/Global_depth/d' ${n_ind}Samples_to_${reffile}_MinQ20_qc.info| awk '$2<0.95 {print $1}'| tail -n1)
	  #Obtained the depth range (5 ~ 95 %tile) from the previous angsd output to pre-reduce the noise in SNP-identification.
	for i in $(seq 1 $n_ind); do echo "1.00000" >> all_1.indF; done
	  #indF; individual inbreeding coefficients are all set as 1 for haploids (i.e. all are set as 'homozygote').
	snp_pval=0.01
	prior=2 #For haploid, uniform prior is assumed.
	angsd -P 24 -b bam_list.txt -ref $reffile.fasta -anc $reffile.fasta\
	 -out ${n_ind}Samples_to_${reffile}_MinQ20_SNPpval${snp_pval}_Prior${prior}_haploid -doGeno 8\
	 -GL 2 -doMaf 1 -doMajorMinor 1 -doPost $prior -doCounts 1 -underFlowProtect 1\
	 -remove_bads 0 -uniqueOnly 0 -only_proper_pairs 1 -SNP_pval $snp_pval\
	 -setMinDepth $mindepth -setMaxDepth $maxdepth -minInd $n_ind -minQ 2 -minMaf 0.01\
	 -doQsDist 1 -doDepth 1 -maxDepth 8000 -doSaf 2 -indF all_1.indF
	  #-doGeno 8 for ngsDist
	  #-doSaf 2 -indF to incorporate indF coefficient
	  #-minInd = $n_ind (i.e. all samples)
	n_sites=$(zcat ${n_ind}Samples_to_${reffile}_MinQ20_SNPpval${snp_pval}_Prior${prior}_haploid.mafs.gz | tail -n +2 | wc -l)
	  #get the number of SNP sites

	if [ $n_sites != 0 ]; then	#Escape if no SNP-sites were identified.
		echo "  Calculating genetic distances based on genotype probabilities (ngsDist v1.0.2) and phylogenies (FastME v2.1.5.1)"
		cut -f1 sample_info.list | sort > sample.list
		bblock=50; if (( $n_sites <= $bblock )); then bblock=$(echo $n_sites | awk '{print int($1/2)}'); fi 
		  #This is to avoid crashing in the raxml step later.
		ngsDist -geno ${n_ind}Samples_to_${reffile}_MinQ20_SNPpval${snp_pval}_Prior${prior}_haploid.geno.gz\
		   -probs -n_ind $n_ind -n_sites $n_sites -labels sample.list -verbose 1 -n_threads 24\
		   -o ${n_ind}Samples_to_${reffile}_MinQ20_SNPpval${snp_pval}_Prior${prior}_haploid_boot.dist\
		   -n_boot_rep 100 -boot_block_size $bblock -seed 2
		sed 's/-0.0000000000/0.0000000000/g' ${n_ind}Samples_to_${reffile}_MinQ20_SNPpval${snp_pval}_Prior${prior}_haploid_boot.dist > tmp.dist
		mv tmp.dist ${n_ind}Samples_to_${reffile}_MinQ20_SNPpval${snp_pval}_Prior${prior}_haploid_boot.dist
		  #This is to fix an issue that can be caused by '-0.0000000000' in the dist matrix being not recognized as numeric.
		fastme -T 24          -i ${n_ind}Samples_to_${reffile}_MinQ20_SNPpval${snp_pval}_Prior${prior}_haploid_boot.dist\
		-D 101 -m B -n B      -o ${n_ind}Samples_to_${reffile}_MinQ20_SNPpval${snp_pval}_Prior${prior}_haploid_boot.trees
		echo "  Bootstrap-supported tree building with raxml v8.2.11"
		head -n 1                ${n_ind}Samples_to_${reffile}_MinQ20_SNPpval${snp_pval}_Prior${prior}_haploid_boot.trees > main.nwk
		tail -n+2                ${n_ind}Samples_to_${reffile}_MinQ20_SNPpval${snp_pval}_Prior${prior}_haploid_boot.trees > boot.nwk
		raxmlHPC -m GTRCAT    -n ${n_ind}Samples_to_${reffile}_MinQ20_SNPpval${snp_pval}_Prior${prior}_haploid_bootstrap.tree\
		   -t main.nwk -z boot.nwk -f b
		mv    RAxML_bipartitions.${n_ind}Samples_to_${reffile}_MinQ20_SNPpval${snp_pval}_Prior${prior}_haploid_bootstrap.tree\
								 ${n_ind}Samples_to_${reffile}_MinQ20_SNPpval${snp_pval}_Prior${prior}_haploid_bootstrap.tree
	else
		echo "  No SNPs were found. Use <06a_snp-calling_prob_repeat.sh> with filtering samples."
	fi
done