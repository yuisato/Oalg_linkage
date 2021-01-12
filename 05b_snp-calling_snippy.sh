#This scripts process metagenome data to infer SNP-phylogeny of symbionts and mitochondria using a deterministic genotyping approach.

for sample in $(cut -f1 sample_info.list); do
	echo "$sample"
	cp $DATAPATH/${sample}_R?.fq.gz . 
	cp ref/Gamma*fasta ref/Delta*.fasta ref/Spiro.fasta ref/Mito.fasta .
	 #copy input files and reference files in the working directory

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
	
	echo "  SNP-calling with a deterministic genotyping (snippy v3.2)"
	for ref in Mito Gamma1 Gamma3 Delta1a Delta1b Delta4 Spiro; do
		echo "    $ref"
		snippy --cpus 2  --outdir ./snp_${ref}/${ref}_${sample} --ref $ref.fasta \
		  --pe1 ${sample}_${ref}_1.fastq --pe2 ${sample}_${ref}_2.fastq --mincov 5
	done
done; ehco

for ref in Mito Gamma1 Gamma3 Delta1a Delta1b Delta4 Spiro; do
	echo "  Identifying core SNPs across all samples with 'snippy-core'"
	cd snp_${ref}
	snippy-core --prefix=core_${ref}_Oalg_all ${ref}_Oalg* --noref
	echo "  IQTree-building of the coreSNP tree (v1.5.5)"
	iqtree -nt 4 -alrt 1000 -bb 1000 -s core_${ref}_Oalg_all.aln -m GTR+ASC -redo
	  #GTR+ASC model for core SNP trees
	cd ..
done
paste snp_*/core_*_Oalg_all.txt > all_horizontal_coverage_snippy.txt 
  #This contains horizontal read coverage (read depths >= 5X) for all samples for all references.

echo "If no/too-few SNPs were found, use <06b_snp-calling_snippy_repeat.sh> with filtering samples."