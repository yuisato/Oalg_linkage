#This script prepares single-copy gene reference sequences from symbiont genome bins using checkM and
# quantifies symbiont composition using Kallisto.

echo "checkM (identifying single-copy genes; v1.0.7)"
mkdir bins; cp ref/Gamma?.fasta ref/Delta*.fasta ref/Spiro.fasta /bins
#symbiont bins <ref/*.fasta> are copied in bins/ folder
checkm lineage_wf -t 4 -f checkm -x fasta --nt bins checkm
checkm bin_qa_plot --image_type pdf -x fasta checkm bins checkm

echo "extracting single marker genes from output of CheckM"
#Sequences of predicted protein-coding genes are in e.g. checkm/bins/Delta1a/genes.fna
#Sequences that matched with single-copy marker genes are included in e.g. checkm/storage/marker_gene_stats.tsv
#Sequences that were detected as contaminations are listed in e.g. checkm/storage/aai_qa/Delta1a/*.masked.faa
for bin in $(ls checkm/bins); do
	grep $bin checkm/storage/marker_gene_stats.tsv | sed "s/NODE/\nNODE/g" | sed "s/': {'/\n/g" | grep NODE > tmp1
	  #All sequence headers in bins need to begin with 'NODE'. This includes headers of all marker genes.
 	cat checkm/storage/aai_qa/$bin/*.masked.faa | sed "s/NODE/\nNODE/g" | grep NODE > tmp2
	  #These are headers of sequences detected as contaminations. This is to be removed.
	filterbyname.sh in=checkm/bins/$bin/genes.fna names=tmp1 include out=stdout.fasta |\
	  filterbyname.sh in=stdin.fasta names=tmp2 include=f out=stdout.fasta |\
	  sed "s/NODE/${bin}_NODE/g" > checkm/$bin.scg.fasta
	  #Filtering protein-coding gene sequences using the lists above to compile single-copy gene sequences per bin.
	  #<filterbyname.sh> is part of BBmap program v36.86..
	  #Headers of the single-copy gene sequences are modified to include the bin name (e.g. Delta1a_NODE*).
done; rm tmp?
cat checkm/*.scg.fasta > Oalg_symbionts.scg.fasta

#echo "cd-hit (checking cross-symbiont-species sequence similarity at 90%; v4.5.4)"
#Remove undifferentiated genes by 90% clustering from the final single-copy gene sequences.
cd-hit-est -i Oalg_symbionts.scg.fasta -o Oalg_symbionts.scg.cdhit90.fasta -c 0.90 -n 10 -d 0
#Clustered sequences were listed to remove from Oalg_symbionts.scg.chdmit90.fasta
grep -e "^1" -B1 Oalg_symbionts.scg.cdhit90.fasta.clstr | grep '>' | tr -s " " | cut -d " " -f2 | sed 's/>//g' | sed 's/\.\.\.//g' > remove.list
grep -e "^2" Oalg_symbionts.scg.cdhit90.fasta.clstr | grep '>' | tr -s " " | cut -d " " -f2 | sed 's/>//g' | sed 's/\.\.\.//g' >> remove.list
grep -e "^3" Oalg_symbionts.scg.cdhit90.fasta.clstr | grep '>' | tr -s " " | cut -d " " -f2 | sed 's/>//g' | sed 's/\.\.\.//g' >> remove.list
   #^More need to be added as required to the highest clustered number. In this study, it was only required to be 1 and 2. 3 is shown above as an example.
filterbyname.sh in=Oalg_symbionts.scg.cdhit90.fasta out=Oalg_symbionts.scg.cdhit90.single.fasta include=f names=remove.list ow; rm remove.list
#Single-copy gene sequences of symbionts (only singletons >90% ID) are saved as <Oalg_symbionts.scg.cdhit90.single.fasta>.

echo "kallisto (making a kallisto index of the singleton single-copy genes of symbionts; v0.44.0)"
kallisto index -i Oalg_symbionts.scg.cdhit90.single.kdb Oalg_symbionts.scg.cdhit90.single.fasta

echo "kallisto (quantifying symbionts' relative abundance based on single-copy genes; v0.44.0)"
#Abundance estimates are based on kallisto counts against the above single copy marker gene collection (200~400 genes per symbiont).
#Raw reads are first QC-trimmed. After kallisto count, a summary composition table based on reported tpm is made.
#To exclude potentially multiple-copied genes and assess the relative abundance reliably, the mean tpm of the inter-quartile range genes
#  (IQR = 25%ile~75%ile genes in ranking of tpm) is reported.
#This part is looped over all samples listed in <sample_info.list>
for sample in $(cut -f1 sample_info.list); do
	echo "  file copy and unziping"
	cp $DATAPATH/${sample}_R1.fq.gz read1.fastq.gz
	cp $DATAPATH/${sample}_R2.fq.gz read2.fastq.gz
	unpigz read?.fastq.gz

	echo "  QC-trimming reads by trimmomatic v0.36"
	java -jar $TOOLPATH/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 2 -phred33\
	  read1.fastq read2.fastq\
	  R1_paired.fastq R1_unpaired.fastq\
	  R2_paired.fastq R2_unpaired.fastq\
	  LEADING:10 SLIDINGWINDOW:4:2 MINLEN:50
	  #first 10-bases off, window of 4 bases to check quality >2, shorter than 50 to dump.
	rm read?.fastq

	echo "  kallisto to check abundance per reference single copy gene sequence"
	kallisto quant -i Oalg_symbionts.scg.cdhit90.single.kdb -o kallisto -t 2 \
	  --plaintext R1_paired.fastq R2_paired.fastq --pseudobam
	mv kallisto/abundance.tsv $sample.kallisto.scg.raw.output.tsv

	echo "  making a summary table..."
	tail $sample.kallisto.scg.raw.output.tsv -n+2 | sort -k1 | cut -f1 | cut -d "_" -f1 > tmp1   #list symbiont species from headers
	tail $sample.kallisto.scg.raw.output.tsv -n+2 | sort -k1 | cut -f5 > tmp2		     #list TPMs from kallisto output
	paste tmp1 tmp2 > tmp3						#list symbiont and TPM
	for symb in $(uniq tmp1); do
		grep -w $symb tmp3 | sort -nk2 > tmp_s			#sort the list of symbiont and TPM, by TPM
		n_gene=$(cat tmp_s | wc -l)					#number of genes listed
		endpos=$(awk -v n=$n_gene 'BEGIN{printf "%d\n", n * 0.75}')	#end line to keep for IQR
		n_line=$(awk -v n=$n_gene 'BEGIN{printf "%d\n", n * 0.5}')	#number of lines to keep for IQR
		head tmp_s -n $endpos | tail -n $n_line > iqr			#interquartile range is kept to calculate mean coverage as TPM
		awk -v sample=$sample '{sum+=$2} END {if (NR > 0) printf "%s\t %s\t %.3f\n", sample, $1, sum/NR}' iqr  #calculate iqr mean per $symb
	done > tmp4   							#list $sample, $symb, $iqr_mean_TPM
	echo "sampleID	symbiont	meanTPM(IQR)	relative_abundance" > $sample.kallisto.scg.symb.composition #header for a composition table
	awk 'FNR==NR{sum+=$3; next;} {printf "%s\t%s\t%s\t%.3f\t\n",$1,$2,$3,$3/sum}' tmp4 tmp4 >> $sample.kallisto.scg.symb.composition
	echo "$sample.kallisto.scg.symb.composition"
	cat $sample.kallisto.scg.symb.composition | column -t
	echo
done; rm tmp? tmp_s iqr
#After all the runs are finished, concatenate all .composition files to output a summary table.
echo "sampleID	symbiont	meanTPM(IQR)	relative_abundance" > all.kallisto.scg.symb.composition.txt
tail -n+2 *.kallisto.scg.symb.composition >> all.kallisto.scg.symb.composition.txt

echo "Additional check for the distribution of relative coverage (mean+-s.d.) across single-copy genes in each symbiont in each sample."
#Samples are analyzed for a given symbiont if the IQR mean TMP > 5.
awk '$3>5 {printf "%s\t%s\t\n", $1, $2}' all.kallisto.scg.symb.composition > sample_symb.list
  #sample and symbiont are listed if TMP>5.
echo "sample	symbiont	gene	tpm	prop" > all.kallisto.scg.prop.txt  #Header for the output.
for sample in $(cut -f1 sample_info.list); do
	for symb in $(grep $sample sample_symb.list | cut -f2); do
		grep -w $symb $sample.kallisto.scg.raw.output.tsv | cut -f1,5 > $sample.$symb.tmp  #recording gene and tmp
		awk -v sample=$sample -v symb=$symb 'FNR==NR{sum+=$2; next;} {printf "%s\t %s\t %s\t %.3f\t %0.4f\n", sample,symb,$1,$2,$2/sum}' \
		   $sample.$symb.tmp $sample.$symb.tmp >> all.kallisto.scg.prop.txt  #adding relative abundance per gene in each symbiont per sample
		rm $sample.$symb.tmp
	done
done