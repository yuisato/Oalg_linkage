#This script surveys relative abundance of symbionts based on 16S rRNA (SSU) gene sequences from NCBI database using Kallisto.
#
echo "making a kallisto index from reference SSU sequences (Kallisto v0.44.0)"
cd ref
kallisto index -i Oalg_symbionts_ssu.kdb Oalg_symbionts_ssu.fasta
 #Oalg_symbionts_ssu.fasta contains AJ620502, AJ620497, AJ620496, AM493254, AF328857, AF328856, KP943802 from NCBI.
 #sequence headers:
 # >Gamma1_in_Olavius_algarvensis_AF328856.1
 # >Gamma3_in_Olavius_algarvensis_AJ620496.1
 # >Delta1_in_Olavius_algarvensis_AF328857.1
 # >Delta3_in_Olavius_algarvensis_AM493254.1
 # >Delta4_in_Olavius_algarvensis_AJ620497.1
 # >Spirochete_in_Olavius_algarvensis_AJ620502.2
 # >x-Oalg_Olavius_algarvensis_18S_KP943802.1
cd ..

echo "kallisto (quantifying symbionts' relative abundance based on SSU genes; v0.44.0)"
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

	echo "  kallisto to check abundance per reference SSU gene sequence"
	kallisto quant -i Oalg_symbionts_ssu.kdb -o kallisto -t 2 \
	  --plaintext R1_paired.fastq R2_paired.fastq --pseudobam
	mv kallisto/abundance.tsv $sample.kallisto.ssu.raw.output.tsv

	echo "  making a summary table..."
	tail $sample.kallisto.ssu.raw.output.tsv -n+2 | sort -k1 | cut -f1 | cut -d "_" -f1 > tmp1   #list symbiont species from headers
	tail $sample.kallisto.ssu.raw.output.tsv -n+2 | sort -k1 | cut -f5 > tmp2		     		 #list TPMs from kallisto output
	paste tmp1 tmp2 > tmp3
	for symb in $(cat tmp1 | uniq); do
		grep $symb tmp3 > tmp_s			#sort the list of symbiont and TPM, by TPM
		awk -v sample=$sample 'FNR==NR{sum+=$2; next;} {printf "%s\t%s\t%.3f\n", sample, $1, sum}' tmp_s tmp_s | head -n1
	done > tmp4							#this part is for multiple SSUs per symbiont. For development...
	echo "sampleID	symbiont_clade	TPM	relative_to_host-SSU	relative_in_symbiont-SSU" > $sample.kallisto.ssu.symb.composition #header
	grep -v 'x-' tmp4 > symb
	host=$(grep 'x-' tmp4 | cut -f3)
	awk -v host=$host 'FNR==NR{sum+=$3; next;} {printf "%s\t%s\t%s\t%.3f\t%.3f\n",$1,$2,$3,$3/host,$3/sum}' symb symb\
		>> $sample.kallisto.ssu.symb.composition
done

#After all the runs are finished, concatenate all .composition files to output a summary table.
echo "sampleID	symbiont_clade	TPM	relative_to_host-SSU	relative_in_symbiont-SSU" > all.kallisto.ssu.symb.composition.txt
tail -n+2 *.kallisto.ssu.symb.composition > all.kallisto.ssu.symb.composition.txt


