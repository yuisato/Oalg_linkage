#This script takes paired-end read files and outputs phyloflash analysis results.
#For setting up a phyloflash database based on SILVA SSU database, please refer to http://hrgv.github.io/phyloFlash/
#Input files are $DATAPATH/${sample}_R1.fq.gz and $DATAPATH/${sample}_R2.fq.gz

echo "phyloFlash (taxonomic composition analysis based on SSU; v3.3.beta1"
#Looped over for all samples listed in <sample_info.list>
for sample in $(cut -f1 sample_info.list); do
	readlen=$(zcat ${sample}_R1.fq.gz| head -n2| tail -n1|wc -c) #read length checked
	phyloFlash.pl -lib $sample\
	 -read1 ${sample}_R1.fq.gz -read2 ${sample}_R2.fq.gz \
	 -readlength $readlen -CPUs 2 -id 63
done
#Check $sample.spades*final.fasta for assembled SSU sequences and closest hit for taxonomic assignment.
