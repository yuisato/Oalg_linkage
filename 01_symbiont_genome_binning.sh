#This script aims to assemble bacterial symbiont genome bins from raw fastq files.

echo "bbduk (adapter trim left; bbmap v36.86)"
#zcat $TOOLPATH/bbmap/resources/phix174_ill.ref.fa.gz > phix-and-adaptors.fa
#cat $TOOLPATH/bbmap/resources/adapters.fa >> phix-and-adaptors.fa    #reference file with PhiX and Illumina adaptors
bbduk.sh in=$DATAPATH/Oalg_R1.fq.gz in2=$DATAPATH/Oalg_R2.fq.gz ref=phix-and-adaptors.fa \
	out=Oalg_temp_R1.fq.gz out2=Oalg_temp_R2.fq.gz \
	mink=11 minlength=36 ktrim=l hdist=1 threads=16
echo "bbduk (adapter trim right + quality filtering)"
bbduk.sh in=Oalg_temp_R1.fq.gz in2=Oalg_temp_R2.fq.gz ref=phix-and-adaptors.fa \
	out=Oalg_bbduk_R1.fq.gz out2=Oalg_bbduk_R2.fq.gz \
	mink=11 minlength=36 ktrim=r hdist=1 trimq=2 qtrim=rl threads=16

echo "bayeshammer in spades (sequencing error correction; spades v3.9.1)"
spades.py -1 Oalg_bbduk_R1.fq.gz -2 Oalg_bbduk_R2.fq.gz -o bayeshammer \
	--only-error-correction --threads 16
mv bayeshammer/corrected/Oalg*R1* Oalg_processed_R1.fq.gz
mv bayeshammer/corrected/Oalg*R2* Oalg_processed_R2.fq.gz

echo "MEGAHIT (assembly; metagenome mode; v1.0.6)"
megahit -1 Oalg_processed_R1.fq.gz -2 Oalg_processed_R2.fq.gz --presets meta \
	--out-prefix Oalg --min-contig-len 500 --num-cpu-threads 24
mv megahit_out/Oalg.contigs.fa Oalg.megahit.contigs.fa
megahit_toolkit contig2fastg 141 megahit_out/intermediate_contigs/k141.contigs.fa \
	> Oalg.megahit.k141.fastg	#assembly graph to visualise in Bandage

echo "pre-mapping for metabat"
bbmap.sh in=Oalg_processed_R1.fq.gz in2=Oalg_processed_R2.fq.gz ref=Oalg.megahit.contigs.fa \
	outm=mapped.bam nodisk ow unpigz pigz threads=4
samtools sort mapped.bam -o Oalg.contigs.mapped.sorted.bam -@ 2

echo "metabat (binning; v0.26.3)"
$TOOLPATH/metabat/jgi_summarize_bam_contig_depths \
	--outputDepth Oalg.contigs.depth.txt Oalg.contigs.mapped.sorted.bam
metabat -i Oalg.megahit.contigs.fa -a Oalg.contigs.depth.txt -o metabat_bin --sensitive -v \
	--saveTNF assembly_bin_saved.tnf --saveDistance assembly_bin_saved.dist --numThreads 24

echo "checkM (bin quality analysis; v1.0.7)"
mkdir bins; mv metabat_bin.*.fa bins #this folder contains bins from metabat
checkm lineage_wf -t 4 -f Oalg.metabat.bin.checkm -x fa --nt bins checkm
checkm bin_qa_plot --image_type pdf -x fa checkm bins .

#metabat bins were taxonomically assigned to symbiont species using 16S rRNA gene sequences.
#Further refinement of the matabat bins was performed by visualizing <Oalg.megahit.k141.fastg> with Bandage.
#checkM analysis was repeated on the final revised bins.

echo "quast (bin contig stats; v5.0.2)"
quast.py $FINALBIN.fasta -m 1
