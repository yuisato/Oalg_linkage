#This script aims to assemble a full circular mitochondrial genome from raw fastq data and a publicly available reference COI-PCR sequence.

echo "fastuniq (PCR-deduplication; v1.1)"
cp $DATAPATH/${sample}_R1.fq.gz R1.fq.gz #sample="OalgSANT_A04" etc. Data file names need this formating.
cp $DATAPATH/${sample}_R2.fq.gz R2.fq.gz
unpigz R?.fq.gz; ls R?.fq > seqfile.list     #make an input-file list of unzipped reads for fastuniq
fastuniq -i seqfile.list -t q -o FastUniq1.fastq -p FastUniq2.fastq
rm seqfile.list R?.fq

echo "trimmomatic (adapter-clipping; v0.36)"
java -jar $TOOLPATH/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 4 -phred33\
 FastUniq1.fastq FastUniq2.fastq\
 AdapterClipped_forward_paired.fastq AdapterClipped_forward_unpaired.fastq\
 AdapterClipped_reverse_paired.fastq AdapterClipped_reverse_unpaired.fastq\
 ILLUMINACLIP:$TOOLPATH/bbmap/resources/adapters.fa:2:30:10:4:true
  seed mitmatch=2; palindrome clip threshold=30; simple clip threshold=10; min adapter length=4
rm FastUniq1.fastq FastUniq2.fastq
rm AdapterClipped_forward_unpaired.fastq AdapterClipped_reverse_unpaired.fastq

echo "trimmomatic (QC-trimming; v0.36; separated here for assessment purpose)"
java -jar $TOOLPATH/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 4 -phred33\
 AdapterClipped_forward_paired.fastq AdapterClipped_reverse_paired.fastq\
 QCed_forward_paired.fastq QCed_forward_unpaired.fastq\
 QCed_reverse_paired.fastq QCed_reverse_unpaired.fastq\
 LEADING:10 SLIDINGWINDOW:4:2 MINLEN:50
 #first 10-bases off, window of 4 bases to check quality >2, shorter than 50 to dump.
rm AdapterClipped_forward_paired.fastq AdapterClipped_reverse_paired.fastq
rm QCed_forward_unpaired.fastq QCed_reverse_unpaired.fastq

echo "MITObim (mapping mode to generate a reference mtDNA scaffold; v1.9)"
interleave-fastqgz-MITOBIM.py QCed_forward_paired.fastq QCed_reverse_paired.fastq > interleaved.fastq
MITObim.pl -start 1 -end 100 -sample Oalg -quick ref/Oalgarvensis_COI_KP943854.fasta -ref refCOI \
	-readpool interleaved.fastq -pair -clean
 #reference COI sequence file <Oalgarvensis_COI_KP943854.fasta> needs to be in ref/ folder.
echo ">Oalg-refCOI_mitobim_mapping" > mitobim_out.fasta #prepare a better sequence header for the resulting fasta.
tail -n+2 */*noIUPAC.fasta >> mitobim_out.fasta         #result sequence copied to the final file.

echo "bbmap (identify mitochondrial reads using the mtDNA scaffold from MITObim; v36.86)"
bbmap.sh in1=QCed_forward_paired.fastq in2=QCed_reverse_paired.fastq ref=mitobim_out.fasta\
 nodisk ow threads=4 minid=0.90 killbadpairs pairedonly mappedonly maxindel=1000 pairlen=1000\
 out1=mito1.fastq out2=mito2.fastq

echo "spades (assembly of mtDNA; v3.9.1)"
spades.py --pe1-1 mito1.fastq --pe1-2 mito2.fastq -o spades -k 21,33,55,77,99,127 -t 4
#check 'spades/assembly_graph.fastg' on Bandage to identify and output a circular DNA sequence.
