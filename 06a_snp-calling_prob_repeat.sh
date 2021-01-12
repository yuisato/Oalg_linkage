#This script repeats building of phylogeny of symbionts using a probabilistic SNP-identification approach 
# when some samples need to be removed from the test set due to no SNPs identified.

#Filtering of samples is based on horizontal coverage (% reference base covered with reads => 1X) for a given symbiont. 
#The Horizontal coverage is first surveyed from results of <05b_snp-calling_snippy.sh>.
#This is done by counting "-" in the snippy output 'snps.aligned.fa' where no base was covered.
#('N' marks bases that were mapped less than the mincov depth 5X).
for folder in snp*; do
	cd $folder
	for subfolder in *; do
		echo $subfolder >> list1			#${ref}_${sample} (e.g. Delta1a_OalgCAVL_A11)
		cat $subfolder/snps.aligned.fa | tr -cd - | wc -c >> list2
	done
	paste list1 list2 > depth0.base.counts  #saved in each snp_${ref} folder, containing counts from all samples.
	rm list1 list2
	cd ..
done
paste snp*/depth0.base.counts > depth0.base.counts  #Summary of the number of 0-depth bases for all samples for all reference bins.
rm    snp*/depth0.base.counts
#Output 'depth0.base.counts' is used to calculate horizontal coverage manually as referring 02b output <all_horizontal_coverage_snippy.txt>.
#Deriving file (one per reference) 'mapped_percent.txt' is saved in separate folders e.g. "snp_gamma1/mapped_percent.txt"
#'mapped_percent.txt' Example
#ID	Unmapped Bases	Reference bases	% Mapped Bases
#alpha_OalgBEDO_H01	633777	633777	0.00
#alpha_OalgCAVL_A11	633777	633777	0.00
#alpha_OalgCAVL_B31	56269	633777	91.12
#alpha_OalgCAVL_D12	584931	633777	7.71
#...
#Based on this table, an individual cutoff for the horizontal coverage is set for each symbiont as needed below.  


#Script below here is run for each symbiont that requires sample-filtering.
reffile=Gamma3 #etc. Change this as needed.
cutoff=65 #etc. Set individually for each of the symbiont set above (%).

echo "Filter samples based on horizontal coverage above $cutoff per cent and <snp_$reffile/mapped_percent.txt>"
sed "s/${reffile}_//g" snp_$reffile/mapped_percent.txt | tail -n+2 | while IFS= read -r line; do
	cover=$(echo $line | cut -d " " -f4)
	echo $line | awk -v cov=$cover -v cutoff=$cutoff 'cov>cutoff {print $1}' >> sample.list
done
n_ind=$(cat sample.list | wc -l)
mkdir mapped/repeat #temporary folder to store bam and bai files for this process.
for sample in $(cat sample.list); do
	cp mapped/$reffile/$sample*realigned.bam mapped/repeat
	cp mapped/$reffile/$sample*realigned.bai mapped/repeat
done
ls *realigned.bam > bam_list.txt
n_ind=$(cat bam_list.txt | wc -l)          #Calculate the number of individual samples
#cp ref/Gamma*fasta ref/Delta*.fasta ref/Spiro.fasta .
#copy input files and reference files in the working directory

echo "ANGSD step1 for depth-count... (v0.929)"
sleep 60; samtools faidx $reffile.fasta  #this has to be freshly made...
angsd -P 24 -b bam_list.txt -out ${n_ind}Samples_to_${reffile}_MinQ20_qc -ref $reffile.fasta\
 -trim 0 -remove_bads 0 -uniqueOnly 0 -only_proper_pairs 1 -doQsDist 1 -doDepth 1 -doCounts 1 -maxDepth 8000

#Obtaining the depth range from ANGSD-1 output for the next ANGSD-2 step
mindepth=$(sed -e '1,/Global_depth/d' ${n_ind}Samples_to_${reffile}_MinQ20_qc.info| awk '$2>0.05 {print $1}'| head -n1) #15 etc.  5 percentile
maxdepth=$(sed -e '1,/Global_depth/d' ${n_ind}Samples_to_${reffile}_MinQ20_qc.info| awk '$2<0.95 {print $1}'| tail -n1) #85 etc. 95 percentile
for i in $(seq 1 $n_ind); do echo "1.00000" >> all_1.indF; done
  #indF; individual inbreeding coefficients are all set as 1 for haploids (i.e. all are set as 'homozygote').
	  
echo "ANGSD step2 for processing with inbreeding adjustment for haploids..."
samtools faidx $reffile.fasta
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
	echo "  No SNPs were found. Adjust the horizontal coverage cutoff value."
fi

rm -r mapped/repeat sample.list bam_list.txt