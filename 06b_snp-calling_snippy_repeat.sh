#This script repeats building of SNP-based phylogeny of symbionts using a deterministic genotyping approach 
# when some samples need to be removed from the test set due to no SNPs identified.

#Filtering of samples is based on horizontal coverage (% reference base covered with reads => 5X) for a given symbiont. 
#Horizontal coverage is checked in <snp_${ref}/core_${ref}_Oalg_all.txt>
#<snp_Gamma3/core_Gamma3_Oalg_all.txt> Example
# ID	Aligned Bases	Reference bases	% Aligned Bases
# gamma3_OalgCAVL_A11	4097589	4272394	95.91
# gamma3_OalgCAVL_A12	4120596	4272394	96.45
# gamma3_OalgCAVL_A13	4142224	4272394	96.95
# ...

ref=Gamma3 #etc.  Change as needed.
cutoff=80 #etc. (Horizontal coverage % with =>5X reads mapped) Individually set for the given symbiont.
 
echo "copying snippy snp-folders of filtered samples based on horizontal coverage above $cutoff per cent"
tail -n+2 snp_${ref}/core_${ref}_Oalg_all.txt | while IFS= read -r line; do
	cover=$(echo $line | cut -d" " -f4)
	echo $line | awk -v cov=$cover -v cutoff=$cutoff 'cov>cutoff {print $1}' >> sample.list
done
n_ind=$(cat sample.list | wc -l) #Number of filtered samples. 
mkdir snippy_repeat				 #Temporal folder to store snp data.
for sample in $(cat sample.list); do
	cp -r snp_${ref}/$sample snippy_repeat
done

echo "snippy-core running"
snippy-core --prefix=core_${ref}_${cutoff}perc snippy_repeat/* --noref

echo "iqtree-buidling for coreSNP tree..."
iqtree -nt 1 -alrt 1000 -bb 1000 -s core_${ref}_${cutoff}perc.aln -m GTR+ASC -redo
 #GTR+ASC model for core SNP trees
echo "If no SNPs were found. Adjust the horizontal coverage cutoff value."

rm -r snippy_repeat sample.list