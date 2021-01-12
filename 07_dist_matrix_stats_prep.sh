#This scripts prepares data tables based on distance matrices from ngsDist output generated in the step 05a and 06a for further statistical testing.
#First, a list of ngsDist results used for testing is made (tab separated).
echo "Mito	80Samples_to_Mito_MinQ20_SNPpval0.01_Prior2_haploid.dist" > result-set
echo "Gamma1	80Samples_to_Gamma1_MinQ20_SNPpval0.01_Prior2_haploid.dist" >> result-set
echo "Gamma3	80Samples_to_Gamma3_MinQ20_SNPpval0.01_Prior2_haploid.dist" >> result-set
echo "Delta1a	37Samples_Delta1a_60perc_MinQ20_SNPpval0.01_Prior2_haploid.dist" >> result-set
echo "Delta1b	46Samples_Delta1b_65perc_MinQ20_SNPpval0.01_Prior2_haploid.dist" >> result-set
echo "Delta4	67Samples_Delta4_50perc_MinQ20_SNPpval0.01_Prior2_haploid.dist" >> result-set
echo "Spiro	41Samples_Spiro_63perc_MinQ20_SNPpval0.01_Prior2_haploid.dist" >> result-set
#Copy these .dist files in the working folder.

for type in $(cut -f1 result-set); do
	file=$(grep $type result-set | cut -f2)
	
	#Pairwise distances of randomized pairs will be extracted without replacement for statistical analyses.
	#Count the number of samples in each host-group.
	n_cav_a=$(cut -f1 $file | tail -n+3 | grep OalgCAVL_A -c)
	n_cav_b=$(cut -f1 $file | tail -n+3 | grep OalgCAVL_B -c)
	n_san_a=$(cut -f1 $file | tail -n+3 | grep OalgSANT_A -c)
	n_san_b=$(cut -f1 $file | tail -n+3 | grep OalgSANT_B -c)

	echo "Group	Distance" > ${type}_pairdist_test-set.tsv
	echo "Within group"
	echo "  A_CAVL"
	if (( n_cav_a > 1 )); then
	  n_pair=$(echo $n_cav_a | awk '{print int($1/2)}')  #Number of pairs within the host group.
	  for i in $(seq 1 $n_pair); do
		echo A_CAVL >> tmp1
	  done
	  for j in $(seq 1 $n_cav_a); do echo "$j	$RANDOM"; done | sort -nk2 | cut -f1 > rand  #randomize numbers
	  for l in $(seq 1 $n_pair); do
		m=$(tail rand -n+$((2 * $l - 1)) | head -n1)
		n=$(tail rand -n+$((2 * $l)) | head -n1)	  #taking two random numbers.
		cut -f $((1 + $m)) $file | tail -n+$((2 + $n)) | head -n1 >> tmp2   
		  #take a value from a random (m, n) position in the A_CAVL x A_CAVL region of the matrix. 
	  done
	  paste tmp1 tmp2 >> ${type}_pairdist_test-set.tsv
	  rm tmp1 tmp2 rand
	fi

	echo "  B_CAVL"
	if (( n_cav_b > 1 )); then
	  n_pair=$(echo $n_cav_b | awk '{print int($1/2)}')
	  for i in $(seq 1 $n_pair); do
		echo B_CAVL >> tmp1
	  done
	  for j in $(seq 1 $n_cav_b); do echo -e "$j\t$RANDOM"; done | sort -nk2 | cut -f1 > rand
	  for l in $(seq 1 $n_pair); do
		m=$(tail rand -n+$((2 * $l - 1)) | head -n1)
		n=$(tail rand -n+$((2 * $l)) | head -n1)
		cut -f $((1 + $n_cav_a + $m)) $file | tail -n+$((2 + $n_cav_a + $n)) | head -n1 >> tmp2
		  #take a value from a random (m, n) position in the B_CAVL x B_CAVL region of the matrix. 
	  done
	  paste tmp1 tmp2 >> ${type}_pairdist_test-set.tsv
	  rm tmp1 tmp2 rand
	fi

	echo "  A_SANT"
	if (( n_san_a > 1 )); then
	  n_pair=$(echo $n_san_a | awk '{print int($1/2)}')
	  for i in $(seq 1 $n_pair); do
		echo A_SANT >> tmp1
	  done
	  for j in $(seq 1 $n_san_a); do echo -e "$j\t$RANDOM"; done | sort -nk2 | cut -f1 > rand
	  for l in $(seq 1 $n_pair); do
		m=$(tail rand -n+$((2 * $l - 1)) | head -n1)
		n=$(tail rand -n+$((2 * $l)) | head -n1)
		cut -f $((1 + $n_cav_a + $n_cav_b + $m)) $file | tail -n+$((2 + $n_cav_a + $n_cav_b + $n)) | head -n1 >> tmp2
		  #take a value from a random (m, n) position in the A_SANT x A_SANT region of the matrix. 
	  done
	  paste tmp1 tmp2 >> ${type}_pairdist_test-set.tsv
	  rm tmp1 tmp2 rand
	fi

	echo "  B_SANT"
	if (( n_san_b > 1 )); then
	  n_pair=$(echo $n_san_b | awk '{print int($1/2)}')
	  for i in $(seq 1 $n_pair); do
		echo B_SANT >> tmp1
	  done
	  for j in $(seq 1 $n_san_b); do echo -e "$j\t$RANDOM"; done | sort -nk2 | cut -f1 > rand
	  for l in $(seq 1 $n_pair); do
		m=$(tail rand -n+$((2 * $l - 1)) | head -n1)
		n=$(tail rand -n+$((2 * $l)) | head -n1)
		cut -f $((1 + $n_cav_a + $n_cav_b + $n_san_a + $m)) $file | tail -n+$((2 + $n_cav_a + $n_cav_b + $n_san_a + $n)) | head -n1 >> tmp2
		  #take a value from a random (m, n) position in the B_SANT x B_SANT region of the matrix. 
	  done
	  paste tmp1 tmp2 >> ${type}_pairdist_test-set.tsv
	  rm tmp1 tmp2 rand
	fi


	echo "Between COI-haplotypes (Same location)"
	echo "  AvsB_CAVL"
	if (( n_cav_a > 0 & n_cav_b > 0 )); then
	  if (( n_cav_a <= n_cav_b )); then
		n_pair=$n_cav_a
		for i in $(seq 1 $n_pair); do
		  echo AvsB_CAVL >> tmp1
		done
		for j in $(seq 1 $n_cav_b); do echo -e "$j\t$RANDOM"; done | sort -nk2 | cut -f1 > rand
		for l in $(seq 1 $n_pair); do
		  m=$(tail rand -n+$l | head -n1)
		  cut -f $((1 + $l)) $file | tail -n+$((2 + $n_cav_a + $m)) | head -n1 >> tmp2
  		    #take l x values from each of random (m) position in the A_CAVL x B_CAVL region of the matrix. 
		done
		paste tmp1 tmp2 >> ${type}_pairdist_test-set.tsv
	  else
		n_pair=$n_cav_b
		for i in $(seq 1 $n_pair); do
		  echo AvsB_CAVL >> tmp1
		done
		for j in $(seq 1 $n_cav_a); do echo -e "$j\t$RANDOM"; done | sort -nk2 | cut -f1 > rand
		for l in $(seq 1 $n_pair); do
		  m=$(tail rand -n+$l | head -n1)
		  cut -f $((1 + $n_cav_a + $l)) $file | tail -n+$((2 + $m)) | head -n1 >> tmp2
  		    #take l x values from each of random (m) position in the A_CAVL x B_CAVL region of the matrix. 
		done
		paste tmp1 tmp2 >> ${type}_pairdist_test-set.tsv
	  fi
	  rm tmp1 tmp2 rand
	fi

	echo "  AvsB_SANT"
	if (( n_san_a > 0 & n_san_b > 0 )); then
	  if (( n_san_a <= n_san_b )); then
		n_pair=$n_san_a
		for i in $(seq 1 $n_pair); do
		  echo AvsB_SANT >> tmp1
		done
		for j in $(seq 1 $n_san_b); do echo -e "$j\t$RANDOM"; done | sort -nk2 | cut -f1 > rand
		for l in $(seq 1 $n_pair); do
		  m=$(tail rand -n+$l | head -n1)
		  cut -f $((1 + $n_cav_a + $n_cav_b + $l)) $file | tail -n+$((2 + $n_cav_a + $n_cav_b + $n_san_a + $m)) | head -n1 >> tmp2
  		    #take l x values from each of random (m) position in the A_SANT x B_SANT region of the matrix. 
		done
		paste tmp1 tmp2 >> ${type}_pairdist_test-set.tsv
	  else
		n_pair=$n_san_b
		for i in $(seq 1 $n_pair); do
		  echo AvsB_SANT >> tmp1
		done
		for j in $(seq 1 $n_san_a); do echo -e "$j\t$RANDOM"; done | sort -nk2 | cut -f1 > rand
		for l in $(seq 1 $n_pair); do
		  m=$(tail rand -n+$l | head -n1)
		  cut -f $((1 + $n_cav_a + $n_cav_b + $n_san_a + $l)) $file | tail -n+$((2 + $n_cav_a + $n_cav_b + $m)) | head -n1 >> tmp2
  		    #take l x values from each of random (m) position in the A_SANT x B_SANT region of the matrix. 
		done
		paste tmp1 tmp2 >> ${type}_pairdist_test-set.tsv
	  fi
	  rm tmp1 tmp2 rand
	fi
	
	echo "Between locations (Same COI-haplotype)"
	echo "  SAvsCA_A"
	if (( n_cav_a > 0 & n_san_a > 0 )); then
	  if (( n_cav_a <= n_san_a )); then
		n_pair=$n_cav_a
		for i in $(seq 1 $n_pair); do
		  echo SAvsCA_A >> tmp1
		done
		for j in $(seq 1 $n_san_a); do echo -e "$j\t$RANDOM"; done | sort -nk2 | cut -f1 > rand
		for l in $(seq 1 $n_pair); do
		  m=$(tail rand -n+$l | head -n1)
		  cut -f $((1 + $l)) $file | tail -n+$((2 + $n_cav_a + $n_cav_b + $m)) | head -n1 >> tmp2
  		    #take l x values from each of random (m) position in the A_SANT x A_CAVL region of the matrix. 
		done
		paste tmp1 tmp2 >> ${type}_pairdist_test-set.tsv
	  else
		n_pair=$n_san_a
		for i in $(seq 1 $n_pair); do
		  echo SAvsCA_A >> tmp1
		done
		for j in $(seq 1 $n_cav_a); do echo -e "$j\t$RANDOM"; done | sort -nk2 | cut -f1 > rand
		for l in $(seq 1 $n_pair); do
		  m=$(tail rand -n+$l | head -n1)
		  cut -f $((1 + $n_cav_a + $n_cav_b + $l)) $file | tail -n+$((2 + $m)) | head -n1 >> tmp2
  		    #take l x values from each of random (m) position in the A_SANT x A_CAVL region of the matrix. 
		done
		paste tmp1 tmp2 >> ${type}_pairdist_test-set.tsv
	  fi
	  rm tmp1 tmp2 rand
	fi

	echo "  SAvsCA_B"
	if (( n_san_b > 0 & n_cav_b > 0 )); then
	  if (( n_cav_b <= n_san_b )); then
		n_pair=$n_cav_b
		for i in $(seq 1 $n_pair); do
		  echo SAvsCA_B >> tmp1
		done
		for j in $(seq 1 $n_san_b); do echo -e "$j\t$RANDOM"; done | sort -nk2 | cut -f1 > rand
		for l in $(seq 1 $n_pair); do
		  m=$(tail rand -n+$l | head -n1)
		  cut -f $((1 + $n_cav_a + $l)) $file | tail -n+$((2 + $n_cav_a + $n_cav_b + $n_san_a + $m)) | head -n1 >> tmp2
  		    #take l x values from each of random (m) position in the B_SANT x B_CAVL region of the matrix. 
		done
		paste tmp1 tmp2 >> ${type}_pairdist_test-set.tsv
	  else
		n_pair=$n_san_b
		for i in $(seq 1 $n_pair); do
		  echo SAvsCA_B >> tmp1
		done
		for j in $(seq 1 $n_cav_b); do echo -e "$j\t$RANDOM"; done | sort -nk2 | cut -f1 > rand
		for l in $(seq 1 $n_pair); do
		  m=$(tail rand -n+$l | head -n1)
		  cut -f $((1 + $n_cav_a + $n_cav_b + $n_san_a + $l)) $file | tail -n+$((2 + $n_cav_a + $m)) | head -n1 >> tmp2
  		    #take l x values from each of random (m) position in the B_SANT x B_CAVL region of the matrix. 
		done
		paste tmp1 tmp2 >> ${type}_pairdist_test-set.tsv
	  fi
	  rm tmp1 tmp2 rand
	fi
done

#Each group is categorized and combined into 'Within', 'AvsB', 'SAvsCA' catergories for all references.
for table in *.tsv; do 
	sed "s/^A_CAVL/Within/g" < $table | sed "s/^B_CAVL/Within/g" | sed "s/^A_SANT/Within/g" | sed "s/^B_SANT/Within/g" | sed "s/_CAVL//g" |\
	 sed "s/_SANT//g" |sed "s/_A//g" | sed "s/_B//g" > $(echo $table | sed 's/.tsv/_combined.tsv/')
done
#New tables <${type}_pairdist_test-set_combined.tsv> are used for statistical tests. 