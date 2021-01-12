# GToTree phylogenomic tree of Oalg symbionts and 88 reference genome assemblies (v1.5.38)

$TOOLPATH/GToTree\
 -a GToTree_ref-ncbi-acc.list\
 -f GToTree_symb-bin-fasta.list\
 -m GToTree_acc-label.map\
 -H Bacteria_and_Archaea\
 -t -j 8\
 -L Phylum,Class,Species\
 -o GToTree
#Reference genomes were selected based on 16S rRNA gene similarity to the symbionts.
#The label map file indicates to which symbionts reference bins are similar (with serial numbers).
#Indicated values in the map files are published mutation rates in Lynch et al. 2016 Nat Rev Genet.
#An archaeon is included as the outgroup.


#Tree-build with IQTree based on conserved single copy gene AA alignment (v1.5.5)
cd GToTree
iqtree -s run_files/Aligned_SCGs.faa -p run_files/Partitions.txt -B 1000 -m MFP+MERGE -nt 2
cd ..
#Model finder set 'MFP+MERGE' is specified for allowing multiple models per genetic regions, indicated in the partition file.