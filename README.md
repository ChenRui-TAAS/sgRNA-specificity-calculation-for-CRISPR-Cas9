sgRNA-specificity-calculation-for-CRISPR-Cas9
=============================================

A simple pipeline for genome-wide calculating sgRNA specificity and off-target number for the CRISPR-Cas9 system

2019.09.03

Rui Chen

chenrui.taas@gmail.com; chenrui.2011@outlook.com


Brief introduction
------------------
Clustered regularly interspaced short palindromic repeats (CRISPR), a class of immune-associated sequences in bacteria, have been developed as a powerful tool for editing eukaryotic genomes in diverse cells and organisms in recent years. The CRISPR-Cas9 system can recognize upstream 20 nt sequence (guide sequence) adjacent to the PAM site, and trigger double strand DNA cleavage and DNA repair mechanisms, which eventually result in knock-out, knock-in or site-specific mutagenesis. However, off-target effect caused by guide sequence misrecognition is the major drawback and restricts its widespread application. In this study, global analysis of specificities of all guide sequences in _Arabidopsis thaliana_, _Oryza sativa_ (rice) and _Glycine max_ (soybean) genomes were performed and three genome-wide databases for these plants were established including intergenic regions. For each target site of CRISPR-Cas9, specificity score and off-target number were calculated and evaluated. The mean values of off-target number for A. thaliana, rice and soybean were determined as 27.5, 57.3 and 174.7, respectively. Comparative analysis among these plants suggested that the frequency of off-target effect was positively correlated to genome size. Our results contributed to the better understanding of CRISPR-Cas9 system in plants and helped to minimize the off-target effect during its applications in the future.


Version of plant genomes:
-------------------------
_Arabidopsis thaliana_ genome (TAIR10) was downloaded from Ensembl (http://www.ensembl.org).

_Oryza sativa_ genome (IRGSP-1.0) was downloaded from NCBI (GCF_001433935.1).

_Glycine max_ genome (Gmax_275_v2.0) and annotation files were downloaded from Phytozome (https://phytozome.jgi.doe.gov, Gmax V10).


Manual for sgRNA analysis
-------------------------
###Step_1, get whole set of sgRNA sequences.

"extract_sgRNAs_across_genome.pl"

Usage:
perl   extract_sgRNAs_across_genome.pl   gff3_file   genome.fa   output_filename

For example:


gff3_file	Arabidopsis_thaliana.TAIR10.37.gff3
genome.fa	Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
out_file	sgRNA_whole_genome_TAIR10.fa

###Step_2, exclude repeat sgRNAs.

"build_orphan_list.pl"

Usage:	
perl   build_orphan_list.pl    input_sgRNA.fa    output_filename

For example:
input_sgRNA.fa	sgRNA_whole_genome_TAIR10.fa
output_filename	sgRNA_whole_genome_TAIR10_orphan.fa


###Step_3, extend sgRNA sequence to 23nt including PAM site.

"reform_sgRNA.pl"

Usage: 
perl   reform_sgRNA.pl   orphan.fa   output_filename

For example:
orphan.fa	sgRNA_whole_genome_TAIR10_orphan.fa
output_filename	sgRNA_whole_genome_TAIR10_orphan_23nt.fa

###Step_4, genomic mapping using BatMis.

build_index	Arabidopsis_thaliana.TAIR10.dna.toplevel.fa

batman \
-g Arabidopsis_thaliana.TAIR10.dna.toplevel.fa \
-q sgRNA_whole_genome_TAIR10_orphan_23nt.fa \
-o batman_out_TAIR10_orphan_23nt \
-n 6 \
-mall

batdecode \
--zip \
-g Arabidopsis_thaliana.TAIR10.dna.toplevel.fa \
-i batman_out_TAIR10_orphan_23nt \
-o batman_out_TAIR10_orphan_23nt.txt.gz


###Step_5, scoring sgRNAs

"scoring_off_targets_and_specificities.pl"

Usage:
perl   scoring_off_targets_and_specificities.pl   batmanOut.txt.gz   scoringOut.txt

For example:
batmanOut.txt.gz        batman_out_TAIR10_orphan_23nt.txt.gz
scoringOut.txt          scoringOut_TAIR10_orphan_23nt.txt

###Step_6, label back to sgRNAs.fa

Usage:	
perl   label_back_sgRNA.pl   sgRNA_orphan.fa   score_table   output_filename

For example:
sgRNA_orphan.fa	sgRNA_whole_genome_TAIR10_orphan.fa
score_table	scoringOut_TAIR10_orphan_23nt.txt
output_filename	sgRNA_whole_genome_TAIR10_orphan_scoreLabeled.fa

Herein, "sgRNA_whole_genome_TAIR10_orphan_scoreLabeled.fa" is the generated DATABASE for sgRNAs.

###Step_7, extract and sort sgRNAs according to appointed genes

"grep_sgRNA_via_locusTable.pl"

Usage:	
perl	grep_sgRNA_via_locusTable.pl   locus_table   score_table	 output_filename

For example:
locus_table	locus_table_ncrna
score_table	sgRNA_whole_genome_TAIR10_orphan_scoreLabeled.fa
output_filename	locus_table_ncrna_sgRNA_out.xls

Content of locus_table_ncrna file: 
AT1G05853.1	Chr1	8525104	8525204
AT2G04455.1	Chr2	852860	853054
AT5G09945.1	Chr5	26966885	26967079
AT4G04185.1	Chr4	1114529	1114629
AT4G06405.1	Chr4	9370787	9370841
AT3G08735.1	Chr3	21052801	21052994
AT5G08810.1	Chr5	23448523	23448717
AT3G56825.1	Chr3	21041727	21042609
AT5G02255.1	Chr5	4690371	4690412
AT4G07120.1	Chr4	12225622	12225786
AT5G04085.1	Chr5	8972466	8972618
AT3G57645.1	Chr3	21346342	21347204
...


At last, you can obtain a sgRNA table according to each gene. 
Ten sgRNAs that have the most high specificities will be listed with detailed information.
Then, you can pick up desired sgRNA candidates for knock-out CRISPR-Cas9 experiment.
