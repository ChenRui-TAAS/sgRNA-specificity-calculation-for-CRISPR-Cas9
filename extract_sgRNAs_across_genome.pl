#! usr/bin/perl -w
use strict;

# Author:			Rui Chen <chenrui.taas@gmail.com; chenrui.2011@outlook.com>
# Program Date:		2014.07.29	
# Modified:			2015.03.25 
# Modified:			2015.03.26	Optimize memory consumption
# Modified:			2015.04.07	Only for each chromosome
# Modified:			2017.09.26	For TAIR10 ensemblgenomes_TAIR10_37	
#############################################################

## Explanation for cleavage site:
## Screen for 3 levels: genome, gene and exon(CDS)
## Finding CCN and NGG form 5' to 3'
## Cleavage site: 			(NNNNNNNNNNNNNNNNNxNNN)NGG		CCN(NNNxNNNNNNNNNNNNNNNNN)
## Extract region:	20bp	5'-(N20)NGG-3'					5'-CCN(N20 reverse complemented)-3'
## Length(chr): 1 initial
## NGG region:	genome	chr_start	+5&+6	cut-site	+0&+1	
##						chr_end 	-0&-1	cut-site	-5&-6
##				gene	gene_start	+5$+6	cut-site	+0&+1
##						gene_end	+4&+5	cut-site	-0&-1
##				exon	exon_start	+5&+6	cut-site	+0&+1
##						exon_end	+4&+5	cut-site	-0&-1
##
## CCN region:	genome	chr_start	+0&+1	cut-site	+5&+6	
##						chr_end 	-5&-6	cut-site	-0&-1
##				gene	gene_start	-4$-5	cut-site	+0&+1
##						gene_end	-5&-6	cut-site	-0&-1
##				exon	exon_start	-4&-5	cut-site	+0&+1
##						exon_end	-5&-6	cut-site	-0&-1
#############################################################

my $gff3_file = shift;
my $genome_fa = shift;
my $out = shift;

die "
Usage:
perl   extract_sgRNAs_across_genome.pl   gff3_file   genome.fa   output_filename

For example:
gff3_file	Arabidopsis_thaliana.TAIR10.37.gff3
genome.fa	Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
out_file	sgRNA_whole_genome_TAIR10.fa

" if !defined $out;

#############################################################
my(
$count,
@t,
$chr,
@w,
$gene_name,
%chromosome_gene,
%gene_coordinate,
$gene_line_count,
$tag,
%seq,
$i,
$j,
$k,
$r,
%coordinate,
@temp,
$pam_site,
$cut_site_left,
$cut_site_right,
$pam_site_real,
$seed_seq,
$seed_location,
$tmp_seq,
%gene_tag,
$loci,
@t_left,
@t_right,
$x,
$y,
$z
);

##Based on the No.3 column of Arabidopsis gff3 file
my @type =(
"CDS", 
"exon", 
"five_prime_UTR", 
"gene", 
"lnc_RNA", 
"miRNA", 
"mRNA", 
"ncRNA_gene", 
"pre_miRNA", 
"RNase_MRP_RNA", 
"rRNA", 
"snoRNA", 
"snRNA", 
"SRP_RNA", 
"three_prime_UTR", 
"transcript", 
"tRNA" 
);																	

my %re = (
"A" => "T",
"T" => "A",
"C" => "G",
"G" => "C",
"N" => "N",
);
#############################################################


## record gff3 file
open (IN, "< $gff3_file") or die $!;
while (<IN>) {
	$count ++; print "Loading gff3_file lines: $count\r" if $count % 1000 == 0;
	
	chomp;
	@t = split(/\t/, $_);
	
	if($t[0] =~ /^\d+/){		## Chromosome name can not be initialed by number, here for arabidopsis. 
		$chr = "Chr".$t[0];
	}else{
		$chr = $t[0];
	}
	
	if ($t[8]){
		@w = split(";", $t[8]);		##Record all lines. 
		
		if($w[0] =~ /:(AT\w+)/|| $w[0] =~ /:(at\w+)/ || $w[0] =~ /:(ENSRNA\w+)/){	## Modifiy the gene ID for rice or soybean genome, here is for arabidopsis. 
			
			if(!defined $gene_name || $gene_name ne $1){
				$gene_name = $1;
				push (@{$chromosome_gene{$chr}}, $gene_name);
				$gene_line_count = 0;
			}else{
				$gene_line_count ++;
			}
			
			${$gene_coordinate{$gene_name}}[$gene_line_count] = "$t[2]\t$t[3]\t$t[4]\t$t[6]";	## Record information for each gene: locus, start_site, end_site, and orientation.
		}
	}
}
print "\n";
$count = 0;
close IN;


## record genome_fa
open (IN, "< $genome_fa") or die $!;
while (<IN>) {
    $_ =~ s/[\r\n]+//g;
	if($_ =~ /^>(\S+)/){
		$tag = $1;
		if($tag =~ /^\d+/){		## Chromosome name can not be initialed by number, here for arabidopsis. 
			$tag = "Chr".$tag;
		}
		print "Recording genome: $tag\r";
	}else{
		$seq{$tag} .= $_;   ## 0 initialed
	}
}
print "\n";
close IN;


## output sgRNA sequences across genome
die "The output file (".$out.") already exist.\n" if (-e $out);
open (OUT, "> $out") or die $!;
foreach $i (sort keys %seq){
	print "Outputting according to chromosome: $i\r";
	
	foreach $j (@{$chromosome_gene{$i}}){
		for $k (0..$#{$gene_coordinate{$j}}){
			@t = split ("\t", ${$gene_coordinate{$j}}[$k] ); 
			
			for $r ($t[1]..$t[2]){
				if (!defined $coordinate{$i."_".$r."_".$t[0]} || $coordinate{$i."_".$r."_".$t[0]} !~ /$j/){	## "_"
					$coordinate{$i."_".$r."_".$t[0]} .= $j.";";
				}
			}		
		}
	}
	
	undef @temp;
	@temp = split('', $seq{$i});	## 0 initialed
	
	for $j (1..$#temp){
		undef $pam_site;
		undef $cut_site_left;
		undef $cut_site_right;
		undef $pam_site_real;
		undef $seed_seq;
		undef $seed_location;
		undef $tmp_seq;
	
		if($temp[$j-1].$temp[$j] eq "CC"){		## 0 initialed
			$pam_site = "CCN";
			$cut_site_left = $j +5;
			$cut_site_right = $j +6;
			
			if($cut_site_right > $#temp-16){
				next;
			}else{
				$tmp_seq = substr($seq{$i}, $j+2, 20);
				if($tmp_seq !~ /[ATCGN]{20}/){
					next;
				}
				
				$pam_site_real = $temp[$j-1].$temp[$j].$temp[$j+1];
				$seed_seq = base_re (substr($seq{$i}, $j+2, 20));		## 0 initialed
				$seed_location = $i."_-_".($j+3)."_".($j+22);
			}
			
		}elsif($temp[$j-1].$temp[$j] eq "GG"){	## 0 initialed
			$pam_site = "NGG";
			$cut_site_left = $j -5;
			$cut_site_right = $j -4;
			
			if($cut_site_left <= 16){
				next;
			}else{
				$tmp_seq = substr($seq{$i}, $j-22, 20);
				if($tmp_seq !~ /[ATCGN]{20}/){
					next;
				}
				
				$pam_site_real = $temp[$j-2].$temp[$j-1].$temp[$j];
				$seed_seq = substr($seq{$i}, $j-22, 20);				## 0 initialed
				$seed_location = $i."_+_".($j-21)."_".($j-2);
			}
		}
		
		if(defined $cut_site_left && defined $cut_site_right && $cut_site_right <= $#temp-16 && $cut_site_left >= 16){
			next if $seed_seq =~ /N/;  ## Discard sgRNA sequence with ambigous N
			undef %gene_tag;
			undef $loci;
			
			##  
			foreach $k (@type){
				undef @t_left;
				undef @t_right;
				if(defined $coordinate{$i."_".$cut_site_left."_".$k} && defined $coordinate{$i."_".$cut_site_right."_".$k}){
					@t_left  = split(";", $coordinate{$i."_".$cut_site_left."_".$k});
					@t_right = split(";", $coordinate{$i."_".$cut_site_right."_".$k});
					
					foreach $x (@t_left) { 
						foreach $y (@t_right){
							if ($x eq $y){
								$gene_tag{$x} .= $k.";"; 
							}
						}
					}
				}
			}
			
			if(%gene_tag){
				foreach $z (keys %gene_tag){
					$loci .= "\t$z\t$gene_tag{$z}";
				}
			}else{
				$loci = "\tIntergenic_Region\t-";
			}

			print OUT  ">".$seed_location."\t".$pam_site."\t".$pam_site_real."\t".$i."_".$cut_site_left."_".$cut_site_right.$loci."\n".$seed_seq."\n";
		}
	}
}
close OUT;
print "\n";


sub base_re {
    my @re = reverse split('', $_[0]);
	my $s0;
	foreach(@re){
		$s0 .= $re{$_};
	}
    $s0;
}

