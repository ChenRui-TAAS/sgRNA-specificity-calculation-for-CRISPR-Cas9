#! usr/bin/perl -w
use strict;

# Author:		Rui Chen <chenrui.taas@gmail.com; chenrui.2011@outlook.com>
# Program Date:		2017.04.13	
# Version:		1.0		
############################################################

my $locus_table = shift;
my $score_table = shift;
my $out = shift;

die "
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

" if !defined $out;

`date` =~ /(\d+).+\s+(\d+).+\s+(\d+).+\s+(\d+):(\d+):(\d+)\s+/;
my $date = $1.$2.$3.$4.$5;

my(@t, @t2, @t3, @t4, $count, $tag, $i, $j, $k, %locus_line, %gene_gRNA, %geneID, %test);
my($chr, $start, $end, $left, $right, $gene_tag, $gRNA_tag, @temp);
my(%gRNA_line, %gRNA_seq, %gRNA_specific, %gRNA_offtargets, %gRNA_PAM_type, %gRNA_PAM_site, %gRNA_cleavage_site, %gRNA_annotation);
my(%start_site, %end_site, %new_gene_gRNA);
############################################################


open (IN, "< $locus_table") or die $!;
while (<IN>) {
	$_ =~ s/[\r\n]+//g;
	
	undef @t;
	@t = split("\t", $_);
	
	$count ++; 
	print "Recording table lines: $count\r" if $count % 10 == 0;
		
	$gene_tag = $t[0];
	$chr = $t[1];
	$start = $t[2];
	$end = $t[3];
	
	$locus_line{$gene_tag} = $chr."_".$start."_".$end;
	
	for $i($start..$end){
		$geneID{$chr."_".$i} .= $gene_tag.";";
	}
}
close IN;
print "\n";
undef $count;


open (IN, "< $score_table") or die $!;
while (<IN>) {
	$_ =~ s/[\r\n]+//g;
	
	if(/^>(\S+)/){
		$count ++; 
		print "Recording score_table: $count\r" if $count % 1000 == 0;
		
		$gRNA_tag = $1;
		$gRNA_line{$gRNA_tag} = $_;
		
		undef @t;
		@t = split("\t", $_);
		
		$gRNA_PAM_type{$gRNA_tag} 		= $t[1];
		$gRNA_PAM_site{$gRNA_tag} 		= $t[2];
		$gRNA_cleavage_site{$gRNA_tag} 	= $t[3];
		$gRNA_annotation{$gRNA_tag} 	= $t[4];
		
		if($t[5] =~ /;/){
			$gRNA_annotation{$gRNA_tag} .= ":".$t[5];
		}
		
		if($gRNA_line{$gRNA_tag} =~ /specific:(\S+)/){
			$gRNA_specific{$gRNA_tag} = $1;
		}
		
		if($gRNA_line{$gRNA_tag} =~ /offtargets:(\S+)/){
			$gRNA_offtargets{$gRNA_tag} = $1;
		}
		
		if($gRNA_specific{$gRNA_tag} && $t[3] =~ /([\w\d]+)_(\d+)_(\d+)/){
			$chr = $1;
			$left = $2;
			$right = $3;
			if(defined $geneID{$chr."_".$left} && defined $geneID{$chr."_".$right}){
				
				undef @t3;
				@t3 = split(";", $geneID{$chr."_".$left});
				undef @t4;
				@t4 = split(";", $geneID{$chr."_".$right});
				
				undef %test;
				foreach $j(@t3){
					$test{$j} ++;
				}
				foreach $j(@t4){
					$test{$j} ++;
				}
				
				foreach $j (keys %test){
					if($test{$j} == 2){
						push (@{$gene_gRNA{$j}}, $gRNA_tag);
					}
				}
			}
		}

	}else{
		$gRNA_seq{$gRNA_tag} = $_;
	}
}
close IN;
print "\n";
undef $count;


open (OUT, "> $out") or die $!;
print OUT "Gene_ID\tLocus";
for $i (1..10){
	print OUT "\t".$i."\tgRNA_locus\tSpecific_score\tOff_targets\tgRNA_seq\tPAM_type\tPAM_site\tCleavage_site\tAnnotation";
}
print OUT "\n";

foreach $i (sort keys %gene_gRNA){
	$count ++; 
	print "Outputting lines: $count\r" if $count % 10 == 0;
	
	@{$gene_gRNA{$i}} = sort { $gRNA_specific{$b} <=> $gRNA_specific{$a} } @{$gene_gRNA{$i}};
	
	print OUT $i."\t".$locus_line{$i};
	
	for $j(0..9){
		if(${$gene_gRNA{$i}}[$j]){
			$tag = ${$gene_gRNA{$i}}[$j];
			$k = $j +1;
			print OUT "\t".$k."\t".$tag;
			print OUT "\t".$gRNA_specific{$tag};
			print OUT "\t".$gRNA_offtargets{$tag};	
			print OUT "\t".$gRNA_seq{$tag};
			print OUT "\t".$gRNA_PAM_type{$tag};
			print OUT "\t".$gRNA_PAM_site{$tag};
			print OUT "\t".$gRNA_cleavage_site{$tag};
			print OUT "\t".$gRNA_annotation{$tag};	
		}	
	}
	print OUT "\n";
}
close OUT;
print "\n";


