#! usr/bin/perl -w
use strict;

# Author:		Rui Chen <chenrui.taas@gmail.com; chenrui.2011@outlook.com>
# Program Date:		2017.04.13	
# Version:		1.0		
############################################################

my $guide_fa = shift;
my $scoreTable = shift;
my $out = shift;

die "
Usage:	
perl   label_back_sgRNA.pl   sgRNA_orphan.fa   score_table   output_filename

For example:
sgRNA_orphan.fa	sgRNA_whole_genome_TAIR10_orphan.fa
score_table	scoringOut_TAIR10_orphan_23nt.txt
output_filename	sgRNA_whole_genome_TAIR10_orphan_scoreLabeled.fa

" if !defined $out;

`date` =~ /(\d+).+\s+(\d+).+\s+(\d+).+\s+(\d+):(\d+):(\d+)\s+/;
my $date = $1.$2.$3.$4.$5;

my (@t, $tag, $count);
my (%specific_score, %offtarget_number);
############################################################


open (IN, "< $scoreTable") or die $!;
while (<IN>) {
	$_ =~ s/[\r\n]+//g;
	
	undef @t;
	@t =split("\t", $_);
	
	if($t[0] =~ /\d+_[+-]_\d+_\d+/){
		$count ++; 
		print "Recording scoreTable: $count\r" if $count % 1000 == 0;
		
		$specific_score{$t[0]} = $t[1];
		if($t[2] eq "0"){
			$offtarget_number{$t[0]} = "None";
		}else{
			$offtarget_number{$t[0]} = $t[2];
		}
	}
}
close IN;
print "\n";
undef $count;


open (IN, "< $guide_fa") or die $!;
open (OUT, "> $out") or die $!;
while (<IN>) {
	$_ =~ s/[\r\n]+//g;

	if($_ =~ /^>(\S+)/){
		$count ++; 
		print "Outputting fasta: $count\r" if $count % 1000 == 0;
		
		$tag = $1;
		
		if($specific_score{$tag}){
			print OUT $_."\tUniq\tspecific:".$specific_score{$tag}."\tofftargets:".$offtarget_number{$tag}."\n";
		}else{
			print OUT $_."\tNot_uniq\n";
		}

	}else{
		print OUT $_."\n";
	}
}
close IN;
close OUT;
print "\n";
