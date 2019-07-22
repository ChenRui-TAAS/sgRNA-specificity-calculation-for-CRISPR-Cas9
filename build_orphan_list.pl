#! usr/bin/perl -w
use strict;

# Author:			Rui Chen <chenrui.taas@gmail.com; chenrui.2011@outlook.com>
# Program Date:		2014.07.30		
# Modified:			2015.03.09 
############################################################

my $in = shift;
my $out = shift;

die "
Usage:	
perl   build_orphan_list.pl    input_sgRNA.fa    output_filename

For example:
input_sgRNA.fa	sgRNA_whole_genome_TAIR10.fa
output_filename	sgRNA_whole_genome_TAIR10_orphan.fa

" if !defined $out;

my(
$count,
$i,
@t,
%orphan,
%name,
$count_out
);
############################################################


## Record fasta reads
open (IN, "< $in") or die $!;
while (<IN>) {
    $_ =~ s/[\r\n]+//g;
	if($_ =~ /^>(\S+)/){
		$count ++;
		print "Recording reads: $count\r" if $count % 10000 == 0;
	    
		$i = 0;
	}else{
		$i ++;
	}
	
	$t[$i] = $_;
	
	if($i == 1){
		$orphan{$t[1]} ++;
		$name{$t[1]} .= $t[0];
	}
}
close IN;
print "\n";


open (IN, "< $in") or die $!;
open (OUT, "> $out") or die $!;
while (<IN>) {
    $_ =~ s/[\r\n]+//g;
	if($_ =~ /^>(\S+)/){
		$i = 0;
	}else{
		$i ++;
	}
	
	$t[$i] = $_;
	
	if($i == 1){
		if($orphan{$t[1]} == 1){
			
			$count_out ++;
			print "Outputting orphans: $count_out\r" if $count_out % 10000 == 0;
	
			print OUT $name{$t[1]}."\n".$t[1]."\n";
		}
	}
}
close IN;
close OUT;

print "\n\nTotal input fasta:\t$count\n";
print "Total output orphan fasta:\t$count_out\n";

