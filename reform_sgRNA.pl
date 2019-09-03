#! usr/bin/perl -w
use strict;

# Author:		Rui Chen <chenrui.taas@gmail.com; chenrui.2011@outlook.com>
# Program Date:		2017.02.20	
# Version:		1.0
		
############################################################

my $in = shift;
my $out = shift;

die "
Usage: 
perl   reform_sgRNA.pl   orphan.fa   output_filename

For example:
orphan.fa	sgRNA_whole_genome_TAIR10_orphan.fa
output_filename	sgRNA_whole_genome_TAIR10_orphan_23nt.fa

" if !defined $out;

my %tr = ( 
"A" => "T",
"T" => "A",
"C" => "G",
"G" => "C",
"N" => "N"
);

my(
@t, 
$i, 
$PAM_tag, 
$PAM_seq, 
$out_seq, 
$guide_ID, 
$count
); 

############################################################


open (IN, "< $in") or die $!;
open (OUT, "> $out") or die $!;
while (<IN>) {
	chomp;
	
	if($_ =~ /^>(\S+)/){
		$count ++; 
		print "Recording and outputting fasta: $count\r" if $count % 1000 == 0;
		
		$i = 0;
		$guide_ID = $1;
		
		undef @t;
		@t =split("\t", $_);
		$PAM_tag = $t[1];
		$PAM_seq = $t[2];
		
		#print OUT $_."\n";
		print OUT ">".$guide_ID."\n";
		
	}else{
		$i ++;
	}
	
	if($i == 1){
		if($PAM_tag eq "CCN"){
			$out_seq = base_reverse ($PAM_seq);
		}elsif($PAM_tag eq "NGG"){
			$out_seq = $PAM_seq;
		}
		
		print OUT $_.$out_seq."\n";
	}
}
close IN;
close OUT;
print "\n";


sub base_reverse{
	my @bre = reverse(split ('', $_[0]));
	my $baseout;
	foreach $i(@bre){
	    $baseout .= $tr{$i};
	}
	$baseout;
}
