#! usr/bin/perl -w
use strict;

# Author:			Rui Chen <chenrui.taas@gmail.com; chenrui.2011@outlook.com>
# Program Date:		2017.02.20	specificity scoring system for each sgRNA
# Modified:			2017.03.02	
# Modified:			2017.03.11	
# Modified:			2017.03.13	
# Modified:			2017.03.20  a.mismatch_number > 10000 
#								b.mismatch_number = 0
#								c.23nt_perfect_match and 20nt_perfect_match, output "0" and "error"

## Explanation for gRNA seeds and cleavage site:
## Finding CCN and NGG from 5' to 3'
## Cleavage site: 						(NNNNNNNNNNNNNNNNNxNNN)NGG		CCN(NNNxNNNNNNNNNNNNNNNNN)
## Extracted gRNA seed region:	20bp	5'-(N20)NGG-3'					5'-CCN(N20 reverse complemented)-3'

## Scoring principles:
## more matched number -> more off-targets -> more close to 0 (Poor specificity)
## d value ranging from 20 -> -2, PAM-proximal = -2,  PAM-distal = 20.
## Reference: "CRISPR-P a web tool for synthetic single-guide RNA design of CRISPR-system in plants"

## For example:
## ====================NGG
## Four conditions: 
## ====================NGG	fully match + NGG, discarded  
## ====*====*=*=*======NGG	<=4mismatch + NGG, scoring [0~100]
## ====================NAG	fully match + NAG, specificity score = 0.1 and labeling "NAG_match"
## ====*====*=*=*======NAG	<=4mismatch + NAG, scoring [0~100]
##############################################################

my $batmanOut = shift;
my $out = shift;

if(-e $out){

}else{
	system ("touch $out");
}

die "
Usage:	
perl   scoring_off_targets_and_specificities.pl   batmanOut.txt.gz   scoringOut.txt

For example:
batmanOut.txt.gz	batman_out_TAIR10_orphan_23nt.txt.gz
scoringOut.txt		scoringOut_TAIR10_orphan_23nt.txt

" if !defined $out;

#`date` =~ /(\d+).+\s+(\d+).+\s+(\d+).+\s+(\d+):(\d+):(\d+)\s+/;
#my $date = $1.$2.$3.$4.$5;

##############################################################

my @effect_value =(
0, 
0, 
0.014, 
0, 
0, 
0.395, 
0.317, 
0, 
0.389, 
0.079, 
0.445, 
0.508, 
0.613, 
0.851, 
0.732, 
0.828, 
0.615, 
0.804, 
0.685, 
0.583
);

@effect_value = reverse @effect_value; 

my %relocate = (	## d value <-> effect_value
"20" => "-2",
"19" => "-1",
"18" => "0",
"17" => "1",
"16" => "2",
"15" => "3",
"14" => "4",
"13" => "5",
"12" => "6",
"11" => "7",
"10" => "8",
"9" => "9",
"8" => "10",
"7" => "11",
"6" => "12",
"5" => "13",
"4" => "14",
"3" => "15",
"2" => "16",
"1" => "17",
"0" => "18",
"-1" => "19",
"-2" => "20",
);

my ($count, @t, @temp, $i, $j, $w, $ww, $x, $y, $z, $N_mm, $max); 
my ($ori, $MM_tag, $base, @MM_position, @MM_base, $PAM_judge);
my ($MM_position_sum, $MM_position_average, $multiple_multiplication, %off_target_score_sum);
my (%guide_ID, %off_target_score, %guide_specificity);
my ($NGG_perfect_20nt_judge, $NAG_perfect_20nt_judge, %NGG_perfect_23nt_judge);
my (%perfect_match20nt_NGG, %perfect_match20nt_NAG, %perfect_match23nt);
my ($tag_judge, %cut_judge, %already_done);
my (%count_tags, $count_lines);
############################################################


open(IN, "gunzip -c $batmanOut |") or die $!;
open(OUT, "> $out") or die $!;
while (<IN>){
    $_ =~ s/[\r\n]+//g;
	
	if($_ !~ /^\@SQ/){
		
		undef @t;
		@t = split("\t", $_);
		
		$count ++; print "Inputting and outputting lines:\t$count\r" if $count % 1000 == 0;
		
		next if $cut_judge{$t[0]};
		
		if($t[1] == "0"){				## orientation
			$ori = "+";
		}elsif($t[1] == "16"){
			$ori = "-";
		}else{
			undef $ori;
		}
		
		if($t[12] =~ /MD:Z:([\d\w]+)/){	## mismatch_position
			$MM_tag = $1;
			
			undef @MM_position;
			undef @MM_base;
			undef $w;
			undef $ww;
			undef $base;
			undef $x;
			undef $y;
			$z = 20;

			while(1){
				if($MM_tag =~/^(\d+[ATCG])/){
					$x = $1;
					$MM_tag = $';
					
					$w ++;
					if($w == 1){
						$ww = 0;
					}elsif($w > 1){
						$ww = 1;
					}
					
					if($x =~ /(\d+)([ATCG])/){
						$y = $1;
						$base = $2;
						$z = $z - $y - $ww;		## mismatch site, 
						
						push @MM_position, $z;
						push @MM_base, $base;
					}
				}else{
					last;
				}
			}
			
			undef @temp;
			if($ori && $ori eq "-"){			##relocate position for "-" match
				@temp = @MM_position;
				undef @MM_position;
				for $i (0..$#temp){
					unshift @MM_position, $relocate{$temp[$i]};
				}
				@MM_base = reverse @MM_base;
			}
			
			undef $PAM_judge;
			undef $NGG_perfect_20nt_judge;
			undef $NAG_perfect_20nt_judge;
			
			if(@MM_position){
				if ($MM_position[$#MM_position] ne "-2"){								
					if($MM_position[$#MM_position] ne "-1"){									## PAM site is "NGG"
						pop @MM_position if $MM_position[$#MM_position] eq "0";					## only for 20 to 1 position
						
						if(!@MM_position){

							$NGG_perfect_20nt_judge ++;
							#print "Error! Same 20nt guide was found in line:\n".$_."\n";
						
						}elsif(scalar @MM_position <= 4){	
							$PAM_judge ++;
						}
					}elsif($MM_position[$#MM_position] eq "-1" && $MM_base[$#MM_base] eq "A"){	## PAM site is "NAG"
						pop @MM_position if $MM_position[$#MM_position] eq "-1";
						if(@MM_position){
							pop @MM_position if $MM_position[$#MM_position] eq "0";				## only for 20 to 1 position
						}
						
						if(!@MM_position){
							$NAG_perfect_20nt_judge ++;
						}elsif(scalar @MM_position <= 4){	
							$PAM_judge ++;
						}
					}
				}
				
			}elsif(!@MM_position && $t[12] eq "MD:Z:23"){		## 23nt perfect match number
	
				$NGG_perfect_23nt_judge{$t[0]} ++;
				
			}
			
			##judge and output
			if    (defined $NGG_perfect_23nt_judge{$t[0]} && $NGG_perfect_23nt_judge{$t[0]} > 1){
				
				$perfect_match23nt{$t[0]} ++;		##23nt_perfect
			
			}elsif($NGG_perfect_20nt_judge){
				
				$perfect_match20nt_NGG{$t[0]} ++;	##'N'GG, 20nt_perfect
				
			}elsif($NAG_perfect_20nt_judge){
				
				$guide_ID{$t[0]} ++;
				if($guide_ID{$t[0]} <= 10000){
					$off_target_score{$t[0]."_".$guide_ID{$t[0]}} = 100; ## NAG_perfect_match
				}
				
			}elsif($PAM_judge){
				
				$guide_ID{$t[0]} ++;
				if($guide_ID{$t[0]} <= 10000){		## Maximal 10000 off-targets (<= 6 mismatch for each)
				
					undef $MM_position_sum;
					undef $MM_position_average;
					$multiple_multiplication = 1;

					for $i(0..$#MM_position){
						$MM_position[$i] = 19 if $MM_position[$i] == 20;
						$MM_position_sum += $MM_position[$i];

						$j = $MM_position[$i] -1;
						$multiple_multiplication = $multiple_multiplication * (1 - $effect_value[$j]);			
					}
					
					$N_mm = scalar @MM_position;
					$MM_position_average = $MM_position_sum/$N_mm;
					$off_target_score{$t[0]."_".$guide_ID{$t[0]}} = $multiple_multiplication / (((19 - $MM_position_average)/19*4 +1) * $N_mm * $N_mm) *100; ## [0-100]
				}
			}
			
		}else{
			undef $MM_tag;
			#print "Error in column 13:\t No mismatch_tag found!\n";
		}
		
		
		if(!$tag_judge){
			$tag_judge = $t[0];
		}else{
			
			if($tag_judge eq $t[0] && $guide_ID{$tag_judge} && $guide_ID{$tag_judge} == 10000){  
				$cut_judge{$tag_judge} ++;
				
				if($perfect_match23nt{$tag_judge}){
					print OUT $tag_judge."\t0\t23nt_perfect\n";	## 23nt_perfect_match =0
					#$tag_judge = $t[0];
				}elsif($perfect_match20nt_NGG{$tag_judge}){
					print OUT $tag_judge."\t0\tNGG_perfect\n";	## NGG_perfect_match =0
					#$tag_judge = $t[0];
				
				}else{
				
					for $j(1..$guide_ID{$tag_judge}){
						$off_target_score_sum{$tag_judge} += $off_target_score{$tag_judge."_".$j};
					}
	
					$guide_specificity{$tag_judge} = 100 /(100 + $off_target_score_sum{$tag_judge});	##0-100
					$guide_specificity{$tag_judge} = int ($guide_specificity{$tag_judge} *10000) /100;
					print OUT $tag_judge."\t".$guide_specificity{$tag_judge}."\t".$guide_ID{$tag_judge}."\n";
				}
				
			
			}elsif($tag_judge ne $t[0]){	## 
				
				if($cut_judge{$tag_judge}){	## 
					$tag_judge = $t[0];
					
				}elsif($perfect_match23nt{$tag_judge}){
					print OUT $tag_judge."\t0\t23nt_perfect\n";	## 23nt_perfect_match =0
					$tag_judge = $t[0];
				}elsif($perfect_match20nt_NGG{$tag_judge}){
					print OUT $tag_judge."\t0\tNGG_perfect\n";	## NGG_perfect_match =0
					$tag_judge = $t[0];
				
				}else{
					
					if(!defined $guide_ID{$tag_judge}){
						print OUT $tag_judge."\t100\t0\n";	## No mismatch
						$tag_judge = $t[0];
					}else{
				
						for $j(1..$guide_ID{$tag_judge}){
							$off_target_score_sum{$tag_judge} += $off_target_score{$tag_judge."_".$j};
						}
	
						$guide_specificity{$tag_judge} = 100 /(100 + $off_target_score_sum{$tag_judge});	## 0-100
						$guide_specificity{$tag_judge} = int ($guide_specificity{$tag_judge} *10000) /100;
						print OUT $tag_judge."\t".$guide_specificity{$tag_judge}."\t".$guide_ID{$tag_judge}."\n";
				
						$tag_judge = $t[0];
					}
				}
			}
		}
	}
}


## last line
if($cut_judge{$tag_judge}){	

}elsif($perfect_match23nt{$tag_judge}){
	print OUT $tag_judge."\t0\t23nt_perfect\n";	## 23nt_perfect_match =0
}elsif($perfect_match20nt_NGG{$tag_judge}){
	print OUT $tag_judge."\t0\tNGG_perfect\n";	## NGG_perfect_match =0	
}else{
	if(!defined $guide_ID{$tag_judge}){
		print OUT $tag_judge."\t100\t0\n";	## No mismatch
	}else{
		
		for $j(1..$guide_ID{$tag_judge}){
			$off_target_score_sum{$tag_judge} += $off_target_score{$tag_judge."_".$j};
		}

		$guide_specificity{$tag_judge} = 100 /(100 + $off_target_score_sum{$tag_judge});	##0-100
		$guide_specificity{$tag_judge} = int ($guide_specificity{$tag_judge} *10000) /100;
		print OUT $tag_judge."\t".$guide_specificity{$tag_judge}."\t".$guide_ID{$tag_judge}."\n";
	}			
}
close IN;
close OUT;
print "\n";

