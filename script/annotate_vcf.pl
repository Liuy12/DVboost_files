#!/usr/bin/perl

# Annotate VCF 
use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;

# script options
my($input, $output, $help, $trio, $MIN_SV_LEN, $RECIP_OVERLAP_FRAC, $MIN_PE_SUPPORT, $MIN_SR_SUPPORT, $BND_BREAKPOINT_BUFFER, $germ_art, $cnvmap_bed, $cnvr_bed, $seg_dup_bed, $refflat, $strict_mask, $bedtools);
# defaults
$MIN_SV_LEN = 50;
$RECIP_OVERLAP_FRAC = 0.5;
$MIN_PE_SUPPORT = 2;
$MIN_SR_SUPPORT = 2;
$BND_BREAKPOINT_BUFFER = 200;
GetOptions("in|i:s" => \$input,
		"out|o:s" => \$output,
		"germ_art|g:s" => \$germ_art,
		"cnvmap_bed|c:s" => \$cnvmap_bed,
		"cnvr_bed|v:s" => \$cnvr_bed,
		"seg_dup_bed|d:s" => \$seg_dup_bed,
		"refflat|f:s" => \$refflat,
		"strict_mask|m:s" => \$strict_mask,
		"bedtools|e:s" => \$bedtools,
		"trio|t:s" => \$trio,
		"min_length|l:i" => \$MIN_SV_LEN,
		"recip|r:i" => \$RECIP_OVERLAP_FRAC,
		"min_pe|p:i" => \$MIN_PE_SUPPORT,
		"min_sr|s:i" => \$MIN_SR_SUPPORT,
		"bnd_buffer|b:i" => \$BND_BREAKPOINT_BUFFER,
		"help|h|?" => \&help); 

if(not defined $input or not defined $output){
	&help();
}

my %trio_samples = ();
if(defined $trio){
	my @trio_arr = split(",",$trio);
	$trio_samples{"child"} = $trio_arr[0];
	$trio_samples{"father"} = $trio_arr[1];
	$trio_samples{"mother"} = $trio_arr[2];
}

open IN, "<$input" or die "opening $input\n";
open OUT, ">$output" or die "opening $output\n";
open GER, "<$germ_art" or die "opening $germ_art\n";
open BED, ">$input.tmp.germline.bed" or die "opening $input.tmp.germline.bed\n";
open REF, "<$refflat" or die "opening $refflat\n";

# load lumpy VCF
my @chrs = ("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y");
my %acceptable_chrs = ();
foreach my $c (@chrs){ $acceptable_chrs{"chr".$c} = 1;}
my @samples = ();
my %sample_index = ();
my @formats_to_use = ("GT","PE","SR");
my $final_header = "ChrA\tStart\tChrB\tEnd\tSVLen\tSVType";
my %additional_sample_info = ();
while(<IN>){
	my $row = $_;
	chomp $row;
	if($row =~ /^#CHROM/){
		my @line = split("\t",$row);
		@samples = @line[9..scalar(@line)-1];
		my $sample_counter = 9;
		foreach my $sample (@samples){
			foreach my $fmt (@formats_to_use){
				$final_header.="\t".$sample."_".$fmt;
			}
			$sample_index{$sample} = $sample_counter;
			$sample_counter++;
		}
	}
	next if $row =~ /^#/;
	my @line = split("\t",$row);
	my $chr1 = $line[0];
	my $start = $line[1];
	# parse INFO
	my $results = $line[7] =~ /SVTYPE=(.*?);/;
	my $svtype = "NA";
	$svtype = $1 if $results;
	my $svlen = "NA";
	$results = $line[7] =~ /;SVLEN=(.*?);/;
	$svlen = $1 if $results;
	my $chr2 = $chr1;
	my $end = "NA";
	if($svtype eq "BND"){
		my @break2 = ();
		$results = $line[4] =~ /\](.*?)\]/;
		@break2 = split(":",$1) if $results;
		$results = $line[4] =~ /\[(.*?)\[/;
		@break2 = split(":",$1) if $results;
		$chr2 = $break2[0];
		$end = $break2[1];
		if($chr1 eq $chr2){
			$svlen = $end-$start;
		}
	}else{
		$results = $line[7] =~ /;END=(.*?);/;
		$end = $1 if $results;
	}
	if($svlen ne "NA"){
		$svlen = $svlen*-1 if $svlen < 0;
		next if $svlen < $MIN_SV_LEN;
        }

	next if !($acceptable_chrs{$chr1} && $acceptable_chrs{$chr2});
	# parse FORMAT
	my @format = split(":",$line[8]);
	my %format_ind = ();
	foreach my $fmt (@formats_to_use){
		for(my $i = 0; $i < scalar(@format); $i++){
			$format_ind{$fmt} = $i if $fmt eq $format[$i];
		}
	}
	my $sample_info = "";
	for(my $i = 0; $i < scalar(@samples); $i++){
		my $sample = $samples[$i];
		my @sample_field = split(":",$line[9+$i]);
		foreach my $fmt (@formats_to_use){
			$sample_info.="\t".$sample_field[$format_ind{$fmt}]
		}
	}
	$additional_sample_info{$chr1.":".$start.":".$chr2.":".$end.":".$svtype} = $chr1."\t".$start."\t".$chr2."\t".$end."\t".$svlen."\t".$svtype.$sample_info;
}

# load germline artifact calls
my @germline_calls = (); #chr1:start:chr2:end:svtype:sample
while(<GER>){
	my $row = $_;
	chomp $row;
	my @line = split("\t",$row);
	my ($chr1,$start,$chr2,$end,$svtype,$sample) = @line;
	next if !($acceptable_chrs{$chr1} && $acceptable_chrs{$chr2});
	push(@germline_calls,join(":",@line));
	if($chr1 eq $chr2 and $svtype ne "BND"){
		if($start <= $end){
			print BED $chr1."\t".$start."\t".$end."\t".$svtype."_".$sample."\n";
		}else{
			print BED $chr1."\t".$end."\t".$start."\t".$svtype."_".$sample."\n";
		}
	}else{
		# give a buffer to each breakpoint of translocations so we can still use intersectBed
		print BED $chr1."\t".($start-$BND_BREAKPOINT_BUFFER)."\t".($start+$BND_BREAKPOINT_BUFFER)."\t".$svtype."_".$sample."\n";
		print BED $chr2."\t".($end-$BND_BREAKPOINT_BUFFER)."\t".($end+$BND_BREAKPOINT_BUFFER)."\t".$svtype."_".$sample."\n";
	}
}
close GER;
close BED;

# create temp output bed
open TMP, ">$input.tmp.bed" or die "opening $input.tmp.bed\n";
open TMP_SM, ">$input.tmp.strict_mask.bed" or die "opening $input.tmp.strict_mask.bed\n";
foreach my $record (keys %additional_sample_info){
	my @rec_arr = split(":",$record);
	if($rec_arr[0] eq $rec_arr[2]){
		if($rec_arr[1] <= $rec_arr[3]){
			print TMP join("\t",$rec_arr[0],$rec_arr[1],$rec_arr[3])."\n";
		}else{
			print TMP join("\t",$rec_arr[0],$rec_arr[3],$rec_arr[1])."\n";
		}
	}else{ # for translocations we only care about annotating the breakpoints
		print TMP join("\t",$rec_arr[0],($rec_arr[1]-$BND_BREAKPOINT_BUFFER),($rec_arr[1]+$BND_BREAKPOINT_BUFFER))."\n";
		print TMP join("\t",$rec_arr[2],($rec_arr[3]-$BND_BREAKPOINT_BUFFER),($rec_arr[3]+$BND_BREAKPOINT_BUFFER))."\n";
	}
	print TMP_SM join("\t",$rec_arr[0],$rec_arr[1]-50,$rec_arr[1]+50,$rec_arr[0]."_".$rec_arr[1])."\n";
	print TMP_SM join("\t",$rec_arr[2],$rec_arr[3]-50,$rec_arr[3]+50,$rec_arr[2]."_".$rec_arr[3])."\n";
}
close TMP;
close TMP_SM;

# create refflat bed files (coding/gene body)
open CODING, ">$input.tmp.coding.bed" or die "opening $input.tmp.coding.bed\n";
open GENE, ">$input.tmp.gene.bed" or die "opening $input.tmp.gene.bed\n";
while(<REF>){
	my $row = $_;
	chomp $row;
	my @line = split("\t",$row);
	print GENE join("\t",$line[2],$line[4],$line[5],$line[0])."\n";
	my @starts = split(",",$line[9]);
	my @stops = split(",",$line[10]);
	for(my $i = 0; $i < scalar(@starts); $i++){
		print CODING join("\t",$line[2],$starts[$i],$stops[$i],$line[0])."\n";
	}
}
close CODING;
close GENE;
close REF;

# Perform intersections with annotations
`cat $input.tmp.bed | sortBed > $input.tmp.sort.bed`;
`cat $input.tmp.strict_mask.bed | sortBed > $input.tmp.strict_mask.sort.bed`;
`$bedtools intersect -a $input.tmp.sort.bed -b $input.tmp.germline.bed -f $RECIP_OVERLAP_FRAC -r -c > $input.tmp.intersect.germ.bed`;
`$bedtools intersect -a $input.tmp.sort.bed -b $cnvmap_bed -f $RECIP_OVERLAP_FRAC -r -c > $input.tmp.intersect.cnvmap.bed`;
`$bedtools intersect -a $input.tmp.sort.bed -b $cnvr_bed -f $RECIP_OVERLAP_FRAC -r -c > $input.tmp.intersect.cnvr.bed`;
`$bedtools intersect -a $input.tmp.sort.bed -b $seg_dup_bed -f $RECIP_OVERLAP_FRAC -r -wo | cut -f1,2,3,7 > $input.tmp.intersect.seg_dup.bed`;
`$bedtools intersect -a $input.tmp.sort.bed -b $input.tmp.gene.bed -wo | cut -f1,2,3,7 > $input.tmp.intersect.gene.bed`;
`$bedtools intersect -a $input.tmp.sort.bed -b $input.tmp.coding.bed -wo | cut -f1,2,3,7 > $input.tmp.intersect.coding.bed`;
`$bedtools getfasta -bed $input.tmp.strict_mask.sort.bed -fi $strict_mask -fo $input.tmp.strict_mask.sort.sequence.txt -tab -name`;

my %intersected = ();
open INTER, "<$input.tmp.intersect.germ.bed" or die "opening $input.tmp.intersect.germ.bed\n";
while(<INTER>){
	my $row = $_;
	chomp $row;
	my @line = split("\t",$row);
	$intersected{$line[0]."_".$line[1]."_".$line[2]} = $line[3];
}
close INTER;

my %intersected_cnvmap = ();
open INTER, "<$input.tmp.intersect.cnvmap.bed" or die "opening $input.tmp.intersect.cnvmap.bed\n";
while(<INTER>){
        my $row = $_;
        chomp $row;
        my @line = split("\t",$row);
        $intersected_cnvmap{$line[0]."_".$line[1]."_".$line[2]} = $line[3];
}
close INTER;

my %intersected_cnvr = ();
open INTER, "<$input.tmp.intersect.cnvr.bed" or die "opening $input.tmp.intersect.cnvr.bed\n";
while(<INTER>){
        my $row = $_;
        chomp $row;
        my @line = split("\t",$row);
        $intersected_cnvr{$line[0]."_".$line[1]."_".$line[2]} = $line[3];
}
close INTER;

my %intersected_seg_dup = ();
open INTER, "<$input.tmp.intersect.seg_dup.bed" or die "opening $input.tmp.intersect.seg_dup.bed\n";
while(<INTER>){
        my $row = $_;
        chomp $row;
        my @line = split("\t",$row);
        if(not exists $intersected_seg_dup{$line[0]."_".$line[1]."_".$line[2]}){
                $intersected_seg_dup{$line[0]."_".$line[1]."_".$line[2]} = $line[3];
        }else{
                $intersected_seg_dup{$line[0]."_".$line[1]."_".$line[2]} .= ",".$line[3];
        }
}
close INTER;


my %intersected_gene = ();
open INTER, "<$input.tmp.intersect.gene.bed" or die "opening $input.tmp.intersect.gene.bed\n";
while(<INTER>){
	my $row = $_;
	chomp $row;
	my @line = split("\t",$row);
	if(not exists $intersected_gene{$line[0]."_".$line[1]."_".$line[2]}){
		$intersected_gene{$line[0]."_".$line[1]."_".$line[2]} = $line[3];
	}else{
		$intersected_gene{$line[0]."_".$line[1]."_".$line[2]} .= ",".$line[3];
	}
}
close INTER;

my %intersected_coding = ();
open INTER, "<$input.tmp.intersect.coding.bed" or die "opening $input.tmp.intersect.coding.bed\n";
while(<INTER>){
	my $row = $_;
	chomp $row;
	my @line = split("\t",$row);
	if(not exists $intersected_coding{$line[0]."_".$line[1]."_".$line[2]}){
		$intersected_coding{$line[0]."_".$line[1]."_".$line[2]} = $line[3];
	}else{
		$intersected_coding{$line[0]."_".$line[1]."_".$line[2]} .= ",".$line[3];
	}
}
close INTER;

my %intersected_strictmask = ();
open INTER, "<$input.tmp.strict_mask.sort.sequence.txt" or die "opening $input.tmp.strict_mask.sort.sequence.txt\n";
while(<INTER>){ 
        my $row = $_;
        chomp $row;
        my @line = split("\t",$row);
        $intersected_strictmask{$line[0]} = $line[1];
}
close INTER;


foreach my $cnv (keys %intersected_gene){
	my @values = split(",",$intersected_gene{$cnv});
	$intersected_gene{$cnv} = join(",",sort(uniq(@values)));
}

foreach my $cnv (keys %intersected_coding){
	my @values = split(",",$intersected_coding{$cnv});
	$intersected_coding{$cnv} = join(",",sort(uniq(@values)));
}


# print output table
print OUT $final_header."\tGERMLINE_ARTIFACT\tPOLYMORPHIC_CNVMAP\tPOLYMORPHIC_CNVR\tINHERITANCE\tCODING\tGENEBODY\tSEG_DUP\tL_START_STRICTMASK\tH_START_STRICTMASK\tZ_START_STRICTMASK\tQ_START_STRICTMASK\tL_END_STRICTMASK\tH_END_STRICTMASK\tZ_END_STRICTMASK\tQ_END_STRICTMASK\n";
foreach my $cnv (sort sort_chr_pos keys %additional_sample_info){
	my $row = $additional_sample_info{$cnv};
	chomp $row;
	my @line = split("\t",$row);

	my $chr1 = $line[0];
	my $start = $line[1];
	my $chr2 = $line[2];
	my $end = $line[3];
	my $svtype = $line[5];

	my $recip_count_germ = 0;
	my $recip_count_cnvmap = 0;
	my $recip_count_cnvr = 0;
	my $gene = ".";
	my $coding = ".";
	my $seg_dup = ".";
	my $strictmask_start_l = ".";
	my $strictmask_start_h = ".";
	my $strictmask_start_z = ".";
	my $strictmask_start_q = ".";
	my $strictmask_end_l = ".";
        my $strictmask_end_h = ".";
        my $strictmask_end_z = ".";
        my $strictmask_end_q = ".";

	# germline artifacts
	if(exists $intersected{$chr1."_".$start."_".$end} and $chr1 eq $chr2){
		$recip_count_germ = $intersected{$chr1."_".$start."_".$end};
	}elsif(exists $intersected{$chr1."_".$end."_".$start} and $chr1 eq $chr2){
		$recip_count_germ = $intersected{$chr1."_".$end."_".$start}
	}elsif($svtype eq "BND"){
		my $new_start_start = $start-$BND_BREAKPOINT_BUFFER;
		my $new_start_end = $start+$BND_BREAKPOINT_BUFFER;
		my $new_end_start = $end-$BND_BREAKPOINT_BUFFER;
		my $new_end_end = $end+$BND_BREAKPOINT_BUFFER;
		if(exists $intersected{$chr1."_".$new_start_start."_".$new_start_end}){
			$recip_count_germ = $intersected{$chr1."_".$new_start_start."_".$new_start_end};
		}
		if(exists $intersected{$chr2."_".$new_end_start."_".$new_end_end}){
			$recip_count_germ += $intersected{$chr2."_".$new_end_start."_".$new_end_end};
		}
	}else{
		$recip_count_germ = 0;
	}

	# polymorphic CNVs
        if(exists $intersected_cnvmap{$chr1."_".$start."_".$end} and $svtype ne "BND" and $svtype ne "INV"){
                $recip_count_cnvmap = $intersected_cnvmap{$chr1."_".$start."_".$end};
        }elsif(exists $intersected_cnvmap{$chr1."_".$end."_".$start} and $svtype ne "BND" and $svtype ne "INV"){
                $recip_count_cnvmap = $intersected_cnvmap{$chr1."_".$end."_".$start}
        }else{
                $recip_count_cnvmap = "NA";
        }

        if(exists $intersected_cnvr{$chr1."_".$start."_".$end} and $svtype ne "BND" and $svtype ne "INV"){
                $recip_count_cnvr = $intersected_cnvr{$chr1."_".$start."_".$end};
        }elsif(exists $intersected_cnvr{$chr1."_".$end."_".$start} and $svtype ne "BND" and $svtype ne "INV"){
                $recip_count_cnvr = $intersected_cnvr{$chr1."_".$end."_".$start}
        }else{
                $recip_count_cnvr = "NA";
        }

	# Gene Body
	if(exists $intersected_gene{$chr1."_".$start."_".$end} and $svtype ne "BND"){
		$gene = $intersected_gene{$chr1."_".$start."_".$end};
	}elsif(exists $intersected_gene{$chr1."_".$end."_".$start} and $svtype ne "BND"){
		$gene = $intersected_gene{$chr1."_".$end."_".$start}
	}elsif($svtype eq "BND"){
		if(exists $intersected_gene{$chr1."_".($start-1)."_".$start}){
			$gene = $intersected_gene{$chr1."_".($start-1)."_".$start};
		}
		if(exists $intersected_gene{$chr1."_".($end-1)."_".$end}){
			if($gene eq "."){
				$gene = $intersected_gene{$chr1."_".($end-1)."_".$end};
			}else{
				$gene = $gene.",".$intersected_gene{$chr1."_".($end-1)."_".$end};
			}
		}
	}else{
		$gene = ".";
	}

	# Coding region
	if(exists $intersected_coding{$chr1."_".$start."_".$end} and $svtype ne "BND"){
		$coding = $intersected_coding{$chr1."_".$start."_".$end};
	}elsif(exists $intersected_coding{$chr1."_".$end."_".$start} and $svtype ne "BND"){
		$coding = $intersected_coding{$chr1."_".$end."_".$start}
	}elsif($svtype eq "BND"){
                if(exists $intersected_coding{$chr1."_".($start-1)."_".$start}){
                        $coding = $intersected_coding{$chr1."_".($start-1)."_".$start};
                }
                if(exists $intersected_coding{$chr1."_".($end-1)."_".$end}){
                        if($coding eq "."){
                                $coding = $intersected_coding{$chr1."_".($end-1)."_".$end};
                        }else{
                                $coding = $coding.",".$intersected_coding{$chr1."_".($end-1)."_".$end};
                        }
                }
	}else{
		$coding = ".";
	}

	# segmental dups
        if(exists $intersected_seg_dup{$chr1."_".$start."_".$end} and $svtype ne "BND"){
                $seg_dup = $intersected_seg_dup{$chr1."_".$start."_".$end};
        }elsif(exists $intersected_seg_dup{$chr1."_".$end."_".$start} and $svtype ne "BND"){
                $seg_dup = $intersected_seg_dup{$chr1."_".$end."_".$start}
	}elsif($svtype eq "BND"){
		if(exists $intersected_seg_dup{$chr1."_".($start-1)."_".$start}){
                        $seg_dup = $intersected_seg_dup{$chr1."_".($start-1)."_".$start};
                }
                if(exists $intersected_seg_dup{$chr1."_".($end-1)."_".$end}){
                        if($seg_dup eq "."){
                                $seg_dup = $intersected_seg_dup{$chr1."_".($end-1)."_".$end};
                        }else{
                                $seg_dup = $seg_dup.",".$intersected_seg_dup{$chr1."_".($end-1)."_".$end};
                        }
                }
        }else{
                $seg_dup = ".";
        }

	# strict mask
	if(exists $intersected_strictmask{$chr1."_".$start}){
		my @sm_l = $intersected_strictmask{$chr1."_".$start} =~ /L/g;
		my @sm_h = $intersected_strictmask{$chr1."_".$start} =~ /H/g;
		my @sm_z = $intersected_strictmask{$chr1."_".$start} =~ /Z/g;
		my @sm_q = $intersected_strictmask{$chr1."_".$start} =~ /Q/g;
		$strictmask_start_l = @sm_l;
		$strictmask_start_h = @sm_h;
		$strictmask_start_z = @sm_z;
		$strictmask_start_q = @sm_q;
	}
	if(exists $intersected_strictmask{$chr2."_".$end}){
                my @sm_l = $intersected_strictmask{$chr2."_".$end} =~ /L/g;
                my @sm_h = $intersected_strictmask{$chr2."_".$end} =~ /H/g;
                my @sm_z = $intersected_strictmask{$chr2."_".$end} =~ /Z/g;
                my @sm_q = $intersected_strictmask{$chr2."_".$end} =~ /Q/g;
                $strictmask_end_l = @sm_l;
                $strictmask_end_h = @sm_h;
                $strictmask_end_z = @sm_z;
                $strictmask_end_q = @sm_q;
        }


	my $inheritance = ".";
	if(defined $trio){
		my $proband_support = 0;
		my $pe_support = 0;
		my $sr_support = 0;
		$pe_support += $line[(($sample_index{$trio_samples{"child"}}-8)*3)+4] if $line[(($sample_index{$trio_samples{"child"}}-8)*3)+4] ne ".";
		$sr_support += $line[(($sample_index{$trio_samples{"child"}}-8)*3)+5] if $line[(($sample_index{$trio_samples{"child"}}-8)*3)+5] ne ".";
		if($pe_support >= $MIN_PE_SUPPORT and $sr_support >= $MIN_SR_SUPPORT){
			$proband_support += $pe_support;
			$proband_support += $sr_support;
		}
		if($proband_support >= 5){
			# figure out inheritence mode
			my $p_gt = $line[(($sample_index{$trio_samples{"child"}}-8)*3)+3];
			my $f_gt = $line[(($sample_index{$trio_samples{"child"}}-8)*3)+3];
			my $m_gt = $line[(($sample_index{$trio_samples{"child"}}-8)*3)+3];
			my $p_support = 0;
			my $f_support = 0;
			my $m_support = 0;
			$p_support += $line[(($sample_index{$trio_samples{"child"}}-8)*3)+4] if $line[(($sample_index{$trio_samples{"child"}}-8)*3)+4] ne ".";
			$p_support += $line[(($sample_index{$trio_samples{"child"}}-8)*3)+5] if $line[(($sample_index{$trio_samples{"child"}}-8)*3)+5] ne ".";
			$f_support += $line[(($sample_index{$trio_samples{"father"}}-8)*3)+4] if $line[(($sample_index{$trio_samples{"father"}}-8)*3)+4] ne ".";
			$f_support += $line[(($sample_index{$trio_samples{"father"}}-8)*3)+5] if $line[(($sample_index{$trio_samples{"father"}}-8)*3)+5] ne ".";
			$m_support += $line[(($sample_index{$trio_samples{"mother"}}-8)*3)+4] if $line[(($sample_index{$trio_samples{"mother"}}-8)*3)+4] ne ".";
			$m_support += $line[(($sample_index{$trio_samples{"mother"}}-8)*3)+5] if $line[(($sample_index{$trio_samples{"mother"}}-8)*3)+5] ne ".";
			# see if svtyper was run
			if($p_gt eq "0/0" or $p_gt eq "0/1" or $p_gt eq "1/1"){
				if($f_gt eq "0/0" and $m_gt eq "0/0" and $p_gt ne "0/0"){
					$inheritance = "denovo";
				}elsif($p_gt eq "1/1" and $f_gt eq "0/1" and $m_gt eq "0/1"){
					$inheritance = "autosomal-recessive";
				}elsif($p_gt eq "1/1" and (($f_gt eq "1/1" and $m_gt eq "0/0") or ($f_gt eq "0/0" and $m_gt eq "1/1"))){
					$inheritance = "autosomal-dominant";
				}elsif($p_gt eq "0/1" and (($f_gt eq "0/1" and $m_gt eq "0/0") or ($f_gt eq "0/0" and $m_gt eq "0/1") or ($f_gt eq "0/1" and $m_gt eq "0/1"))){
					$inheritance = "autosomal-dominant";
				}else{
					$inheritance = ".";
				}
			}else{
				# figure out basic inheritence without svtyper genotypes
				my $parent_support = $f_support + $m_support;
				if($proband_support >= 5 and $parent_support <= 2){
					$inheritance = "denovo";
				}
			}
		}
	}
	print OUT join("\t",$row,$recip_count_germ,$recip_count_cnvmap,$recip_count_cnvr,$inheritance,$coding,$gene,$seg_dup,$strictmask_start_l,$strictmask_start_h,$strictmask_start_z,$strictmask_start_q,$strictmask_end_l,$strictmask_end_h,$strictmask_end_z,$strictmask_end_q)."\n";

}

close IN;
close OUT;

# remove tmp files from bedtools
`rm $input.tmp.intersect.germ.bed $input.tmp.intersect.cnvmap.bed $input.tmp.intersect.cnvr.bed $input.tmp.intersect.gene.bed $input.tmp.intersect.coding.bed $input.tmp.bed $input.tmp.strict_mask.bed $input.tmp.intersect.seg_dup.bed $input.tmp.strict_mask.sort.bed $input.tmp.strict_mask.sort.sequence.txt $input.tmp.coding.bed $input.tmp.gene.bed $input.tmp.sort.bed $input.tmp.germline.bed`;


sub sort_chr_pos {
	my @a_cnv_key = split(":",$a);
	my @b_cnv_key = split(":",$b);
	my %chr_order = (1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13,14,14,15,15,16,16,17,17,18,18,19,19,20,20,21,21,22,22,"X",23,"Y",24);
	$chr_order{substr($a_cnv_key[0],3)} <=> $chr_order{substr($b_cnv_key[0],3)} or $a_cnv_key[1] <=> $b_cnv_key[1];
}


sub uniq{
	my %seen;
	grep !$seen{$_}++, @_;
}

sub help{

        print '
DESCRIPTION:
        annotate_lumpy_vcf.pl will add various annotations to lumpy SV calls and output the
        results in a table format for easy filtering and interpretation.

USAGE:
        annotate_vcf.pl -i path/to/input.vcf -o path/to/output_table.txt -g path/to/germline_artifacts.txt -c path/to/CNVMAP.bed -v path/to/CNVR.bed -d path/to/segmental_dups.bed -f path/to/Homo_sapiens.refFlat -m path/to/strictmask.fa -e path/to/bedtools

OPTIONS:
        --in,-i            Required path to uncompressed VCF file 

        --out,-o           Required path to output annotate file
		
		--germ_art,-g	   Required path to germline_artifacts.txt
		
		--cnvmap_bed,-c	   Required path to CNVMAP.bed
		
		--cnvr_bed,-v	   Required path to CNVR.bed
		
		--seg_dup_bed,-d   Required path to segmental_dups.bed
		
		--refflat,-f	   Required path to Homo_sapiens.refFlat
		
		--strict_mask,-m   Required path to strictmask.fa
		
		--bedtools,-e	   Required path to bedtools

        --trio,-t          Optional parameter to define trio samples in the format: child,father,mother
                           If defined then an inheritance mode annotation will be added

        --min_length,-l    Optional value defining the min SV length to output. Default value: 50

        --recip,-r         Optional value defining the fraction of reciprical overlap required when
                           intersecting SV calls with annotations. Default value: 0.5

        --min_pe,-p        Optional value defining the min number of paired-end support for an SV call
                           Default value: 2

        --min_sr,-s        Optional value defining the min number of split-read support for an SV call
                           Default value: 2

        --bnd_buffer,-b    Optional value defining the size of the buffer to put around BND (translocation)
                           breakpoints when intersecting with annotations. Default value: 200
						   
        --help,-h,-?       Display this help documentation.

';
exit;
}


