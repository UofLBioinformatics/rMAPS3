#!/usr/bin/perl -w
use strict;



##############################################################
sub parseIsoforms { 
  my($iso) = @_;
  my @wds = split(/\@/, $iso);
  my $numWds = @wds;

#  if($numWds != 11) { 
#     die("Error in parseIsoforms: 11 fields expected; $numWds encountered\n");
#  }
  my $first = $wds[0];
  my $second = $wds[1];
  my $third = $wds[2];
  my $fourth = $wds[3];
  

  
  my @wds2 = split(/\:/, $first);
  my $firstS = $wds2[1];
  my $firstE = $wds2[2];
  
  @wds2 = split(/\:/, $second);
  my $secondS = $wds2[1];
  my $secondE = $wds2[2];
  
  @wds2 = split(/\:/, $third);
  my $thirdS = $wds2[1];
  my $thirdE = $wds2[2];
  
  @wds2 = split(/\:/, $fourth);
  my $fourthS = $wds2[1];
  my $fourthE = $wds2[2];

 
  

  my @retArr = ($firstS, $firstE, $secondS, $secondE, 
                $thirdS, $thirdE,$fourthS, $fourthE);
  return(@retArr);
}
##############################################################


##############################################################
sub parseAssignedCounts {
   my ($ac) = @_;

   my $IC = 0;
   my $SC = 0;

   my @wds = split(/\,/, $ac);
   my $numWds = @wds;
   for(my $i = 0; $i < $numWds; $i++) { 
      my @wds2 = split(/\:/, $wds[$i]);
      my $currIso = $wds2[0];
      my $currIsoCnt = $wds2[1];
      if($currIso == 0) { 
         $IC = $currIsoCnt; 
      }
      else {
         if($currIso == 1) { 
            $SC = $currIsoCnt;
         }
         else {
            die("ERROR in parseAssignedCounts: Expecting 0 or 1 for isoform; received $currIso\n");
         }
      }
   }

   my @retArr = ($IC, $SC);
   return(@retArr);
}
##############################################################

###################
## MAIN FUNCTION ##
###################
my $argc = @ARGV;
if($argc != 4) { 
   die("USAGE: perl miso2rMATS_v2.pl <bayes cutoff low> <bayes cutoff high> <inFile> <outfile>\n");
}

my $bayesLOW = $ARGV[0];
my $bayesHIGH = $ARGV[1];
my $inFN = $ARGV[2];
my $outFN = ">" . $ARGV[3];

open(INFILE, $inFN) || die ("Cannot open $inFN for reading");
open(OUTFILE, $outFN) || die("Cannot open $outFN for writing");

my $hdrLine = <INFILE>; ## ASSUME MISO INPUT HAS HEADER
my @wds = split(/\t/, $hdrLine);
my $numCategories = @wds;
if($numCategories != 18) { 
   die("ERROR IN LINE: $hdrLine: theere are $numCategories columns when 18 are expected");
}

print OUTFILE "ID\tGeneID\tgeneSymbol\tchr\tstrand\t1stExonStart_0base\t1stExonEnd\t";
print OUTFILE "2ndExonStart_0base\t2ndExonEnd\tupstreamES\tupstreamEE\tdownstreamES\tdownstreamEE\tIC_SAMPLE_1\t";
print OUTFILE "SC_SAMPLE_1\tIC_SAMPLE_2\tSC_SAMPLE_2\tIncFormLen\tSkipFormLen\t";
print OUTFILE "PValue\tFDR\tIncLevel1\tIncLevel2\tIncLevelDifference\n";

while(defined(my $line = <INFILE>)) { 
   chomp($line);
   my @tmpWds = split(/\t/, $line);
   my $numCats = @tmpWds;
   if($numCats != 18) { 
      die("ERROR IN LINE: $line: theere are $numCats columns when 18 are expected");
   }

   my ($event_name, $sample1_posterior_mean, $sample1_ci_low, $sample1_ci_high,
       $sample2_posterior_mean, $sample2_ci_low, $sample2_ci_high, $diff, $bayes_factor,
       $isoforms, $sample1_counts, $sample1_assigned_counts, $sample2_counts, 
       $sample2_assigned_counts, $chrom,
       $strand, $mRNA_starts, $mRNA_ends) = split(/\t/, $line);

   my $rMATS_chr    = $chrom;
   my $rMATS_strand = $strand;
   my $rMATS_psi1   = $sample1_posterior_mean;
   my $rMATS_psi2   = $sample2_posterior_mean;
   my $rMATS_deltaPSI = $diff;
   my $rMATS_ID         = "NA";
   my $rMATS_GeneID     = "NA";
   my $rMATS_geneSymbol = "NA";
   my $rMATS_firstS = 0;
   my $rMATS_firstE         = 0;
   my $rMATS_secondS      = 0;
   my $rMATS_secondE      = 0;
   my $rMATS_thirdS    = 0;
   my $rMATS_thirdE    = 0; 
   my $rMATS_fourthS    = 0;
   my $rMATS_fourthE    = 0;
   my $rMATS_IC_SAMPLE_1 = 0;
   my $rMATS_SC_SAMPLE_1 = 0;
   my $rMATS_IC_SAMPLE_2 = 0;
   my $rMATS_SC_SAMPLE_2 = 0;  
   my $rMATS_IncFormLen = 0;
   my $rMATS_SkipFormLen = 0;
   my $rMATS_PValue = 1;
   my $rMATS_FDR    = 1;
   if($bayes_factor >= $bayesHIGH) { 
      $rMATS_FDR = 0.01;
   }
   else {
      if($bayes_factor <= $bayesLOW) { 
         $rMATS_FDR = 1;
      }
      else {
         $rMATS_FDR = 0.4;
      }
   }
   $rMATS_PValue = 10**(-1*$bayes_factor);

   ($rMATS_IC_SAMPLE_1, $rMATS_SC_SAMPLE_1) = parseAssignedCounts($sample1_assigned_counts);
   ($rMATS_IC_SAMPLE_2, $rMATS_SC_SAMPLE_2) = parseAssignedCounts($sample2_assigned_counts);

   ($rMATS_firstS, $rMATS_firstE, $rMATS_secondS,$rMATS_secondE,
    $rMATS_thirdS, $rMATS_thirdE,$rMATS_fourthS, $rMATS_fourthE) = parseIsoforms($isoforms, $strand);

   print OUTFILE "$rMATS_ID\t$rMATS_GeneID\t$rMATS_geneSymbol\t$rMATS_chr\t$rMATS_strand\t";

   if ($strand eq "+"){
      $rMATS_secondS=$rMATS_secondS-1;
	  $rMATS_thirdS=$rMATS_thirdS-1;
      $rMATS_firstS=$rMATS_firstS-1;
	  $rMATS_fourthS=$rMATS_fourthS-1;
      print OUTFILE "$rMATS_secondS\t$rMATS_secondE\t";
      print OUTFILE "$rMATS_thirdS\t$rMATS_thirdE\t$rMATS_firstS\t$rMATS_firstE\t$rMATS_fourthS\t$rMATS_fourthE\t";
   }
   if ($strand eq "-"){
      $rMATS_thirdS=$rMATS_thirdS-1;
	  $rMATS_secondS=$rMATS_secondS-1;
	  $rMATS_fourthS=$rMATS_fourthS-1;
	  $rMATS_firstS=$rMATS_firstS-1;
      print OUTFILE "$rMATS_thirdS\t$rMATS_thirdE\t";
      print OUTFILE "$rMATS_secondS\t$rMATS_secondE\t$rMATS_fourthS\t$rMATS_fourthE\t$rMATS_firstS\t$rMATS_firstE\t";
   }
   print OUTFILE "$rMATS_ID\t$rMATS_IC_SAMPLE_1\t$rMATS_SC_SAMPLE_1\t$rMATS_IC_SAMPLE_2\t";
   print OUTFILE "$rMATS_SC_SAMPLE_2\t$rMATS_IncFormLen\t$rMATS_SkipFormLen\t$rMATS_PValue\t";
   print OUTFILE "$rMATS_FDR\t$rMATS_psi1\t$rMATS_psi2\t$rMATS_deltaPSI\n";
}

close(INFILE);
close(OUTFILE);
