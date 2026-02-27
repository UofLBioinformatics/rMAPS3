#!/usr/bin/perl -w
use strict;



##############################################################
sub parseIsoforms { 
  my($iso) = @_;
  my @wds = split(/\@/, $iso);
  my $numWds = @wds;

  if($numWds != 11) { 
     die("Error in parseIsoforms: 11 fields expected; $numWds encountered\n");
  }
  my $upstream = $wds[0];
  my $skipped = $wds[1];
  my $downstream = $wds[2];
  
  my @wds2 = split(/\:/, $upstream);
  my $strand = $wds2[3];
  my $upstreamStart = $wds2[1];
  my $upstreamEnd   = $wds2[2];
 
  @wds2 = split(/\:/, $skipped);
  my $skippedStart = $wds2[1];
  my $skippedEnd   = $wds2[2];

  @wds2 = split(/\:/, $downstream);

  my $downstreamStart = $wds2[1];
  my $downstreamEnd   = $wds2[2];

  my @retArr = ($upstreamStart, $upstreamEnd, $skippedStart, $skippedEnd, 
                $downstreamStart, $downstreamEnd);
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

print OUTFILE "ID\tGeneID\tgeneSymbol\tchr\tstrand\texonStart_0base\texonEnd\t";
print OUTFILE "upstreamES\tupstreamEE\tdownstreamES\tdownstreamEE\tID\tIC_SAMPLE_1\t";
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
   my $rMATS_exonStart_0base = 0;
   my $rMATS_exonEnd         = 0;
   my $rMATS_upstreamES      = 0;
   my $rMATS_upstreamEE      = 0;
   my $rMATS_downstreamES    = 0;
   my $rMATS_downstreamEE    = 0; 
   my $rMATS_IC_SAMPLE_1 = 0;
   my $rMATS_SC_SAMPLE_1 = 0;
   my $rMATS_IC_SAMPLE_2 = 0;
   my $rMATS_SC_SAMPLE_2 = 0;  
   my $rMATS_IncFormLen = 0;
   my $rMATS_SkipFormLen = 0;
   my $rMATS_PValue = 1;
   my $rMATS_FDR    = 1;
   if($bayes_factor >= $bayesHIGH) { 
      $rMATS_PValue = 0.01;
   }
   else {
      if($bayes_factor <= $bayesLOW) { 
         $rMATS_PValue = 1;
      }
      else {
         $rMATS_PValue = 0.4;
      }
   }

   $rMATS_FDR = $rMATS_PValue;
   ($rMATS_IC_SAMPLE_1, $rMATS_SC_SAMPLE_1) = parseAssignedCounts($sample1_assigned_counts);
   ($rMATS_IC_SAMPLE_2, $rMATS_SC_SAMPLE_2) = parseAssignedCounts($sample2_assigned_counts);

   ($rMATS_exonStart_0base, $rMATS_exonEnd, $rMATS_upstreamES,$rMATS_upstreamEE,
    $rMATS_downstreamES, $rMATS_downstreamEE) = parseIsoforms($isoforms, $strand);

   print OUTFILE "$rMATS_ID\t$rMATS_GeneID\t$rMATS_geneSymbol\t$rMATS_chr\t$rMATS_strand\t";

   if ($strand eq "+"){
      print OUTFILE "$rMATS_upstreamES\t$rMATS_upstreamEE\t";
      print OUTFILE "$rMATS_exonStart_0base\t$rMATS_exonEnd\t$rMATS_downstreamES\t$rMATS_downstreamEE\t";
   }
   if ($strand eq "-"){
      print OUTFILE "$rMATS_upstreamES\t$rMATS_upstreamEE\t";
      print OUTFILE "$rMATS_downstreamES\t$rMATS_downstreamEE\t$rMATS_exonStart_0base\t$rMATS_exonEnd\t";
   }
   print OUTFILE "$rMATS_ID\t$rMATS_IC_SAMPLE_1\t$rMATS_SC_SAMPLE_1\t$rMATS_IC_SAMPLE_2\t";
   print OUTFILE "$rMATS_SC_SAMPLE_2\t$rMATS_IncFormLen\t$rMATS_SkipFormLen\t$rMATS_PValue\t";
   print OUTFILE "$rMATS_FDR\t$rMATS_psi1\t$rMATS_psi2\t$rMATS_deltaPSI\n";
}

close(INFILE);
close(OUTFILE);
