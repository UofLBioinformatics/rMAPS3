use strict;



##############################################################
sub parseIsoforms { 
  my($iso) = @_;
  my @wds = split(/\@/, $iso);
  my $numWds = @wds;

#  if($numWds != 11) { 
#     die("Error in parseIsoforms: 11 fields expected; $numWds encountered\n");
#  }
  my $longexon = $wds[0];
  my $flankingexon = $wds[1];

  my @wds2 = split(/\:/, $longexon);
  my $strand = $wds2[-1];
  my $longExonStart_0base = $wds2[-3];
  my $shortES = $wds2[-3];
  my $var   = $wds2[-2];
 
  @wds2 = split(/\|/, $var);
  my $longExonEnd = $wds2[-1];
  my $shortEE   = $wds2[-2];
  
  @wds2 = split(/\:/, $flankingexon);
#  my $strand = $wds2[3];
  my $flankingES = $wds2[1];
  my $flankingEE   = $wds2[2];

  


  my @retArr = ($longExonStart_0base, $longExonEnd, $shortES, $shortEE, 
                $flankingES, $flankingEE);
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

print OUTFILE "ID\tGeneID\tgeneSymbol\tchr\tstrand\tlongExonStart_0base\tlongExonEnd\t";
print OUTFILE "shortES\tshortEE\tflankingES\tflankingEE\tID\tIC_SAMPLE_1\t";
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
   my $rMATS_longExonStart_0base = 0;
   my $rMATS_longExonEnd         = 0;
   my $rMATS_shortES      = 0;
   my $rMATS_shortEE      = 0;
   my $rMATS_flankingES    = 0;
   my $rMATS_flankingEE    = 0; 
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

   ($rMATS_longExonStart_0base, $rMATS_longExonEnd, $rMATS_shortES,$rMATS_shortEE,
    $rMATS_flankingES, $rMATS_flankingEE) = parseIsoforms($isoforms, $strand);

	
	
				
   print OUTFILE "$rMATS_ID\t$rMATS_GeneID\t$rMATS_geneSymbol\t$rMATS_chr\t$rMATS_strand\t";

   if ($strand eq "+"){
      $rMATS_longExonStart_0base=$rMATS_longExonStart_0base-1;
	  $rMATS_shortES=$rMATS_shortES-1;
	  $rMATS_flankingES=$rMATS_flankingES-1;
      print OUTFILE "$rMATS_longExonStart_0base\t$rMATS_longExonEnd\t";
      print OUTFILE "$rMATS_shortES\t$rMATS_shortEE\t$rMATS_flankingES\t$rMATS_flankingEE\t";
   }
   if ($strand eq "-"){
      $rMATS_shortEE=$rMATS_shortEE-1;
	  $rMATS_longExonEnd=$rMATS_longExonEnd-1;
	  $rMATS_flankingES=$rMATS_flankingES-1;
      print OUTFILE "$rMATS_shortEE\t$rMATS_shortES\t";
      print OUTFILE "$rMATS_longExonEnd\t$rMATS_shortES\t$rMATS_flankingES\t$rMATS_flankingEE\t";
   }
   print OUTFILE "$rMATS_ID\t$rMATS_IC_SAMPLE_1\t$rMATS_SC_SAMPLE_1\t$rMATS_IC_SAMPLE_2\t";
   print OUTFILE "$rMATS_SC_SAMPLE_2\t$rMATS_IncFormLen\t$rMATS_SkipFormLen\t$rMATS_PValue\t";
   print OUTFILE "$rMATS_FDR\t$rMATS_psi1\t$rMATS_psi2\t$rMATS_deltaPSI\n";
}

close(INFILE);
close(OUTFILE);
