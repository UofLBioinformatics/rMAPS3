#
## this program generates RNA Map using CLIP-seq and set of exons
#
### import necessary libraries
import re, os, sys, logging, time, datetime, csv, subprocess, numpy
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pyx import *
from scipy import stats
from rmaps_core.drawutils import export_canvas_outputs, export_text_rnamap_fallback, suppress_pyx_text_cleanup
from rmaps_core.stat_utils import compute_locus_pvalue, normalize_stat_method, pvalue_header_label
#import fisher,mne;  ## for FDR calculation
#
#
### checking out the number of arguments
if (len(sys.argv) < 10): 
  print('Not enough arguments!!');
  print ('It takes at least 9 arguments.');
  print ('Usage:\n\tProgramName exonPath peakCallerOutput intronLength exonLength windowSize stepSize FDR_cutoff proteinName outputFolder');
  print ('Example\n\tProgramName /exon/folder PIPE-CLIP.Clusters.bed 250 50 10 5 0.05 PTB RNA_map.PTB');
  sys.exit();

def listToString(x): ## log command
  rVal = '';
  for a in x:
    rVal += a+' ';
  return rVal;

## OUTPUT, create an output folder
outDir = sys.argv[9]
os.makedirs(outDir, exist_ok=True)
#
outPath = os.path.abspath(outDir) ## absolute output path

### setting up the logging format 
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(message)s',
                    filename=outPath+'/log.RNA.map.'+'.txt',
                    filemode='w')

##### Getting Start Time ######
logging.debug('Start the program with [%s]\n', listToString(sys.argv));
startTime = time.time();

#
### assigning parameters
r_exonPath = sys.argv[1];
exonPath = os.path.abspath(r_exonPath);
pFile = open(sys.argv[2]); ## peak caller output from CLIP seq
iLen = int(sys.argv[3]); ## length of intron to examine
eLen = int(sys.argv[4]); ## exon length to examine
wLen = int(sys.argv[5]); ## window length
sLen = int(sys.argv[6]) ## step size
nu = int(sys.argv[10]) if len(sys.argv) > 10 else 0
nd = int(sys.argv[11]) if len(sys.argv) > 11 else 0
nb = int(sys.argv[12]) if len(sys.argv) > 12 else 0
separate = int(sys.argv[13]) if len(sys.argv) > 13 else 0
event_type = sys.argv[14].upper() if len(sys.argv) > 14 else 'MXE'
stat_method = normalize_stat_method(os.environ.get('RMAPS_STAT_METHOD', 'fisher'))
#
## OUTPUT, create an output folder
#outDir = sys.argv[9]
#os.makedirs(outDir, exist_ok=True)
#
#outPath = os.path.abspath(outDir) ## absolute output path
#
ueFile = open(exonPath+'/up.coord.txt');
deFile = open(exonPath+'/dn.coord.txt');
beFile = open(exonPath+'/bg.coord.txt');

###### global variables ###################
#
u={};d={};b={}; ## dictionaries for upregulate, downregulated, and background exons
rmats={}; peaks={}; ## for rmats and peak caller output 
chunk=10000; ## to speed up
upCounts={}; dnCounts={}; bgCounts={}; ## dictionary for RNA map making
fdr_cutoff=float(sys.argv[7]);
proteinLabel=sys.argv[8];
boxHeight=iLen; ## height of the box for the plot
cdist_up={}; cdist_dn={}; cdist_bg={}; ## count distribution for wilcoxon rank-sum test 
test_up={}; test_dn={}; ## test stats (p-values, wilcoxon's rank sum)
#
#
###########################################

#
### functions here.. ##########
#
def readPeaks(pF): ### read peaks and height
  global peaks;
  totalPeaks=0; 
  for line in pF: ## process each line
    if line.strip()[0]=='#': ##comment, skip this line
      continue; ## process next line
    ele = line.strip().split('\t'); 
    chr=ele[0]; pStart=int(ele[1])+1; pEnd=int(ele[2]); height=int(float(ele[4]));  strand=ele[5]; 
    for pChunk in range(pStart//chunk, (pEnd//chunk)+1): ## for all possible chung the peak covers 
      if chr in peaks: ## already processed this chromosome
        if pChunk in peaks[chr]: ## we already have this chunk in this chr
          peaks[chr][pChunk].append([pStart,pEnd,height]); ## add it to the list
          totalPeaks += 1; 
        else: ## this chunk is not in this chr
          peaks[chr][pChunk] = [[pStart,pEnd,height]];
          totalPeaks += 1; 
      else: ## first time processing this chromosome
        peaks[chr]={};
        peaks[chr][pChunk]=[[pStart,pEnd,height]];
        totalPeaks += 1; 
  logging.debug("Done populating peaks dictionary with %d items" % totalPeaks);

def getExons(eFile, edic):
  line = eFile.readline(); ## header
  for line in eFile:
    ele=line.strip().split('\t');
    key = ':'.join([ele[0],ele[1],ele[2],ele[3],ele[4],ele[5]]);
    edic[key] = ['NA','NA','NA','NA','NA',ele[6],ele[7],ele[8],ele[9]];
  logging.debug("getExons function populated dictionary with %d items" % len(edic));


def read_OS_File(): ## read Order Statistics output. find up, down, background exons
  global u; global d; global b; #global rmats;
  global ueFile, deFile, beFile

  getExons(ueFile,u);
  getExons(deFile,d);
  getExons(beFile,b);

  logging.debug("There are %d, %d, %d exons in up, down, and background-exons" % (len(u), len(d), len(b)));


def countPeaks(ch,st,region,rv,ws,ss,countDist,ei): ## fill up rv and countDist with proper values, exon index was passed for wilcoxon test

  global peaks; global chunk;

  if ch not in peaks: ## no peak information about this chromosome. return to the program
    return;

  tC=[0]*len(rv); ## temp count

  if st=='+': ### positive strand
    rStart=region[0];
    for i in range(len(rv)): ## we already know how many times we need to do
      td={}; ## temp dictionary to avoid duplicate counting of peaks
      rEnd=rStart+ws-1;
      for rChunk in (rStart//chunk,rEnd//chunk): ## it examines all chunks
        if rChunk in peaks[ch]: ## there are peaks in this range
          for tP in peaks[ch][rChunk]: ## examine all tempPeaks
            if tP[0]<=rEnd and tP[1]>=rStart: ## overlap here
              td[str(tP[0])+':'+str(tP[1])]=tP[2]; ## height, it's okay to overwrite it
      for tk in td: ## for all non-duplicate peak counts, add it to the position
        tC[i]+=td[tk];  
      rStart += ss;   
    for i in range(len(rv)): ## now add it to the original value
      rv[i] += tC[i];
      countDist[i][ei] += tC[i];

  else: ## negative strand
    rEnd=region[1];
    for i in range(len(rv)): ## we already know how many times we need to do
      td={}; ## temp dictionary to avoid duplicate counting of peaks
      rStart=rEnd-ws+1;
      for rChunk in (rStart//chunk,rEnd//chunk): ## it examines all chunks
        if rChunk in peaks[ch]: ## there are peaks in this range
          for tP in peaks[ch][rChunk]: ## examine all tempPeaks
            if tP[0]<=rEnd and tP[1]>=rStart: ## overlap here
              td[str(tP[0])+':'+str(tP[1])]=tP[2]; ## height, it's okay to overwrite it
      for tk in td: ## for all non-duplicate peak counts, add it to the position
        tC[i]+=td[tk];
      rEnd -= ss;
    for i in range(len(rv)): ## now add it to the original value
      rv[i] += tC[i];
      countDist[i][ei] += tC[i];

def processExons(exons,peaks,CD,intronLen=250,exonLen=50,winSize=50,step=10): ## now process exons in each group

  #global cdist_up, cdist_dn, cdist_bg; ## get count distribution dictionaries

  eP=1+(exonLen-winSize)//step; ## exon points
  iP=1+(intronLen-winSize)//step; ##  intron points
  exonNum=len(exons);

  rVal={}; ## count for 8 regions at each step
  rVal[0] = [0] * eP;
  rVal[1] = [0] * iP;
  rVal[2] = [0] * iP;
  rVal[3] = [0] * eP;
  rVal[4] = [0] * eP;
  rVal[5] = [0] * iP;
  rVal[6] = [0] * iP;
  rVal[7] = [0] * eP;
  rVal[8] = [0] * eP;
  rVal[9] = [0] * iP;
  rVal[10] = [0] * iP;
  rVal[11] = [0] * eP;

  ## initializing the count distribution for wilcoxon rank sum
  CD[0]={}; 
  CD[1]={}; 
  CD[2]={}; 
  CD[3]={}; 
  CD[4]={}; 
  CD[5]={}; 
  CD[6]={}; 
  CD[7]={};
  CD[8]={};
  CD[9]={};
  CD[10]={};
  CD[11]={};
##  CD[0]=[[0]*exonNum]*eP;
##  CD[1]=[[0]*exonNum]*iP;
##  CD[2]=[[0]*exonNum]*iP;
##  CD[3]=[[0]*exonNum]*eP;
##  CD[4]=[[0]*exonNum]*eP;
##  CD[5]=[[0]*exonNum]*iP;
##  CD[6]=[[0]*exonNum]*iP;
##  CD[7]=[[0]*exonNum]*eP;

  for epep in range(eP): ## for 0,3,4,7
    for epindex in [0,3,4,7,8,11]: ## that has eP points
      CD[epindex][epep]=[0]*exonNum;
  for ipip in range(iP): ## for 1,2,5,6
    for ipindex in [1,2,5,6,9,10]: ## that has iP points
      CD[ipindex][ipip]=[0]*exonNum;

  exonIndex=0;
  for aa in exons: ## process each exon in this group
    ## 8 regions in 1-base coord
    eEle=aa.split(':');
    chr = eEle[0]; strand = eEle[1];
    firstES=int(eEle[2])+1;
    firstEE=int(eEle[3]);
    secondES = int(eEle[4]) + 1;
    secondEE = int(eEle[5]);
    upES=int(exons[aa][5])+1;
    upEE = int(exons[aa][6]) ;
    dnES = int(exons[aa][7]) + 1;
    dnEE = int(exons[aa][8]) ;

    r1 = [upEE - exonLen + 1, upEE-1];
    r2 = [upEE, upEE + intronLen ];
    r3 = [firstES - intronLen + 1, firstES - 1];
    r4 = [firstES, firstES + exonLen ];
    r5 = [firstEE - exonLen + 1, firstEE-1];
    r6 = [firstEE, firstEE + intronLen ];
    r7 = [secondES - intronLen + 1, secondES - 1];
    r8 = [secondES, secondES + exonLen];
    r9 = [secondEE - exonLen + 1, secondEE - 1];
    r10 = [secondEE, secondEE + intronLen];
    r11 = [dnES - intronLen + 1, dnES - 1];
    r12 = [dnES, dnES + exonLen];



    if strand=='-': ## for negative strand


      r12= [upEE - exonLen + 1, upEE - 1];
      r11= [upEE, upEE + intronLen];
      r10 = [firstES - intronLen + 1, firstES - 1];
      r9 = [firstES, firstES + exonLen];
      r8 = [firstEE - exonLen + 1, firstEE - 1];
      r7 = [firstEE, firstEE + intronLen];
      r6 = [secondES - intronLen + 1, secondES - 1];
      r5 = [secondES, secondES + exonLen];
      r4 = [secondEE - exonLen + 1, secondEE - 1];
      r3 = [secondEE, secondEE + intronLen];
      r2 = [dnES - intronLen + 1, dnES - 1];
      r1 = [dnES, dnES + exonLen];

    #1st region (upstream exon end)
    countPeaks(chr,strand,r1,rVal[0],winSize,step,CD[0],exonIndex);
    countPeaks(chr,strand,r2,rVal[1],winSize,step,CD[1],exonIndex);
    countPeaks(chr,strand,r3,rVal[2],winSize,step,CD[2],exonIndex);
    countPeaks(chr,strand,r4,rVal[3],winSize,step,CD[3],exonIndex);
    countPeaks(chr,strand,r5,rVal[4],winSize,step,CD[4],exonIndex);
    countPeaks(chr,strand,r6,rVal[5],winSize,step,CD[5],exonIndex);
    countPeaks(chr,strand,r7,rVal[6],winSize,step,CD[6],exonIndex);
    countPeaks(chr,strand,r8,rVal[7],winSize,step,CD[7],exonIndex);
    countPeaks(chr, strand, r9, rVal[8], winSize, step, CD[8], exonIndex);
    countPeaks(chr, strand, r10, rVal[9], winSize, step, CD[9], exonIndex);
    countPeaks(chr, strand, r11, rVal[10], winSize, step, CD[10], exonIndex);
    countPeaks(chr, strand, r12, rVal[11], winSize, step, CD[11], exonIndex);
    
    exonIndex += 1;

  logging.debug("Done processing exons in a group");
  return rVal; ## list with total count for each region

def computePValues(cdist_one, cdist_two, test_p): ## count p value for one vs. two
  myFactor = 1.0 * wLen
  for zz in sorted(cdist_one):
    test_p[zz] = {}
    for locus in range(len(cdist_one[zz])):
      test_p[zz][locus] = compute_locus_pvalue(
        cdist_one[zz][locus],
        cdist_two[zz][locus],
        stat_method,
        fisher_scale=myFactor,
      )


def computeWilcoxonP(cdist_one, cdist_two, test_p):
  computePValues(cdist_one, cdist_two, test_p)


def printCountDist(cds, desFile, exNum): ## print count distribution on the destination file, exNum is the nubmer of exons in the group
  rName={0:'R1', 1:'R2', 2:'R3', 3:'R4', 4:'R5', 5:'R6', 6:'R7', 7:'R8',8:'R9', 9:'R10', 10:'R11', 11:'R12'};
  desFile.write('Region\tposition\tsum\tvalues\n');
  if exNum==0: ## no exons in the group
    return;
  for zz in range(12): ## for 8 regions
    for locus in range(len(cds[zz])):
      desFile.write(rName[zz]+'\t'+str(locus)+'\t'+str(sum(cds[zz][locus]))+'\t[');
      desFile.write(','.join([str(ccc) for ccc in cds[zz][locus]])+']\n');


def printPval(pdic, dFile, eNum): ## print p values per position
  rName={0:'R1', 1:'R2', 2:'R3', 3:'R4', 4:'R5', 5:'R6', 6:'R7', 7:'R8',8:'R9', 9:'R10', 10:'R11', 11:'R12'};
  header = pvalue_header_label(stat_method)
  dFile.write('Region\tposition\t' + header + '\n');
  if eNum==0: ## no exons in the group
    return;
  for zz in range(12): ## for 8 regions
    for locus in range(len(pdic[zz])):
      dFile.write(rName[zz]+'\t'+str(locus)+'\t'+str(pdic[zz][locus])+'\n');


def printCounts(counts, dFile, eNum): ## print counts on the destination file, eNum is the nubmer of exons in the group
  rName={0:'R1', 1:'R2', 2:'R3', 3:'R4', 4:'R5', 5:'R6', 6:'R7', 7:'R8',8:'R9', 9:'R10', 10:'R11', 11:'R12'};
  dFile.write('Region\tPosition\tRawPeakCounts\tAverageRawPeakCount\n');
  if eNum==0: ## no exons in the group
    return;
  for zz in range(12): ## for 8 regions
    for ind in range(len(counts[zz])):
      dFile.write('\t'.join([rName[zz],str(ind),str(counts[zz][ind]),str(float(counts[zz][ind])/float(eNum))])+'\n');
    #dFile.write('\n');


def printCombinedCounts(upc,dnc,bgc,destFile,uNum,dNum,bNum): ## print combined counts
  rName={0:'R1', 1:'R2', 2:'R3', 3:'R4', 4:'R5', 5:'R6', 6:'R7', 7:'R8',8:'R9', 9:'R10', 10:'R11', 11:'R12'};
  destFile.write('Region\tPosition\tRawPeakCounts_up\tRawPeakCounts_down\tRawPeakCounts_background\tAverageRawPeakCount_up\tAverageRawPeakCount_down\tAverageRawPeakCount_background\n');
  if uNum*dNum*bNum==0: ## no exons in the group
    logging.debug("No exons in one of the upregulated, downregulated, or background exon groups. No combined result available.");
    return;
  for zz in range(12): ## for 8 regions
    for ind in range(len(upc[zz])): ## everything has the same items. It's okay to use this.
      outString = '\t'.join([rName[zz],str(ind),str(upc[zz][ind]),str(dnc[zz][ind]),str(bgc[zz][ind])])+'\t' ;
      outString += '\t'.join([str(float(upc[zz][ind])/float(uNum)),str(float(dnc[zz][ind])/float(dNum)),str(float(bgc[zz][ind])/float(bNum))]);
      destFile.write(outString+'\n');
    #dFile.write('\n');

def drawNode(c,eH,eW,iW,indent,gap,sGap,scale): ## draw node


  mls = style.linewidth.THIck;  ## my line style

  upstreamExon = path.path(path.moveto(0, 0), path.lineto(eW * scale, 0), path.lineto(scale * eW, scale * eH),
                           path.lineto(0, scale * eH), path.lineto(scale * indent, scale * eH / 2), path.lineto(0, 0));
  intron_1 = path.line(scale * eW, scale * eH / 2, scale * (eW + iW), scale * eH / 2);
  divider_1 = [
    path.line((eW + iW - sGap) * scale, (eH / 2 - gap) * scale, (eW + iW + sGap) * scale, (eH / 2 + gap) * scale),
    path.line((eW + iW - sGap + gap) * scale, (eH / 2 - gap) * scale, (eW + iW + sGap + gap) * scale,
              (eH / 2 + gap) * scale)];
  intron_2 = path.line((eW + iW + gap) * scale, scale * eH / 2, (eW + iW + gap + iW) * scale, scale * eH / 2);

  firstExon = path.rect((eW + iW + gap + iW) * scale, 0, (eW + gap + gap + eW) * scale, eH * scale);

  newX = eW + iW + gap + iW + eW + gap + gap + eW;  ## new start of x
  intron_3 = path.line(newX * scale, (eH / 2) * scale, (newX + iW) * scale, scale * eH / 2);
  divider_2 = [
    path.line((newX + iW - sGap) * scale, (eH / 2 - gap) * scale, (newX + iW + sGap) * scale, (eH / 2 + gap) * scale),
    path.line((newX + iW - sGap + gap) * scale, (eH / 2 - gap) * scale, (newX + iW + sGap + gap) * scale,
              (eH / 2 + gap) * scale)];
  intron_4 = path.line((newX + iW + gap) * scale, scale * eH / 2, (newX + iW + gap + iW) * scale, (eH / 2) * scale);
  secondExon = path.rect((newX + iW + gap + iW) * scale, 0, (eW + gap + gap + eW) * scale, eH * scale);
  newX = newX + iW + gap + iW + eW + gap + gap + eW;  ## new start of x
  intron_5 = path.line(newX * scale, (eH / 2) * scale, (newX + iW) * scale, scale * eH / 2);
  divider_3 = [
    path.line((newX + iW - sGap) * scale, (eH / 2 - gap) * scale, (newX + iW + sGap) * scale, (eH / 2 + gap) * scale),
    path.line((newX + iW - sGap + gap) * scale, (eH / 2 - gap) * scale, (newX + iW + sGap + gap) * scale,
              (eH / 2 + gap) * scale)];
  intron_6 = path.line((newX + iW + gap) * scale, scale * eH / 2, (newX + iW + gap + iW) * scale, (eH / 2) * scale);

  newX = newX + iW + gap + iW;  ## newX
  downstreamExon = path.path(path.moveto(newX * scale, 0), path.lineto((newX + eW) * scale, 0),
                             path.lineto((newX + eW - indent) * scale, scale * eH / 2),
                             path.lineto((newX + eW) * scale, eH * scale), path.lineto(newX * scale, eH * scale),
                             path.lineto(newX * scale, 0));

  c.stroke(upstreamExon, [mls, deco.filled([color.gray(0.8)])]);
  c.stroke(firstExon, [mls, deco.filled([color.rgb.green])]);
  c.stroke(secondExon, [mls, deco.filled([color.rgb.white])]);

  c.stroke(intron_1, [mls]);
  c.stroke(divider_1[0], [mls]);
  c.stroke(divider_1[1], [mls]);
  c.stroke(intron_2, [mls]);
  c.stroke(intron_3, [mls]);
  c.stroke(divider_2[0], [mls]);
  c.stroke(divider_2[1], [mls]);
  c.stroke(intron_4, [mls]);
  c.stroke(intron_5, [mls]);
  c.stroke(divider_3[0], [mls]);
  c.stroke(divider_3[1], [mls]);
  c.stroke(intron_6, [mls]);

  c.stroke(downstreamExon, [mls, deco.filled([color.gray(0.8)])]);

  logging.debug("Done drawing node structure");


def drawBox(c,mpv,eH,eW,iW,indent,gap,sGap,scale,el,il,maxNegP, separate=False): ## draw plot box
  # global motifLabel;
  global boxHeight;

  mls = style.linewidth.THIck;  ## my line style
  bgmls = style.linewidth.Thick;  ## background line style
  ldash = style.linestyle.dashed;
  largeLab = text.size.huge
  titleLab = text.size.Huge
  yLabAtt = [largeLab, text.halign.boxright, text.valign.middle];
  ypLabAtt = [largeLab, text.halign.boxleft, text.valign.middle];

  upColor = color.rgb.red;
  dnColor = color.rgb.blue;
  bgColor = color.rgb.black;

  boxY = eH + eH;
  boxHeight = iW * 1.2;
  bH = boxHeight;
  bW = eW + iW;

  myMax = mpv * 1.01;
  yLab = ['0.0', str("%.3f" % (myMax * 0.25)), str("%.3f" % (myMax * 0.5)), str("%.3f" % (myMax * 0.75)),
          str("%.3f" % float(myMax))];
  ypLab = ['0.0', str("%.3f" % (maxNegP * 0.25)), str("%.3f" % (maxNegP * 0.5)), str("%.3f" % (maxNegP * 0.75)),
           str("%.3f" % float(maxNegP))];

  xLab_1_3 = [str(-el), '0', str(il / 2), str(il)];  ## for 1st and 3rd boxes
  xLab_2_4 = [str(-il), str(-il / 2), '0', str(el)];  ## for the 2nd and 4th boxes

  xS = 0;
  xE = xS + bW;
  yl_offset = 5;
  xl_offset = 23;
  rect1 = path.rect(xS * scale, boxY * scale, bW * scale, bH * scale);
  r1_line = [path.line(xS * scale, (boxY + bH / 4.0) * scale, xE * scale, (boxY + bH / 4.0) * scale),
             path.line(xS * scale, (boxY + bH / 2.0) * scale, xE * scale, (boxY + bH / 2.0) * scale),
             path.line(xS * scale, (boxY + 3 * bH / 4.0) * scale, xE * scale, (boxY + 3 * bH / 4.0) * scale)];
  r1_ss = path.line((xS + eW) * scale, boxY * scale, (xS + eW) * scale, (boxY + boxHeight) * scale);
  ### print labels
  # y-axis
  c.text((xS - yl_offset) * scale, boxY * scale, yLab[0], yLabAtt);
  c.text((xS - yl_offset) * scale, (boxY + bH / 4.0) * scale, yLab[1], yLabAtt);
  c.text((xS - yl_offset) * scale, (boxY + bH / 2.0) * scale, yLab[2], yLabAtt);
  c.text((xS - yl_offset) * scale, (boxY + 3 * bH / 4.0) * scale, yLab[3], yLabAtt);
  c.text((xS - yl_offset) * scale, (boxY + bH) * scale, yLab[4], yLabAtt);
  ## y-label
  c.text((xS - yl_offset * 15) * scale, (boxY + bH / 2.0) * scale, "Average CLIP Peak Height",
         [largeLab, text.halign.boxcenter, trafo.rotate(90)]);
  # x-axis
  c.text(xS * scale, (boxY - xl_offset) * scale, xLab_1_3[0], [largeLab]);
  c.text((xS + eW) * scale, (boxY - xl_offset) * scale, xLab_1_3[1], [largeLab]);
  c.text((xS + eW + iW / 2.0) * scale, (boxY - xl_offset) * scale, xLab_1_3[2], [largeLab]);
  c.text((xS + eW + iW) * scale, (boxY - xl_offset) * scale, xLab_1_3[3], [largeLab, text.halign.boxright]);

  xS = xE + gap;
  xE = xS + bW;
  rect2 = path.rect(xS * scale, boxY * scale, bW * scale, bH * scale);
  r2_line = [path.line(xS * scale, (boxY + bH / 4.0) * scale, xE * scale, (boxY + bH / 4.0) * scale),
             path.line(xS * scale, (boxY + bH / 2.0) * scale, xE * scale, (boxY + bH / 2.0) * scale),
             path.line(xS * scale, (boxY + 3 * bH / 4.0) * scale, xE * scale, (boxY + 3 * bH / 4.0) * scale)];
  r2_ss = path.line((xS + iW) * scale, boxY * scale, (xS + iW) * scale, (boxY + boxHeight) * scale);

  # x-axis
  c.text(xS * scale, (boxY - xl_offset) * scale, xLab_2_4[0], [largeLab]);
  c.text((xS + iW / 2.0) * scale, (boxY - xl_offset) * scale, xLab_2_4[1], [largeLab]);
  c.text((xS + iW) * scale, (boxY - xl_offset) * scale, xLab_2_4[2], [largeLab]);
  c.text((xS + eW + iW) * scale, (boxY - xl_offset) * scale, xLab_2_4[3], [largeLab, text.halign.boxright]);

  ## Title and legends
  # c.text((xE+gap)*scale, (boxY+bH+45)*scale, motifLabel + " Motif MAP", [titleLab, text.halign.boxcenter]);
  c.text((xE + gap) * scale, (boxY + bH + 85) * scale, "RNA MAP FOR " + proteinLabel,
         [titleLab, text.halign.boxcenter]);

  c.stroke(
    path.line((xE + gap * 2) * scale, (boxY + bH + 25) * scale, (xE + gap * 3) * scale, (boxY + bH + 25) * scale),
    [upColor, mls]);
  # c.text((xE+gap*4)*scale, (boxY+bH+20)*scale, "Upregulated ("+str(totalExonCount['up'])+')', [upColor,largeLab]);
  c.text((xE + gap * 4) * scale, (boxY + bH + 20) * scale, "Upregulated(" + str(nu) + ")", [upColor, largeLab]);
  c.stroke(path.line((xE + gap + gap + iW * 2 / 3 + 3 * gap) * scale, (boxY + bH + 25) * scale,
                     (xE + gap + gap + iW * 2 / 3 + gap + 3 * gap) * scale, (boxY + bH + 25) * scale), [dnColor, mls]);
  # c.text((xE+gap+gap+iW*2/3+gap+gap)*scale, (boxY+bH+20)*scale, "Downregulated ("+str(totalExonCount['dn'])+')', [dnColor,largeLab]);
  c.text((xE + gap + gap + iW * 2 / 3 + gap + gap + 3 * gap) * scale, (boxY + bH + 20) * scale,
         "Downregulated(" + str(nd) + ")",
         [dnColor, largeLab]);
  c.stroke(path.line((xE + gap + iW + gap * 9 + 6 * gap) * scale, (boxY + bH + 25) * scale,
                     (xE + gap + iW + gap * 10 + 6 * gap) * scale,
                     (boxY + bH + 25) * scale), [bgColor, bgmls]);
  # c.text((xE+gap+iW+gap*11)*scale, (boxY+bH+20)*scale, "Background ("+str(totalExonCount['bg'])+')', [bgColor,largeLab]);
  c.text((xE + gap + iW + gap * 11 + 6 * gap) * scale, (boxY + bH + 20) * scale, "Background(" + str(nb) + ")",
         [bgColor, largeLab]);

  c.stroke(
    path.line((xE + gap * 2) * scale, (boxY + bH + 55) * scale, (xE + gap * 3) * scale, (boxY + bH + 55) * scale),
    [upColor, ldash]);
  c.text((xE + gap * 4) * scale, (boxY + bH + 50) * scale, "-log(pVal) up vs. bg", [upColor, largeLab]);
  c.stroke(
    path.line((xE + 8 * gap + iW * 2 / 3) * scale, (boxY + bH + 55) * scale, (xE + 8 * gap + iW * 2 / 3 + gap) * scale,
              (boxY + bH + 55) * scale), [dnColor, ldash]);
  c.text((xE + 8 * gap + iW * 2 / 3 + gap + gap) * scale, (boxY + bH + 50) * scale, "-log(pVal) dn vs. bg",
         [dnColor, largeLab]);

  xS = xE + gap + gap;
  xE = xS + bW;
  rect3 = path.rect(xS * scale, boxY * scale, bW * scale, bH * scale);
  r3_line = [path.line(xS * scale, (boxY + bH / 4.0) * scale, xE * scale, (boxY + bH / 4.0) * scale),
             path.line(xS * scale, (boxY + bH / 2.0) * scale, xE * scale, (boxY + bH / 2.0) * scale),
             path.line(xS * scale, (boxY + 3 * bH / 4.0) * scale, xE * scale, (boxY + 3 * bH / 4.0) * scale)];
  r3_ss = path.line((xS + eW) * scale, boxY * scale, (xS + eW) * scale, (boxY + boxHeight) * scale);
  # x-axis
  c.text(xS * scale, (boxY - xl_offset) * scale, xLab_1_3[0], [largeLab]);
  c.text((xS + eW) * scale, (boxY - xl_offset) * scale, xLab_1_3[1], [largeLab]);
  c.text((xS + eW + iW / 2.0) * scale, (boxY - xl_offset) * scale, xLab_1_3[2], [largeLab]);
  c.text((xS + eW + iW) * scale, (boxY - xl_offset) * scale, xLab_1_3[3], [largeLab, text.halign.boxright]);

  xS = xE + gap;
  xE = xS + bW;
  rect4 = path.rect(xS * scale, boxY * scale, bW * scale, bH * scale);
  r4_line = [path.line(xS * scale, (boxY + bH / 4.0) * scale, xE * scale, (boxY + bH / 4.0) * scale),
             path.line(xS * scale, (boxY + bH / 2.0) * scale, xE * scale, (boxY + bH / 2.0) * scale),
             path.line(xS * scale, (boxY + 3 * bH / 4.0) * scale, xE * scale, (boxY + 3 * bH / 4.0) * scale)];
  r4_ss = path.line((xS + iW) * scale, boxY * scale, (xS + iW) * scale, (boxY + boxHeight) * scale);
  # x-axis
  c.text(xS * scale, (boxY - xl_offset) * scale, xLab_2_4[0], [largeLab]);
  c.text((xS + iW / 2.0) * scale, (boxY - xl_offset) * scale, xLab_2_4[1], [largeLab]);
  c.text((xS + iW) * scale, (boxY - xl_offset) * scale, xLab_2_4[2], [largeLab]);
  c.text((xS + eW + iW) * scale, (boxY - xl_offset) * scale, xLab_2_4[3], [largeLab, text.halign.boxright]);

  xS = xE + gap + gap;
  xE = xS + bW;
  rect5 = path.rect(xS * scale, boxY * scale, bW * scale, bH * scale);
  r5_line = [path.line(xS * scale, (boxY + bH / 4.0) * scale, xE * scale, (boxY + bH / 4.0) * scale),
             path.line(xS * scale, (boxY + bH / 2.0) * scale, xE * scale, (boxY + bH / 2.0) * scale),
             path.line(xS * scale, (boxY + 3 * bH / 4.0) * scale, xE * scale, (boxY + 3 * bH / 4.0) * scale)];
  r5_ss = path.line((xS + eW) * scale, boxY * scale, (xS + eW) * scale, (boxY + boxHeight) * scale);
  # x-axis
  c.text(xS * scale, (boxY - xl_offset) * scale, xLab_1_3[0], [largeLab]);
  c.text((xS + eW) * scale, (boxY - xl_offset) * scale, xLab_1_3[1], [largeLab]);
  c.text((xS + eW + iW / 2.0) * scale, (boxY - xl_offset) * scale, xLab_1_3[2], [largeLab]);
  c.text((xS + eW + iW) * scale, (boxY - xl_offset) * scale, xLab_1_3[3], [largeLab, text.halign.boxright]);

  xS = xE + gap;
  xE = xS + bW;
  rect6 = path.rect(xS * scale, boxY * scale, bW * scale, bH * scale);
  r6_line = [path.line(xS * scale, (boxY + bH / 4.0) * scale, xE * scale, (boxY + bH / 4.0) * scale),
             path.line(xS * scale, (boxY + bH / 2.0) * scale, xE * scale, (boxY + bH / 2.0) * scale),
             path.line(xS * scale, (boxY + 3 * bH / 4.0) * scale, xE * scale, (boxY + 3 * bH / 4.0) * scale)];
  r6_ss = path.line((xS + iW) * scale, boxY * scale, (xS + iW) * scale, (boxY + boxHeight) * scale);
  # x-axis
  c.text(xS * scale, (boxY - xl_offset) * scale, xLab_2_4[0], [largeLab]);
  c.text((xS + iW / 2.0) * scale, (boxY - xl_offset) * scale, xLab_2_4[1], [largeLab]);
  c.text((xS + iW) * scale, (boxY - xl_offset) * scale, xLab_2_4[2], [largeLab]);
  c.text((xS + eW + iW) * scale, (boxY - xl_offset) * scale, xLab_2_4[3], [largeLab, text.halign.boxright]);

  c.stroke(rect1, [mls]);
  c.stroke(r1_line[0], [ldash]);
  c.stroke(r1_line[1], [ldash]);
  c.stroke(r1_line[2], [ldash]);
  c.stroke(rect2, [mls]);
  c.stroke(r2_line[0], [ldash]);
  c.stroke(r2_line[1], [ldash]);
  c.stroke(r2_line[2], [ldash]);
  c.stroke(rect3, [mls]);
  c.stroke(r3_line[0], [ldash]);
  c.stroke(r3_line[1], [ldash]);
  c.stroke(r3_line[2], [ldash]);
  c.stroke(rect4, [mls]);
  c.stroke(r4_line[0], [ldash]);
  c.stroke(r4_line[1], [ldash]);
  c.stroke(r4_line[2], [ldash]);
  c.stroke(rect5, [mls]);
  c.stroke(r5_line[0], [ldash]);
  c.stroke(r5_line[1], [ldash]);
  c.stroke(r5_line[2], [ldash]);
  c.stroke(rect6, [mls]);
  c.stroke(r6_line[0], [ldash]);
  c.stroke(r6_line[1], [ldash]);
  c.stroke(r6_line[2], [ldash]);

  ## lines for splice sites
  c.stroke(r1_ss);
  c.stroke(r2_ss);
  c.stroke(r3_ss);
  c.stroke(r4_ss);
  c.stroke(r5_ss);
  c.stroke(r6_ss);
  ### print pvalue labels
  # yp-axis
  xE = xE + eW;
  c.text((xE - 5 * yl_offset - gap) * scale, boxY * scale, ypLab[0], ypLabAtt);
  c.text((xE - 5 * yl_offset - gap) * scale, (boxY + bH / 4.0) * scale, ypLab[1], ypLabAtt);
  c.text((xE - 5 * yl_offset - gap) * scale, (boxY + bH / 2.0) * scale, ypLab[2], ypLabAtt);
  c.text((xE - 5 * yl_offset - gap) * scale, (boxY + 3 * bH / 4.0) * scale, ypLab[3], ypLabAtt);
  c.text((xE - 5 * yl_offset - gap) * scale, (boxY + bH) * scale, ypLab[4], ypLabAtt);
  ## y-label
  c.text((xE + yl_offset * 6) * scale, (boxY + bH / 2.0) * scale, "Negative log10(pValue)",
         [largeLab, text.halign.boxcenter, trafo.rotate(270)]);

  logging.debug("Done drawing plot box");


def oldplotRegions(tRegion_1,tRegion_2,xS,boxY,bH,fW,sW,mpv,scale,howmany=3): ## plot given regions, fW, sW tells the width of the first and the second region
  pup=[]; pdn=[]; pbg=[]; ## path for up, dn, and bg
  pvalup=[]; pvaldn=[]; ## path for pval up, pval dn
  y1=[];y2=[];

  for jj in range(len(tRegion_1)-1): ## for each point in the first region (upstream exon)
      x1=xS+jj*float(fW)/len(tRegion_1)+(float(fW)/len(tRegion_1))/2.0; x2=xS+(jj+1)*float(fW)/len(tRegion_1)+(float(fW)/len(tRegion_1))/2.0;
      y1=[boxY+tRegion_1[jj][0]*bH/mpv, boxY+tRegion_1[jj][1]*bH/mpv, boxY+tRegion_1[jj][2]*bH/mpv];
      y2=[boxY+tRegion_1[jj+1][0]*bH/mpv, boxY+tRegion_1[jj+1][1]*bH/mpv, boxY+tRegion_1[jj+1][2]*bH/mpv];

      pup.append([x1*scale, y1[0]*scale]);
      pdn.append([x1*scale, y1[1]*scale]);
      pbg.append([x1*scale, y1[2]*scale]);

      if jj==len(tRegion_1)-2: ## second to the last element. add coords for the last element
        pup.append([x2*scale, y2[0]*scale]);
        pdn.append([x2*scale, y2[1]*scale]);
        pbg.append([x2*scale, y2[2]*scale]);

  xS=xS+fW;
  for jj in range(len(tRegion_2)-1): ## for each point in the second region (intron)
      x1=xS+jj*float(sW)/len(tRegion_2)+(float(sW)/len(tRegion_2))/2.0; x2=xS+(jj+1)*float(sW)/len(tRegion_2)+(float(sW)/len(tRegion_2))/2.0;
      y1=[boxY+tRegion_2[jj][0]*bH/mpv, boxY+tRegion_2[jj][1]*bH/mpv, boxY+tRegion_2[jj][2]*bH/mpv];
      y2=[boxY+tRegion_2[jj+1][0]*bH/mpv, boxY+tRegion_2[jj+1][1]*bH/mpv, boxY+tRegion_2[jj+1][2]*bH/mpv];

      pup.append([x1*scale, y1[0]*scale]);
      pdn.append([x1*scale, y1[1]*scale]);
      pbg.append([x1*scale, y1[2]*scale]);

      if jj==len(tRegion_2)-2: ## second to the last element. add coords for the last element
        pup.append([x2*scale, y2[0]*scale]);
        pdn.append([x2*scale, y2[1]*scale]);
        pbg.append([x2*scale, y2[2]*scale]);

  logging.debug("Done plotRegions function");
  return pup, pdn, pbg;
#


def plotRegions(tRegion_1,tRegion_2,xS,boxY,bH,fW,sW,mpv,scale,howmany=3): ## plot given regions, fW, sW tells the width of the first and the second region
  pup=[]; pdn=[]; pbg=[]; ## path for up, dn, and bg
  pvalup=[]; pvaldn=[]; ## path for pval up, pval dn
  y1=[];y2=[];

  for jj in range(len(tRegion_1)-1): ## for each point in the first region (upstream exon)
      x1=xS+jj*float(fW)/len(tRegion_1)+(float(fW)/len(tRegion_1))/2.0; x2=xS+(jj+1)*float(fW)/len(tRegion_1)+(float(fW)/len(tRegion_1))/2.0;
      if howmany==3: ## actual points
        y1=[boxY+tRegion_1[jj][0]*bH/mpv, boxY+tRegion_1[jj][1]*bH/mpv, boxY+tRegion_1[jj][2]*bH/mpv];
        y2=[boxY+tRegion_1[jj+1][0]*bH/mpv, boxY+tRegion_1[jj+1][1]*bH/mpv, boxY+tRegion_1[jj+1][2]*bH/mpv];

        pup.append([x1*scale, y1[0]*scale]);
        pdn.append([x1*scale, y1[1]*scale]);
        pbg.append([x1*scale, y1[2]*scale]);

        if jj==len(tRegion_1)-2: ## second to the last element. add coords for the last element
          pup.append([x2*scale, y2[0]*scale]);
          pdn.append([x2*scale, y2[1]*scale]);
          pbg.append([x2*scale, y2[2]*scale]);

      elif howmany==2: ## pvalue points
        y1=[boxY+tRegion_1[jj][0]*bH/mpv, boxY+tRegion_1[jj][1]*bH/mpv];
        y2=[boxY+tRegion_1[jj+1][0]*bH/mpv, boxY+tRegion_1[jj+1][1]*bH/mpv];

        pvalup.append([x1*scale, y1[0]*scale]);
        pvaldn.append([x1*scale, y1[1]*scale]);

        if jj==len(tRegion_1)-2: ## second to the last element. add coords for the last element
          pvalup.append([x2*scale, y2[0]*scale]);
          pvaldn.append([x2*scale, y2[1]*scale]);

  xS=xS+fW;
  for jj in range(len(tRegion_2)-1): ## for each point in the second region (intron)
      x1=xS+jj*float(sW)/len(tRegion_2)+(float(sW)/len(tRegion_2))/2.0; x2=xS+(jj+1)*float(sW)/len(tRegion_2)+(float(sW)/len(tRegion_2))/2.0;
      if howmany==3: ## actual points
        y1=[boxY+tRegion_2[jj][0]*bH/mpv, boxY+tRegion_2[jj][1]*bH/mpv, boxY+tRegion_2[jj][2]*bH/mpv];
        y2=[boxY+tRegion_2[jj+1][0]*bH/mpv, boxY+tRegion_2[jj+1][1]*bH/mpv, boxY+tRegion_2[jj+1][2]*bH/mpv];

        pup.append([x1*scale, y1[0]*scale]);
        pdn.append([x1*scale, y1[1]*scale]);
        pbg.append([x1*scale, y1[2]*scale]);

        if jj==len(tRegion_2)-2: ## second to the last element. add coords for the last element
          pup.append([x2*scale, y2[0]*scale]);
          pdn.append([x2*scale, y2[1]*scale]);
          pbg.append([x2*scale, y2[2]*scale]);
      elif howmany==2: ## pvalue points
        y1=[boxY+tRegion_2[jj][0]*bH/mpv, boxY+tRegion_2[jj][1]*bH/mpv];
        y2=[boxY+tRegion_2[jj+1][0]*bH/mpv, boxY+tRegion_2[jj+1][1]*bH/mpv];

        pvalup.append([x1*scale, y1[0]*scale]);
        pvaldn.append([x1*scale, y1[1]*scale]);

        if jj==len(tRegion_2)-2: ## second to the last element. add coords for the last element
          pvalup.append([x2*scale, y2[0]*scale]);
          pvaldn.append([x2*scale, y2[1]*scale]);

  logging.debug("Done plotRegions function");
  if howmany==3:
    return pup, pdn, pbg;
  elif howmany==2:
    return pvalup, pvaldn;

#

def fillUpPath(tPoints):
  rPath=path.path(path.moveto(tPoints[0][0], tPoints[0][1]));
  for pp in range(1,len(tPoints)): ## for each point from the 2nd element
    rPath.append(path.lineto(tPoints[pp][0], tPoints[pp][1]));
  return rPath;
  logging.debug("Done fillUpPath function");

def drawAcutalPlot(c,dp,mpv,eH,eW,iW,gap,sGap,scale,wl,sl,nPP,mNP, separate=False): ## draw plots now..nPP is negativePvaluePoins, mNP is maximumNegPval
  global boxHeight;

  mls = style.linewidth.THIck;  ## my line style
  bgmls = style.linewidth.Thick;  ## my line style
  ldash = style.linestyle.dashed;
  largeLab = text.size.huge

  boxY = eH + eH;
  bH = boxHeight;
  bW = eW + iW;

  upColor = color.rgb.red;
  dnColor = color.rgb.blue;
  bgColor = color.rgb.black;

  ## 1st and 2nd region
  xS = 0;
  points_up, points_dn, points_bg = plotRegions(dp[0], dp[1], xS, boxY, bH, eW, iW, mpv, scale, 3);  ## return 3 values
  path_up = fillUpPath(points_up);
  path_dn = fillUpPath(points_dn);
  path_bg = fillUpPath(points_bg);
  c.stroke(path_up, [upColor, mls]);
  c.stroke(path_dn, [dnColor, mls]);
  c.stroke(path_bg, [bgColor, bgmls]);
  points_pup, points_pdn = plotRegions(nPP[0], nPP[1], xS, boxY, bH, eW, iW, mNP, scale, 2);  ## return 2 values
  path_pup = fillUpPath(points_pup);
  path_pdn = fillUpPath(points_pdn);
  c.stroke(path_pup, [upColor, ldash]);
  c.stroke(path_pdn, [dnColor, ldash]);

  ## 3rd and 4th region
  xS = xS + eW + iW + gap;
  points_up, points_dn, points_bg = plotRegions(dp[2], dp[3], xS, boxY, bH, iW, eW, mpv, scale);
  path_up = fillUpPath(points_up);
  path_dn = fillUpPath(points_dn);
  path_bg = fillUpPath(points_bg);
  c.stroke(path_up, [upColor, mls]);
  c.stroke(path_dn, [dnColor, mls]);
  c.stroke(path_bg, [bgColor, bgmls]);
  points_pup, points_pdn = plotRegions(nPP[2], nPP[3], xS, boxY, bH, iW, eW, mNP, scale, 2);  ## return 2 values
  path_pup = fillUpPath(points_pup);
  path_pdn = fillUpPath(points_pdn);
  c.stroke(path_pup, [upColor, ldash]);
  c.stroke(path_pdn, [dnColor, ldash]);

  ## 5th and 6th region
  xS = xS + iW + eW + gap + gap;
  points_up, points_dn, points_bg = plotRegions(dp[4], dp[5], xS, boxY, bH, eW, iW, mpv, scale);
  path_up = fillUpPath(points_up);
  path_dn = fillUpPath(points_dn);
  path_bg = fillUpPath(points_bg);
  c.stroke(path_up, [upColor, mls]);
  c.stroke(path_dn, [dnColor, mls]);
  c.stroke(path_bg, [bgColor, bgmls]);
  points_pup, points_pdn = plotRegions(nPP[4], nPP[5], xS, boxY, bH, eW, iW, mNP, scale, 2);  ## return 2 values
  path_pup = fillUpPath(points_pup);
  path_pdn = fillUpPath(points_pdn);
  c.stroke(path_pup, [upColor, ldash]);
  c.stroke(path_pdn, [dnColor, ldash]);

  ## 7th and8 th region
  xS = xS + eW + iW + gap;
  points_up, points_dn, points_bg = plotRegions(dp[6], dp[7], xS, boxY, bH, iW, eW, mpv, scale);
  path_up = fillUpPath(points_up);
  path_dn = fillUpPath(points_dn);
  path_bg = fillUpPath(points_bg);
  c.stroke(path_up, [upColor, mls]);
  c.stroke(path_dn, [dnColor, mls]);
  c.stroke(path_bg, [bgColor, bgmls]);
  points_pup, points_pdn = plotRegions(nPP[6], nPP[7], xS, boxY, bH, iW, eW, mNP, scale, 2);  ## return 2 values
  path_pup = fillUpPath(points_pup);
  path_pdn = fillUpPath(points_pdn);
  c.stroke(path_pup, [upColor, ldash]);
  c.stroke(path_pdn, [dnColor, ldash]);

  ## 9th and 10th region
  xS = xS + iW + eW + gap + gap;
  points_up, points_dn, points_bg = plotRegions(dp[8], dp[9], xS, boxY, bH, eW, iW, mpv, scale);
  path_up = fillUpPath(points_up);
  path_dn = fillUpPath(points_dn);
  path_bg = fillUpPath(points_bg);
  c.stroke(path_up, [upColor, mls]);
  c.stroke(path_dn, [dnColor, mls]);
  c.stroke(path_bg, [bgColor, bgmls]);
  points_pup, points_pdn = plotRegions(nPP[8], nPP[9], xS, boxY, bH, eW, iW, mNP, scale, 2);  ## return 2 values
  path_pup = fillUpPath(points_pup);
  path_pdn = fillUpPath(points_pdn);
  c.stroke(path_pup, [upColor, ldash]);
  c.stroke(path_pdn, [dnColor, ldash]);

  ## 11th and12 th region
  xS = xS + eW + iW + gap;
  points_up, points_dn, points_bg = plotRegions(dp[10], dp[11], xS, boxY, bH, iW, eW, mpv, scale);
  path_up = fillUpPath(points_up);
  path_dn = fillUpPath(points_dn);
  path_bg = fillUpPath(points_bg);
  c.stroke(path_up, [upColor, mls]);
  c.stroke(path_dn, [dnColor, mls]);
  c.stroke(path_bg, [bgColor, bgmls]);
  points_pup, points_pdn = plotRegions(nPP[10], nPP[11], xS, boxY, bH, iW, eW, mNP, scale, 2);  ## return 2 values
  path_pup = fillUpPath(points_pup);
  path_pdn = fillUpPath(points_pdn);
  c.stroke(path_pup, [upColor, ldash]);
  c.stroke(path_pdn, [dnColor, ldash]);

  logging.debug("Done drawing acutal plots");


def drawMap(upc,dnc,bgc,iL,eL,wL,sL,uNum,dNum,bNum, pdic_up, pdic_dn, separate=False): ## draw RNA map here. pdics are pvalue.up.vs.bg, pvalue.dn.vs.bg

  global outPath;

  ## setting up a canvas
  canv = canvas.canvas(); ## set canvas
  drawPoints=[]; 
  negPvalPoints=[];
  maxPointValue=0.0; maxNegPval=0.0;
  for zz in range(12): ## for 8 regions
    drawPoints.append([]);
    negPvalPoints.append([]);
    for ind in range(len(upc[zz])): ## everything has the same items. It's okay to use this.
       drawPoints[zz].append([float(upc[zz][ind])/float(uNum),float(dnc[zz][ind])/float(dNum),float(bgc[zz][ind])/float(bNum)]);
       negPvalPoints[zz].append([-numpy.log10(pdic_up[zz][ind]),-numpy.log10(pdic_dn[zz][ind])]);
       for yy in drawPoints[zz][ind]: ## examine three values
         if yy > maxPointValue: ## new max here
           maxPointValue=yy;
       for yy in negPvalPoints[zz][ind]: ## examine two values
         if yy>maxNegPval: ## new max here
           maxNegPval=yy;
  logging.debug("Max average height is: %f" % maxPointValue);
  logging.debug("Max negPval is: %f" % maxNegPval);

  exonHeight=30; exonWidth=50; intronWidth=200; indentation=5; dividerGap=10; slopeGap=20; Scale=0.05;
  drawNode(canv,exonHeight,exonWidth,intronWidth,indentation,dividerGap,slopeGap,Scale);
  drawBox(canv,maxPointValue,exonHeight,exonWidth,intronWidth,indentation,dividerGap,slopeGap,Scale,iL,eL,maxNegPval,separate);
  drawAcutalPlot(canv,drawPoints,maxPointValue,exonHeight,exonWidth,intronWidth,dividerGap,slopeGap,Scale,wL,sL,negPvalPoints,maxNegPval,separate);
  pdf_path = outPath + '/' + event_type + '.RNA.map.pdf'
  png_path = outPath + '/' + event_type + '.RNA.map.png'
  pdf_ok, png_ok = export_canvas_outputs(
    canv,
    pdf_path,
    png_path,
    png_resolution=100,
    logger=logging,
  )
  if pdf_ok:
    logging.debug("PDF export successful: %s", pdf_path)
  if png_ok:
    logging.debug("PNG export successful: %s", png_path)
  if not pdf_ok and not png_ok:
    logging.warning("Neither PDF nor PNG export succeeded for RNA map")
  logging.debug("Done drawing RNA map. Please find RNA.map output in the output folder.");


##################### main process #######################
#
### fillup peaks dictionary
try:
  readPeaks(pFile);  ## read peaks
except:
  logging.debug("There is an exception in readPeaks function");
  logging.debug("Exception: %s" % sys.exc_info()[0]);
  logging.debug("Detail: %s" % sys.exc_info()[1]);
  print("There is an exception. Please check your log file")
  sys.exit(-9);
try:
  #read_rMATS_File(rFile); ## read rMATS output to find upregulated, downregulated, and background exons
  read_OS_File(); ## read OrderStatisitics output to find upregulated, downregulated, and background exons
  
  if len(u)==0: ## no upregulated exon
    logging.debug("No upregulated exons. Cannot create RNA Map. Exiting..");  
    print("No upregulated exons. Cannot create RNA Map. Exiting..")
    sys.exit(-10);
  if len(d)==0: ## no downregulated exon
    logging.debug("No downregulated exons. Cannot create RNA Map. Exiting..");  
    print("No downregulated exons. Cannot create RNA Map. Exiting..")
    sys.exit(-10);
  if len(b)==0: ## no background exon
    logging.debug("No background exons. Cannot create RNA Map. Exiting..");  
    print("No background exons. Cannot create RNA Map. Exiting..")
    sys.exit(-10);

except:
  #logging.debug("There is an exception in read_rMATS_File function");
  logging.debug("There is an exception in read_OS_File function");
  logging.debug("Exception: %s" % sys.exc_info()[0]);
  logging.debug("Detail: %s" % sys.exc_info()[1]);
  print("There is an exception. Please check your log file")
  sys.exit(-9);

## now process exons in each group
try:
  upCounts=processExons(u,peaks,cdist_up,iLen,eLen,wLen,sLen);
except:
  logging.debug("There is an exception in processExons function with u");
  logging.debug("Exception: %s" % sys.exc_info()[0]);
  logging.debug("Detail: %s" % sys.exc_info()[1]);
  print("There is an exception. Please check your log file")
  sys.exit(-9);
try:
  dnCounts=processExons(d,peaks,cdist_dn,iLen,eLen,wLen,sLen);
except:
  logging.debug("There is an exception in processExons function with d");
  logging.debug("Exception: %s" % sys.exc_info()[0]);
  logging.debug("Detail: %s" % sys.exc_info()[1]);
  print("There is an exception. Please check your log file")
  sys.exit(-9);
try:
  bgCounts=processExons(b,peaks,cdist_bg,iLen,eLen,wLen,sLen); 
except:
  logging.debug("There is an exception in processExons function with b");
  logging.debug("Exception: %s" % sys.exc_info()[0]);
  logging.debug("Detail: %s" % sys.exc_info()[1]);
  print("There is an exception. Please check your log file")
  sys.exit(-9);
#

upFile = open(outPath+'/'+'upregulated.RNAmap.txt', 'w');
dnFile = open(outPath+'/'+'downregulated.RNAmap.txt', 'w');
bgFile = open(outPath+'/'+'background.RNAmap.txt', 'w');
cbFile = open(outPath+'/'+'combined.RNAmap.txt', 'w');

try:
  printCounts(upCounts, upFile, len(u)); ## print counts
  printCounts(dnCounts, dnFile, len(d)); ## print counts
  printCounts(bgCounts, bgFile, len(b)); ## print counts
  printCombinedCounts(upCounts,dnCounts,bgCounts,cbFile,len(u),len(d),len(b)); ## print combined counts
except:
  logging.debug("There is an exception in printCounts function");
  logging.debug("Exception: %s" % sys.exc_info()[0]);
  logging.debug("Detail: %s" % sys.exc_info()[1]);
  print("There is an exception. Please check your log file")
  sys.exit(-9);


updistFile = open(outPath+'/'+'count.dist.upregulated.RNAmap.txt', 'w');
dndistFile = open(outPath+'/'+'count.dist.downregulated.RNAmap.txt', 'w');
bgdistFile = open(outPath+'/'+'count.dist.background.RNAmap.txt', 'w');

try:
  printCountDist(cdist_up, updistFile, len(u)); ## print counts
  printCountDist(cdist_dn, dndistFile, len(d)); ## print counts
  printCountDist(cdist_bg, bgdistFile, len(b)); ## print counts
except:
  logging.debug("There is an exception in printCountDist function");
  logging.debug("Exception: %s" % sys.exc_info()[0]);
  logging.debug("Detail: %s" % sys.exc_info()[1]);
  print("There is an exception in printCountDist. Please check your log file")
  sys.exit(-9);


###### calculate Wilcoxon's rank sum p-value

try:
  computePValues(cdist_up, cdist_bg, test_up); ## count p value for up vs. bg
  computePValues(cdist_dn, cdist_bg, test_dn); ## count p value for dn vs. bg
except:
  logging.debug("There is an exception in computePValues function");
  logging.debug("Exception: %s" % sys.exc_info()[0]);
  logging.debug("Detail: %s" % sys.exc_info()[1]);
  print("There is an exception in computePValues. Please check your log file")
  sys.exit(-95);



upbgPFile = open(outPath+'/'+'pVal.up.vs.bg.RNAmap.txt', 'w');
dnbgPFile = open(outPath+'/'+'pVal.dn.vs.bg.RNAmap.txt', 'w');

try:
  printPval(test_up, upbgPFile, len(u)); ## print p value, up vs bg
  printPval(test_dn, dnbgPFile, len(d)); ## print p value, dn vs bg
except:
  logging.debug("There is an exception in printPval function");
  logging.debug("Exception: %s" % sys.exc_info()[0]);
  logging.debug("Detail: %s" % sys.exc_info()[1]);
  print("There is an exception in printPval. Please check your log file")
  sys.exit(-9);


try:
  drawMap(upCounts,dnCounts,bgCounts,eLen,iLen,wLen,sLen,len(u),len(d),len(b),test_up, test_dn,separate);  ## now draw RNA map using counts
  logging.info("RNA map visualization successfully generated")
except Exception:
  logging.warning("Could not generate RNA map visualization (PNG/PDF files)")
  logging.warning("This may be due to missing or blocked TeX/Ghostscript installation")
  logging.warning("Exception: %s" % sys.exc_info()[0])
  logging.warning("Detail: %s" % sys.exc_info()[1])
  suppress_pyx_text_cleanup(logging)
  fallback_pdf = os.path.join(outDir, event_type + '.RNA.map.pdf')
  fallback_png = os.path.join(outDir, event_type + '.RNA.map.png')
  pdf_ok, png_ok = export_text_rnamap_fallback(
      outDir, event_type, proteinLabel, fallback_pdf, fallback_png, logging)
  if pdf_ok or png_ok:
    logging.info("Generated Pillow RNA map fallback: pdf=%s png=%s", pdf_ok, png_ok)
    print("Generated fallback visualization files from RNA map text output")
  logging.info("RNA map data files were successfully generated even though visualization failed")
  print("Warning: Could not generate visualization files (TeX/Ghostscript may not be available)")
  print("RNA map data files were successfully generated in the output directory")


updistFile.close(); dndistFile.close(); bgdistFile.close();
upbgPFile.close(); dnbgPFile.close();
upFile.close(); dnFile.close(); bgFile.close(); cbFile.close();
ueFile.close(); deFile.close(); beFile.close();
pFile.close();
#rFile.close();

#############
## calculate total running time
#############
logging.debug("Program ended");
currentTime = time.time();
runningTime = currentTime-startTime; ## in seconds
logging.debug("Program ran %.2d:%.2d:%.2d" % (runningTime/3600, (runningTime%3600)/60, runningTime%60));

sys.exit(0);
