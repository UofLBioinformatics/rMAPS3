#
## this program generates motif Map for three lists of exons provided
#
### import necessary libraries
from pyx import *
## for making plot
import re, os, sys, logging, time, datetime, argparse, csv, numpy, subprocess, shutil
import itertools as it
import glob
import operator
#import mne;  ## for FDR calculation
from scipy import stats
import threading
from threading import Thread
import timeit
from multiprocessing import Process, Manager
from math import exp

import drawutils
from genome_access import load_genome, fetch_seq as fetch_seq_from_fasta
#pygr_station_PATH='/u/nobackup/yxing/jwpark/data/pygr/hg19/';
#pygr_station_PATH='/u/nobackup/yxing/jwpark/data/pygr/';
#
############ need to implement Wilcoxon rank-sum test, scipy.stats.ranksums(x,y) ###########


def run_command(cmd):
    completed = subprocess.run(cmd, capture_output=True, text=True)
    status = completed.returncode
    output = (completed.stdout or "") + (completed.stderr or "")
    return status, output


def copy_file(src, dst):
    try:
        shutil.copy2(src, dst)
        return 0, ""
    except Exception as exc:
        return 1, str(exc)


#
if __name__ == "__main__":
    # Legacy script entry point; prefer `python cli.py motif-map se ...` instead.
    pass

### checking out the number of arguments
parser = argparse.ArgumentParser(
    description=
    'Making motif map from a list of motifs and three lists of exon coordinates'
)
parser.add_argument('-m',
                    '--motif',
                    dest='motif',
                    required=True,
                    help='List of motifs. A file with only one header')
parser.add_argument(
    '-k',
    '--knownMotifs',
    dest='knownMotifs',
    required=True,
    help=
    'A tab delimited list of motifs with the following headers. MotifName regularExpression'
)

parser.add_argument('-p',
                    '--pygrStation',
                    dest='pygrStation',
                    required=True,
                    help='pygr station path. e.g., /home/data/pygr/')
#parser.add_argument('-db', '--blastdb', dest='blastdb', required=True, help='Name of the the blastdb. e.g., /u/home/t/ting0514/downloads/ncbi-blast-2.2.30+/blast_db/mm10')
parser.add_argument(
    '-g',
    '--genome',
    dest='genome',
    required=True,
    help='The UCSC genome build name, hg38, hg19, mm10, or dm3')
parser.add_argument('-o',
                    '--output',
                    dest='output',
                    required=True,
                    help='output directory')
parser.add_argument('-r',
                    '--rMATS',
                    dest='rMATS',
                    required=True,
                    help='an rMATS output file from SE event')
parser.add_argument('-mi',
                    '--miso',
                    dest='miso',
                    required=True,
                    help='an miso output file from SE event')
start = timeit.default_timer()

### these two should be mutually exclusive ####
#mxeGroup = parser.add_mutually_exclusive_group();
#mxeGroup.add_argument('-r', '--rMATS', dest='rMATS', required=True, help='an rMATS output file from SE event')
#mxeGroup.add_argument('-mi', '--miso', dest='miso', required=True, help='an miso output file from SE event')
parser.add_argument(
    '-u',
    '--up',
    dest='up',
    required=True,
    help=
    'a tab delimited file containing upregulated exons with the following headers. GeneID geneSymbol chr strand exonStart_0base exonEnd upstreamES upstreamEE downstreamES downstreamEE'
)
parser.add_argument(
    '-d',
    '--down',
    dest='dn',
    required=True,
    help=
    'a tab delimited file containing downregulated exons with the following headers. GeneIDgene Symbol chr strand exonStart_0base exonEnd upstreamES upstreamEE downstreamES downstreamEE'
)
parser.add_argument(
    '-b',
    '--background',
    dest='bg',
    required=True,
    help=
    'a tab delimited file containing background exons with the following headers. GeneID geneSymbol chr strand exonStart_0base exonEnd upstreamES upstreamEE downstreamES downstreamEE'
)
#
parser.add_argument('--label',
                    type=str,
                    dest='label',
                    default="RBP",
                    help='Label of motif. e.g., TG-rich')
parser.add_argument(
    '--intron',
    type=int,
    dest='intron',
    default=250,
    help='grab sequence up to a this number of NTs away from the exon junction'
)
parser.add_argument(
    '--exon',
    type=int,
    dest='exon',
    default=50,
    help=
    'grab sequence up to this number of NTs from both ends of the exon body')
parser.add_argument(
    '--window',
    type=int,
    dest='window',
    default=50,
    help='number of NTs examined for the given motif at a time')
parser.add_argument('--step',
                    type=int,
                    dest='step',
                    default=1,
                    help='slide window by this number of NTs at a time')
parser.add_argument('--sigFDR',
                    type=float,
                    dest='sigFDR',
                    default=0.05,
                    help='FDR cutoff for significant events.')
parser.add_argument(
    '--sigDeltaPSI',
    type=float,
    dest='sigDeltaPSI',
    default=0.05,
    help='inclusion level difference cutoff for significant events.')
parser.add_argument('--separate',
                    dest='separate',
                    default=False,
                    action='store_const',
                    const=True)
#parser.add_argument('--bgFDR', type=float, dest='sigFDR', default=0.5, help='FDR cutoff for background events.')
#parser.add_argument('-5', '--5prime', type=int, dest='fiveprime', default=6, help='do not use sequence in the first x number of nucleotides from the 5 prime splice junction')
#parser.add_argument('-3', '--3prime', type=int, dest='threeprime', default=20, help='do not use sequence in the first x number of nucleotides from the 3 prime splice junction')
args = parser.parse_args()
#


def listToString(x):  ## log command
    rVal = ''
    for a in x:
        rVal += a + ' '
    return rVal


## OUTPUT, create an output folder
outDir = args.output
os.makedirs(outDir, exist_ok=True)
outPath = os.path.abspath(outDir)
## absolute output path
#
exonDir = args.output + '/exon'
os.makedirs(exonDir, exist_ok=True)
exonPath = os.path.abspath(exonDir)
## absolute exon path (storing coord)
#
fastaDir = args.output + '/fasta'
os.makedirs(fastaDir, exist_ok=True)
fastaPath = os.path.abspath(fastaDir)
## absolute fasta path (21 fasta files will be in here)
#
mapsDir = args.output + '/maps'
os.makedirs(mapsDir, exist_ok=True)
mapsPath = os.path.abspath(mapsDir)
## absolute maps path (motif maps will be in here)
#
tempDir = args.output + '/temp'
os.makedirs(tempDir, exist_ok=True)
tempPath = os.path.abspath(tempDir)
## absolute temp path
#
scriptPath = os.path.abspath(os.path.dirname(__file__))
## absolute script path
binPath = scriptPath + '/bin'
## absolute bin path

#
VER = "2.0.1"
#
### setting up the logging format
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s %(message)s',
    #filename=outPath+'/log.motifMap.'+str(datetime.datetime.now())+'.txt',
    filename=outPath + '/log.motifMap' + '.txt',
    filemode='w')

##### Getting Start Time ######
logging.debug('Start the program with [%s]\n', listToString(sys.argv))
startTime = time.time()
logging.debug('motifTools version: %s', VER)

#
############### global variables ###############
#
pygr_station_PATH = args.pygrStation
## pygr station path

#FIVE_PRIME_LENGTH = args.fiveprime  # defines the 5' splice site length (Do not want get splice site sequences)
#THREE_PRIMER_LENGTH = args.threeprime  # defines the 3' splice site length (Do not grab splice site sequences)
rMATS = args.rMATS
miso = args.miso
## one of them are NA
up = args.up
dn = args.dn
bg = args.bg
cdist_up = {}
cdist_dn = {}
cdist_bg = {}
## count distribution for wilcoxon rank-sum test
global mFile
#
### need to have either rMATS or coordinates
#
if rMATS == "NA" and miso == "NA" and (up == "NA" or dn == "NA"
                                       or bg == "NA"):  ## not going to work
    print("Incorrect Input! Need to have rMATS input or miso input or coordinates for all three types of exons")
    sys.exit(-99)

#
iLen = args.intron  # intron length to examine
eLen = args.exon
## exon length to examine
wLen = args.window
## window length
sLen = args.step
## step size
#

sigFDR = args.sigFDR
## fdr must be smaller than this to be significant
sigDeltaPSI = args.sigDeltaPSI
## abs(deltaPSI) must be greater than or equal to this to be significant
bgFDR = 0.5
## fdr must be greater than this to be background
#bgFPKM=5.0; ## FPKM must be greater than this to be background
minPSI = 0.85
maxPSI = 0.15
u = {}
d = {}
b = {}

#uFile=null; dFile=null; bFile=null;
#if rMATS=='NA': ## we got exon coords
#  uFile=open(up);dFile=open(dn);bFile=open(bg);
if args.motif == 'NA':  ## user didn't give motif
    pass
else:  ## optional motif was given
    mFile = open(args.motif)
    ## motif(s)
kFile = open(args.knownMotifs)
## known motifs
#
region = [
    "UpstreamExon", "UpstreamExonIntron", "UpstreamIntron", "TargetExon",
    "DownstreamIntron", "DownstreamExonIntron", "DownstreamExon"
]
rName = [
    "UpstreamExon", "UpstreamExonIntron", "UpstreamIntron", "TargetExon",
    "TargetExon", "DownstreamIntron", "DownstreamExonIntron", "DownstreamExon"
]
pRegion = ["UpstreamExonIntron", "DownstreamIntron", "DownstreamExon"]
nRegion = ["UpstreamExon", "UpstreamIntron", "DownstreamExonIntron"]
#mRegion=['r1','r2','r3','r4','r5','r6','r7','r8'];
#

## motif counts in here
#
totalExonCount = {}
#
motifLabel = args.label
boxHeight = 1.0
#
uNum = 0
dNum = 0
bNum = 0
## dictionary for MAP plotting. They will contain combined counts per window
#
################################################
#
### miso conversion here
#
if miso != "NA":  ## got miso input here
    rMATS = tempPath + '/converted.rMATS.se.txt'
    ## output from conversion script, written in Perl
    convcmd = ['perl', binPath + '/miso2rMATS.SE.pl', '1', '100', miso, rMATS]
    status, output = run_command(convcmd)
    logging.debug("Conversion from miso to rMATS is done. Status: %s" % status)
    if (int(status) != 0):  ## it did not go well
        logging.debug("error in converting miso file")
        logging.debug("error detail: %s" % output)
        sys.exit(-101)
    logging.debug(output)


#
### functions here.. ##########
#
def fetch_seq(genome, strand, chr, start, end):
    """
    Fetch sequence from genome using the shared pyfaidx-based implementation.
    """
    return fetch_seq_from_fasta(genome, strand, chr, start, end)


def makeInputFiles(
    rmats
):  ## make input coords for pygr. Don't need to this if rMATS was "NA"
    global up, dn, bg
    global totalExonCount
    global nu, nd, nb
    logging.debug("Making input files from rMATS output")
    if rmats == "NA":  ## do not need to make new up,dn,bg. Just copy it over
        status, output = copy_file(up, exonPath + '/up.coord.txt')
        logging.debug("Copying up file is done. Status: %s" % status)
        if (int(status) != 0):  ## it did not go well
            logging.debug("error in copying up file")
            logging.debug("error detail: %s" % output)
            raise Exception()
        logging.debug(output)

        status, output = copy_file(dn, exonPath + '/dn.coord.txt')
        logging.debug("Copying dn file is done. Status: %s" % status)
        if (int(status) != 0):  ## it did not go well
            logging.debug("error in copying dn file")
            logging.debug("error detail: %s" % output)
            raise Exception()
        logging.debug(output)

        status, output = copy_file(bg, exonPath + '/bg.coord.txt')
        logging.debug("Copying bg file is done. Status: %s" % status)
        if (int(status) != 0):  ## it did not go well
            logging.debug("error in copying bg file")
            logging.debug("error detail: %s" % output)
            raise Exception()
        logging.debug(output)

        nb = wccount(exonPath + '/bg.coord.txt') - 1
        ## number of exons in Background without a header
        nd = wccount(exonPath + '/dn.coord.txt') - 1
        ## number of exons in Downregulated without a header
        nu = wccount(exonPath + '/up.coord.txt') - 1
        ## number of exons in Upregulated without a header
        totalExonCount = {'up': nu, 'dn': nd, 'bg': nb}
        ## total counts

        return
        ## it's done. Do not need to do more..

    rFile = open(rmats, 'r')

    uFile = open(exonPath + '/up.coord.txt', 'w')
    dFile = open(exonPath + '/dn.coord.txt', 'w')
    bFile = open(exonPath + '/bg.coord.txt', 'w')

    header = 'chr\tstrand\texonStart\texonEnd\tfirstExonStart\tfirstExonEnd\tsecondExonStart\tsecondExonEnd'
    uFile.write(header + '\n')
    dFile.write(header + '\n')
    bFile.write(header + '\n')

    line = rFile.readline()
    for line in rFile:  ## process each line

        ele = line.strip().split('\t')
        key = ':'.join(ele[3:7])
        ###print key; sys.exit();
        value = [1, '\t'.join(ele[3:11])]
        fdr = float(ele[-4])
        if fdr < sigFDR:  ## it could be significant
            deltaPSI = float(ele[-1])
            if deltaPSI >= sigDeltaPSI:  ## it's upregulated. high in sample 1
                u[key] = value
            elif deltaPSI <= -sigDeltaPSI:  ## it's downregulated. high in sample 2
                d[key] = value
        #elif fdr>bgFDR and min(float(ele[-6]),float(ele[-5]))>bgFPKM: ## it could be background
        elif fdr > bgFDR:  ## it could be background
            psi = ele[-3].replace('"', '').split(',')
            t = 0
            sum = 0.0
            PSI_1 = 0.0
            for p in psi:  ##
                if p != "NA":  ## number here
                    t += 1
                    sum += float(p)
            if t > 0:
                PSI_1 = sum / t

            psi = ele[-2].replace('"', '').split(',')
            t = 0
            sum = 0.0
            PSI_2 = 0.0
            for p in psi:  ##
                if p != "NA":  ## number here
                    t += 1
                    sum += float(p)
            if t > 0:
                PSI_2 = sum / t

            if min(PSI_1, PSI_2) < minPSI and max(
                    PSI_1, PSI_2) > maxPSI:  ## it is a background event
                b[key] = value

    logging.debug(
        "Done populating initial dictionaries with possible duplicates")
    logging.debug("Number of up, down, and background exons are: %d, %d, %d" %
                  (len(u), len(d), len(b)))

    logging.debug("Removing exons included in more than one dictionaries..")

    for key in u:  ## going through u
        if key in d:  ## same exon in d
            d[key][0] += 1
            u[key][0] += 1
        if key in b:  ## same exon in b
            b[key][0] += 1
            u[key][0] += 1

    for key in d:  ## going through d
        if key in b:  ## same exon in b
            b[key][0] += 1
            d[key][0] += 1

    ### now write to up,down,background only if count value is 1
    nu = 0
    nd = 0
    nb = 0
    ## number of up,down,background
    for key in u:
        if u[key][0] == 1:  ## it is unique
            nu += 1
            uFile.write(u[key][1] + '\n')
    for key in d:
        if d[key][0] == 1:  ## it is unique
            nd += 1
            dFile.write(d[key][1] + '\n')
    for key in b:
        if b[key][0] == 1:  ## it is unique
            nb += 1
            bFile.write(b[key][1] + '\n')

    logging.debug("Number of upregulated (high in sample_1) exons: %d" % nu)
    logging.debug("Number of downregulated (high in sample_2) exons: %d" % nd)
    logging.debug("Number of background exons: %d" % nb)

    rFile.close()
    uFile.close()
    dFile.close()
    bFile.close()

    totalExonCount = {'up': nu, 'dn': nd, 'bg': nb}
    ## total counts

    logging.debug("Done making input file from rMATS")


def getFasta(t):  ## get fasta for the given exon group

    logging.debug("Making fasta files for: %s" % t)
    tFile = []
    for rg in region:  ## region name
        tFile.append(open(fastaPath + '/' + t + '.' + rg + '.fasta', 'w'))

    f = open(exonPath + '/' + t + '.coord.txt')
    header = f.readline()
    ## don't need one
    for line in f:  ## process each event
        r = []
        ## coordinate of each region
        ele = line.strip().split('\t')
        chr = ele[0]
        strand = '+'
        r1 = [int(ele[4]), int(ele[5])]
        r2 = [int(ele[5]), int(ele[5]) + iLen + wLen]
        r3 = [int(ele[2]) - iLen - wLen, int(ele[2])]
        r4 = [int(ele[2]), int(ele[3])]
        r5 = [int(ele[3]), int(ele[3]) + iLen + wLen]
        r6 = [int(ele[6]) - iLen - wLen, int(ele[6])]
        r7 = [int(ele[6]), int(ele[7])]

        if ele[1] == '-' or ele[1] == '-1':  ## negative strand
            r7 = [int(ele[4]), int(ele[5])]
            r6 = [int(ele[5]), int(ele[5]) + iLen + wLen]
            r5 = [int(ele[2]) - iLen - wLen, int(ele[2])]
            r4 = [int(ele[2]), int(ele[3])]
            r3 = [int(ele[3]), int(ele[3]) + iLen + wLen]
            r2 = [int(ele[6]) - iLen - wLen, int(ele[6])]
            r1 = [int(ele[6]), int(ele[7])]
            strand = '-'

        r = [r1, r2, r3, r4, r5, r6, r7]
        ##    if t=='bg':
        ##      print r;
        fastaID = '>' + ':'.join(ele)
        ##    if t=='bg':
        ##      print fastaID;

        for i in range(len(r)):  ### start writing each region
            fastaSeq = fetch_seq(genome, strand, chr, r[i][0], r[i][1])
            tFile[i].write(fastaID + '\n' + fastaSeq + '\n')


##    if t=='bg':
##      print fastaID;
    logging.debug("Done making fasta files for: %s" % t)


def old_findAll(sub, a_str):
    start = 0
    while True:
        start = a_str.find(sub, start)
        if start == -1: return
        yield start
        start += 1


def findAll(re_motif, s_seq):
    return [[m.start(0), m.end(0)] for m in re.finditer(re_motif, s_seq)
            ], [m.start(0) for m in re.finditer(re_motif, s_seq)]


def initMotif(mF, mc):  ## initialize motif counts

    motifs = []
    logging.debug("Initializing motif counts")
    header = mF.readline()
    ## don't use it
    for line in mF:  ## process each motif
        motifs.append(line.strip().split('\t')[1])
    logging.debug("The number of motifs: %d" % len(motifs))
    return motifs


def countMotif(mc, ttName, ttMotif,
               motifs):  ## counting motifs, motif name, motif reg expression
    global totalExonCount
    global uNum, dNum, bNum
    global nu, nd, nb
    global iLen, eLen, wLen, sLen, region

    cdist = {'up': {}, 'dn': {}, 'bg': {}}
    eP = eLen // sLen
    ## exon points
    iP = iLen // sLen
    ##  intron points

    ###### init cdist
    for kkkk in ['up', 'dn', 'bg']:
        for iiii in range(8):
            cdist[kkkk][iiii] = {}

    for kkkk in ['up', 'dn', 'bg']:
        for epep in range(eP):  ## for 0,3,4,7
            for epindex in [0, 3, 4, 7]:  ## that has eP points
                cdist[kkkk][epindex][epep] = [0] * totalExonCount[kkkk]
        for ipip in range(iP):  ## for 1,2,5,6
            for ipindex in [1, 2, 5, 6]:  ## that has iP points
                cdist[kkkk][ipindex][ipip] = [0] * totalExonCount[kkkk]

    ##### check each of 8 motif regions now
    tempPosition = []
    for jj in ['up', 'dn', 'bg']:
        for rg in region:
            fFile = open(fastaPath + '/' + jj + '.' + rg + '.fasta')
            exonNum = 0
            for dummy in fFile:
                seq = next(fFile).strip()
                seqLen = len(seq)
                for mm in motifs:
                    dd, cc = findAll(mm, seq)
                    if len(cc) > 0:  ## mapped somewhere
                        mLen = dd[0][1] - dd[0][0]
                        for pInd in cc:  ## index
                            nInd = pInd + mLen - 1 - seqLen

                            if pInd >= (max(iLen, eLen) + wLen) or nInd < -(
                                    max(iLen, eLen) + wLen
                            ):  ## index out of range, due to the exon body
                                continue
                            else:  ## incread count for both mc and countdist

                                ### get position, add 1 at the exonNum index
                                enInd = nInd + eLen + wLen - 1
                                ## convert nInd to pInd form for exon
                                inInd = nInd + iLen + wLen - 1
                                ## convert nInd to pInd form for intron
                                if rg == "UpstreamExon":  ## do it negative way
                                    #tempPosition=range(eP+max(-eLen, nInd-wLen+1)// sLen,min(eP,nInd// sLen)-1);
                                    if nInd < (-eLen - wLen):
                                        continue
                                    tempPosition = range(
                                        max(0, enInd - wLen + 1) // sLen,
                                        min(eP, enInd // sLen + 1))
                                    ## print tempPosition; print nInd,eP,sLen,eLen,wLen;
                                    for tptp in tempPosition:
                                        cdist[jj][0][tptp][exonNum] += 1
                                elif rg == "UpstreamExonIntron":  ### do it positive way
                                    if pInd >= (iLen + wLen):
                                        continue
                                    tempPosition = range(
                                        max(0, pInd - wLen + 1) // sLen,
                                        min(iP, pInd // sLen) - 1)
                                    for tptp in tempPosition:
                                        cdist[jj][1][tptp][exonNum] += 1
                                elif rg == "UpstreamIntron":  ## do it negative way
                                    if nInd < (-iLen - wLen):
                                        continue
                                    tempPosition = range(
                                        max(0, inInd - wLen + 1) // sLen,
                                        min(iP, inInd // sLen + 1))
                                    for tptp in tempPosition:
                                        cdist[jj][2][tptp][exonNum] += 1
                                elif rg == "TargetExon":  ### do it both ways
                                    if (pInd < (eLen + wLen)):  ## valid
                                        tempPosition = range(
                                            max(0, pInd - wLen + 1) // sLen,
                                            min(eP, pInd // sLen) - 1)
                                        ## postive way
                                        for tptp in tempPosition:
                                            cdist[jj][3][tptp][exonNum] += 1
                                    if nInd >= (-eLen - wLen):  ## valid
                                        tempPosition = range(
                                            max(0, enInd - wLen + 1) // sLen,
                                            min(eP, enInd // sLen + 1))
                                        ## negative way
                                        for tptp in tempPosition:
                                            cdist[jj][4][tptp][exonNum] += 1
                                elif rg == "DownstreamIntron":  ### do it positive way
                                    if pInd >= (iLen + wLen):
                                        continue
                                    tempPosition = range(
                                        max(0, pInd - wLen + 1) // sLen,
                                        min(iP, pInd // sLen) - 1)
                                    for tptp in tempPosition:
                                        cdist[jj][5][tptp][exonNum] += 1
                                elif rg == "DownstreamExonIntron":  ## do it negative way
                                    if nInd < (-iLen - wLen):
                                        continue
                                    tempPosition = range(
                                        max(0, inInd - wLen + 1) // sLen,
                                        min(iP, inInd // sLen + 1))
                                    for tptp in tempPosition:
                                        cdist[jj][6][tptp][exonNum] += 1
                                elif rg == "DownstreamExon":  ### do it positive way
                                    if pInd >= (eLen + wLen):
                                        continue
                                    tempPosition = range(
                                        max(0, pInd - wLen + 1) // sLen,
                                        min(eP, pInd // sLen) - 1)
                                    for tptp in tempPosition:
                                        cdist[jj][7][tptp][exonNum] += 1
                #print exonNum;
                exonNum += 1
            fFile.close()
    return cdist


def drawNode(c, eH, eW, iW, indent, gap, sGap, scale):  ## draw node

    mls = style.linewidth.THIck
    ## my line style

    upstreamExon = path.path(path.moveto(0, 0), path.lineto(eW * scale, 0),
                             path.lineto(scale * eW, scale * eH),
                             path.lineto(0, scale * eH),
                             path.lineto(scale * indent, scale * eH / 2),
                             path.lineto(0, 0))
    intron_1 = path.line(scale * eW, scale * eH / 2, scale * (eW + iW),
                         scale * eH / 2)
    divider_1 = [
        path.line((eW + iW - sGap) * scale, (eH / 2 - gap) * scale,
                  (eW + iW + sGap) * scale, (eH / 2 + gap) * scale),
        path.line((eW + iW - sGap + gap) * scale, (eH / 2 - gap) * scale,
                  (eW + iW + sGap + gap) * scale, (eH / 2 + gap) * scale)
    ]
    intron_2 = path.line((eW + iW + gap) * scale, scale * eH / 2,
                         (eW + iW + gap + iW) * scale, scale * eH / 2)

    targetExon = path.rect((eW + iW + gap + iW) * scale, 0,
                           (eW + gap + gap + eW) * scale, eH * scale)

    newX = eW + iW + gap + iW + eW + gap + gap + eW
    ## new start of x
    intron_3 = path.line(newX * scale, (eH / 2) * scale, (newX + iW) * scale,
                         scale * eH / 2)
    divider_2 = [
        path.line((newX + iW - sGap) * scale, (eH / 2 - gap) * scale,
                  (newX + iW + sGap) * scale, (eH / 2 + gap) * scale),
        path.line((newX + iW - sGap + gap) * scale, (eH / 2 - gap) * scale,
                  (newX + iW + sGap + gap) * scale, (eH / 2 + gap) * scale)
    ]
    intron_4 = path.line((newX + iW + gap) * scale, scale * eH / 2,
                         (newX + iW + gap + iW) * scale, (eH / 2) * scale)
    newX = newX + iW + gap + iW
    ## newX
    downstreamExon = path.path(
        path.moveto(newX * scale, 0), path.lineto((newX + eW) * scale, 0),
        path.lineto((newX + eW - indent) * scale, scale * eH / 2),
        path.lineto((newX + eW) * scale, eH * scale),
        path.lineto(newX * scale, eH * scale), path.lineto(newX * scale, 0))

    c.stroke(upstreamExon, [mls, deco.filled([color.gray(0.8)])])
    c.stroke(intron_1, [mls])
    c.stroke(divider_1[0], [mls])
    c.stroke(divider_1[1], [mls])
    c.stroke(intron_2, [mls])
    c.stroke(targetExon, [mls, deco.filled([color.rgb.green])])
    c.stroke(intron_3, [mls])
    c.stroke(divider_2[0], [mls])
    c.stroke(divider_2[1], [mls])
    c.stroke(intron_4, [mls])
    c.stroke(downstreamExon, [mls, deco.filled([color.gray(0.8)])])

    logging.debug("Done drawing node structure")


def drawBox(c,
            max_point_value,
            exon_height,
            exon_width,
            intron_width,
            indent,
            divider_gap,
            slope_gap,
            scale,
            exon_length,
            intron_length,
            map_name,
            maxNegP,
            separate=False):  ## draw plot box
    #global motifLabel;
    global boxHeight

    mls = style.linewidth.THIck
    ## my line style
    ldash = style.linestyle.dashed
    largeLab = text.size.huge
    yLabAtt = [largeLab, text.halign.boxright, text.valign.middle]
    ypLabAtt = [largeLab, text.halign.boxleft, text.valign.middle]

    boxY = exon_height + exon_height
    boxHeight = intron_width * 1.2
    bW = exon_width + intron_width

    second_box_offset = -(boxHeight + 2 * boxY)

    myMax = max_point_value * 1.01
    yLab = drawutils.make_label(myMax)
    ypLab = drawutils.make_label(maxNegP)

    units_per_bp = float(exon_width) / exon_length

    xLab_1_3 = [-exon_length, 0, intron_length / 2, intron_length]
    ## for 1st and 3rd boxes
    xLab_2_4 = [-intron_length, -intron_length / 2, 0, exon_length]
    ## for the 2nd and 4th boxes

    xS = 0
    xE = xS + bW
    rect1, r1_line, r1_ss = drawutils.boxes(xS, bW, scale, boxY, boxHeight,
                                            exon_width)

    ### print labels

    # y-axis
    drawutils.draw_y_axis(c, scale, yLab, xS, boxY, boxHeight)

    x_axis_canvas = canvas.canvas()

    #x-axis
    drawutils.draw_x_axis_segment(x_axis_canvas, scale, xS, boxY, xLab_1_3, units_per_bp)

    xS = xE + divider_gap
    xE = xS + bW
    rect2, r2_line, r2_ss = drawutils.boxes(xS, bW, scale, boxY, boxHeight,
                                            intron_width)

    #x-axis
    drawutils.draw_x_axis_segment(x_axis_canvas, scale, xS, boxY, xLab_2_4, units_per_bp)

    ## Title and legends
    drawutils.title_and_legend(c, scale, xE, divider_gap, boxY, boxHeight,
                               exon_width, intron_width, map_name, nu, nd, nb)

    xS = xE + divider_gap + divider_gap
    xE = xS + bW
    rect3, r3_line, r3_ss = drawutils.boxes(xS, bW, scale, boxY, boxHeight,
                                            exon_width)

    #x-axis
    drawutils.draw_x_axis_segment(x_axis_canvas, scale, xS, boxY, xLab_1_3, units_per_bp)

    xS = xE + divider_gap
    xE = xS + bW
    rect4, r4_line, r4_ss = drawutils.boxes(xS, bW, scale, boxY, boxHeight,
                                            intron_width)

    #x-axis
    drawutils.draw_x_axis_segment(x_axis_canvas, scale, xS, boxY, xLab_2_4, units_per_bp)

    c.insert(x_axis_canvas)
    if separate:
        c.insert(x_axis_canvas, [trafo.translate(0, -boxY * scale)])

    box_canvas = canvas.canvas()

    box_canvas.stroke(rect1, [mls])
    box_canvas.stroke(rect2, [mls])
    box_canvas.stroke(rect3, [mls])
    box_canvas.stroke(rect4, [mls])

    for i in range(3):
        box_canvas.stroke(r1_line[i], [ldash])
        box_canvas.stroke(r2_line[i], [ldash])
        box_canvas.stroke(r3_line[i], [ldash])
        box_canvas.stroke(r4_line[i], [ldash])

    ## lines for splice sites
    box_canvas.stroke(r1_ss)
    box_canvas.stroke(r2_ss)
    box_canvas.stroke(r3_ss)
    box_canvas.stroke(r4_ss)

    c.insert(box_canvas)
    if separate:
        c.insert(box_canvas, [trafo.translate(0, second_box_offset * scale)])

    ### print pvalue labels
    # yp-axis
    yp_y_offset = second_box_offset if separate else 0
    xE = xE + exon_width
    drawutils.draw_yp_axis(c, scale, ypLab, yp_y_offset, xE, boxY, boxHeight,
                           divider_gap)

    logging.debug("Done drawing plot box")


def plotRegions(
    tRegion_1,
    tRegion_2,
    xS,
    boxY,
    bH,
    fW,
    sW,
    mpv,
    scale,
    howmany=3
):  ## plot given regions, fW, sW tells the width of the first and the second region
    pup = []
    pdn = []
    pbg = []
    ## path for up, dn, and bg
    pvalup = []
    pvaldn = []
    ## path for pval up, pval dn
    y1 = []
    y2 = []

    for jj in range(len(tRegion_1) -
                    1):  ## for each point in the first region (upstream exon)
        x1 = xS + jj * float(fW) / len(tRegion_1) + (float(fW) /
                                                     len(tRegion_1)) / 2.0
        x2 = xS + (jj + 1) * float(fW) / len(tRegion_1) + (
            float(fW) / len(tRegion_1)) / 2.0
        if howmany == 3:  ## actual points
            y1 = [
                boxY + tRegion_1[jj][0] * bH / mpv,
                boxY + tRegion_1[jj][1] * bH / mpv,
                boxY + tRegion_1[jj][2] * bH / mpv
            ]
            y2 = [
                boxY + tRegion_1[jj + 1][0] * bH / mpv,
                boxY + tRegion_1[jj + 1][1] * bH / mpv,
                boxY + tRegion_1[jj + 1][2] * bH / mpv
            ]

            pup.append([x1 * scale, y1[0] * scale])
            pdn.append([x1 * scale, y1[1] * scale])
            pbg.append([x1 * scale, y1[2] * scale])

            if jj == len(
                    tRegion_1
            ) - 2:  ## second to the last element. add coords for the last element
                pup.append([x2 * scale, y2[0] * scale])
                pdn.append([x2 * scale, y2[1] * scale])
                pbg.append([x2 * scale, y2[2] * scale])

        elif howmany == 2:  ## pvalue points
            y1 = [
                boxY + tRegion_1[jj][0] * bH / mpv,
                boxY + tRegion_1[jj][1] * bH / mpv
            ]
            y2 = [
                boxY + tRegion_1[jj + 1][0] * bH / mpv,
                boxY + tRegion_1[jj + 1][1] * bH / mpv
            ]

            pvalup.append([x1 * scale, y1[0] * scale])
            pvaldn.append([x1 * scale, y1[1] * scale])

            if jj == len(
                    tRegion_1
            ) - 2:  ## second to the last element. add coords for the last element
                pvalup.append([x2 * scale, y2[0] * scale])
                pvaldn.append([x2 * scale, y2[1] * scale])

    xS = xS + fW
    for jj in range(len(tRegion_2) -
                    1):  ## for each point in the second region (intron)
        x1 = xS + jj * float(sW) / len(tRegion_2) + (float(sW) /
                                                     len(tRegion_2)) / 2.0
        x2 = xS + (jj + 1) * float(sW) / len(tRegion_2) + (
            float(sW) / len(tRegion_2)) / 2.0
        if howmany == 3:  ## actual points
            y1 = [
                boxY + tRegion_2[jj][0] * bH / mpv,
                boxY + tRegion_2[jj][1] * bH / mpv,
                boxY + tRegion_2[jj][2] * bH / mpv
            ]
            y2 = [
                boxY + tRegion_2[jj + 1][0] * bH / mpv,
                boxY + tRegion_2[jj + 1][1] * bH / mpv,
                boxY + tRegion_2[jj + 1][2] * bH / mpv
            ]

            pup.append([x1 * scale, y1[0] * scale])
            pdn.append([x1 * scale, y1[1] * scale])
            pbg.append([x1 * scale, y1[2] * scale])

            if jj == len(
                    tRegion_2
            ) - 2:  ## second to the last element. add coords for the last element
                pup.append([x2 * scale, y2[0] * scale])
                pdn.append([x2 * scale, y2[1] * scale])
                pbg.append([x2 * scale, y2[2] * scale])
        elif howmany == 2:  ## pvalue points
            y1 = [
                boxY + tRegion_2[jj][0] * bH / mpv,
                boxY + tRegion_2[jj][1] * bH / mpv
            ]
            y2 = [
                boxY + tRegion_2[jj + 1][0] * bH / mpv,
                boxY + tRegion_2[jj + 1][1] * bH / mpv
            ]

            pvalup.append([x1 * scale, y1[0] * scale])
            pvaldn.append([x1 * scale, y1[1] * scale])

            if jj == len(
                    tRegion_2
            ) - 2:  ## second to the last element. add coords for the last element
                pvalup.append([x2 * scale, y2[0] * scale])
                pvaldn.append([x2 * scale, y2[1] * scale])

    logging.debug("Done plotRegions function")
    if howmany == 3:
        return pup, pdn, pbg
    elif howmany == 2:
        return pvalup, pvaldn


#


def fillUpPath(tPoints):
    rPath = path.path(path.moveto(tPoints[0][0], tPoints[0][1]))
    for pp in range(1, len(tPoints)):  ## for each point from the 2nd element
        rPath.append(path.lineto(tPoints[pp][0], tPoints[pp][1]))
    return rPath
    logging.debug("Done fillUpPath function")


def drawAcutalPlot(
    c,
    dp,
    mpv,
    eH,
    eW,
    iW,
    gap,
    sGap,
    scale,
    wl,
    sl,
    nPP,
    mNP,
    separate=False
):  ## draw plots now.. nPP is negativePvaluePoins, mNP is maximumNegPval
    global boxHeight

    mls = style.linewidth.THIck
    ## my line style
    bgmls = style.linewidth.Thick
    ## my line style
    ldash = style.linestyle.dashed
    largeLab = text.size.huge

    boxY = eH + eH
    bH = boxHeight
    bW = eW + iW
    second_box_offset = -(bH + 2 * boxY)
    p_transform = trafo.translate(0,
                                  second_box_offset * scale if separate else 0)

    upColor = color.rgb.red
    dnColor = color.rgb.blue
    bgColor = color.rgb.black

    ## 1st and 2nd region
    xS = 0
    points_up, points_dn, points_bg = plotRegions(dp[0], dp[1], xS, boxY, bH,
                                                  eW, iW, mpv, scale, 3)
    ## return 3 values
    path_up = fillUpPath(points_up)
    path_dn = fillUpPath(points_dn)
    path_bg = fillUpPath(points_bg)
    c.stroke(path_up, [upColor, mls])
    c.stroke(path_dn, [dnColor, mls])
    c.stroke(path_bg, [bgColor, bgmls])
    points_pup, points_pdn = plotRegions(nPP[0], nPP[1], xS, boxY, bH, eW, iW,
                                         mNP, scale, 2)
    ## return 2 values
    path_pup = fillUpPath(points_pup)
    path_pdn = fillUpPath(points_pdn)
    c.stroke(path_pup, [upColor, ldash, p_transform])
    c.stroke(path_pdn, [dnColor, ldash, p_transform])

    ## 3rd and 4th region
    xS = xS + eW + iW + gap
    points_up, points_dn, points_bg = plotRegions(dp[2], dp[3], xS, boxY, bH,
                                                  iW, eW, mpv, scale)
    path_up = fillUpPath(points_up)
    path_dn = fillUpPath(points_dn)
    path_bg = fillUpPath(points_bg)
    c.stroke(path_up, [upColor, mls])
    c.stroke(path_dn, [dnColor, mls])
    c.stroke(path_bg, [bgColor, bgmls])
    points_pup, points_pdn = plotRegions(nPP[2], nPP[3], xS, boxY, bH, iW, eW,
                                         mNP, scale, 2)
    ## return 2 values
    path_pup = fillUpPath(points_pup)
    path_pdn = fillUpPath(points_pdn)
    c.stroke(path_pup, [upColor, ldash, p_transform])
    c.stroke(path_pdn, [dnColor, ldash, p_transform])

    ## 5th and 6th region
    xS = xS + iW + eW + gap + gap
    points_up, points_dn, points_bg = plotRegions(dp[4], dp[5], xS, boxY, bH,
                                                  eW, iW, mpv, scale)
    path_up = fillUpPath(points_up)
    path_dn = fillUpPath(points_dn)
    path_bg = fillUpPath(points_bg)
    c.stroke(path_up, [upColor, mls])
    c.stroke(path_dn, [dnColor, mls])
    c.stroke(path_bg, [bgColor, bgmls])
    points_pup, points_pdn = plotRegions(nPP[4], nPP[5], xS, boxY, bH, eW, iW,
                                         mNP, scale, 2)
    ## return 2 values
    path_pup = fillUpPath(points_pup)
    path_pdn = fillUpPath(points_pdn)
    c.stroke(path_pup, [upColor, ldash, p_transform])
    c.stroke(path_pdn, [dnColor, ldash, p_transform])

    ## 7th and8 th region
    xS = xS + eW + iW + gap
    points_up, points_dn, points_bg = plotRegions(dp[6], dp[7], xS, boxY, bH,
                                                  iW, eW, mpv, scale)
    path_up = fillUpPath(points_up)
    path_dn = fillUpPath(points_dn)
    path_bg = fillUpPath(points_bg)
    c.stroke(path_up, [upColor, mls])
    c.stroke(path_dn, [dnColor, mls])
    c.stroke(path_bg, [bgColor, bgmls])
    points_pup, points_pdn = plotRegions(nPP[6], nPP[7], xS, boxY, bH, iW, eW,
                                         mNP, scale, 2)
    ## return 2 values
    path_pup = fillUpPath(points_pup)
    path_pdn = fillUpPath(points_pdn)
    c.stroke(path_pup, [upColor, ldash, p_transform])
    c.stroke(path_pdn, [dnColor, ldash, p_transform])

    logging.debug("Done drawing acutal plots")


def countPerWindow(cpw, ic, sign, rLen):
    global wLen, sLen
    sInd = 0
    cInd = 0
    if sign == -1:  ## need to to this in reverse way
        sInd = -rLen
        cInd = -1
    for k in range(len(cpw)):  ## do it this many times
        for mm in motifs:  ## examine all motifs
            for n in range(sInd, sInd + sign * (wLen), sign):
                cpw[k] += ic[mm][n][cInd]
        sInd += sLen


def getCombinedCounts(exon):
    global iLen, eLen, wLen, sLen, region

    #eP=1+(eLen-wLen)// sLen; ## exon points
    #iP=1+(iLen-wLen)// sLen; ##  intron points
    eP = eLen // sLen
    ## exon points
    iP = iLen // sLen
    ##  intron points

    rVal = {}
    ## count for 8 regions at each step
    rVal[0] = [0] * eP
    rVal[1] = [0] * iP
    rVal[2] = [0] * iP
    rVal[3] = [0] * eP
    rVal[4] = [0] * eP
    rVal[5] = [0] * iP
    rVal[6] = [0] * iP
    rVal[7] = [0] * eP

    countPerWindow(rVal[0], exon[region[0]], -1, eLen)
    #print exon[region[3]]['TGGTGG']
    countPerWindow(rVal[1], exon[region[1]], 1, iLen)
    countPerWindow(rVal[2], exon[region[2]], -1, iLen)
    countPerWindow(rVal[3], exon[region[3]], 1, eLen)
    countPerWindow(rVal[4], exon[region[3]], -1, eLen)
    countPerWindow(rVal[5], exon[region[4]], 1, iLen)
    countPerWindow(rVal[6], exon[region[5]], -1, iLen)
    countPerWindow(rVal[7], exon[region[6]], 1, eLen)

    logging.debug("Done getting combined motif counts")
    return rVal
    ## list with total count for each region


def ccc(ic):
    global iLen, eLen, wLen, sLen, region
    eP = eLen // sLen
    ## exon points
    iP = iLen // sLen
    rVal = {}
    ## count for 8 regions at each step
    rVal[0] = [0] * eP
    rVal[1] = [0] * iP
    rVal[2] = [0] * iP
    rVal[3] = [0] * eP
    rVal[4] = [0] * eP
    rVal[5] = [0] * iP
    rVal[6] = [0] * iP
    rVal[7] = [0] * eP
    for i in range(8):
        for k in range(len(rVal[i])):
            rVal[i][k] += sum(ic[i][k])
    return rVal


def plotMotifs_finale(mapName,
                      maxPointValue,
                      maxNegPval,
                      drawPoints,
                      negPvalPoints,
                      separate=False):

    exonHeight = 30
    exonWidth = eLen
    intronWidth = iLen
    indentation = 5
    dividerGap = 12
    slopeGap = 20
    Scale = 0.03
    canv = canvas.canvas()

    drawNode(canv, exonHeight, exonWidth, intronWidth, indentation, dividerGap,
             slopeGap, Scale)
    drawBox(canv, maxPointValue, exonHeight, exonWidth, intronWidth,
            indentation, dividerGap, slopeGap, Scale, eLen, iLen, mapName,
            maxNegPval, separate)
    drawAcutalPlot(canv, drawPoints, maxPointValue, exonHeight, exonWidth,
                   intronWidth, dividerGap, slopeGap, Scale, wLen, sLen,
                   negPvalPoints, maxNegPval, separate)

    safe_map_name = re.sub(r'[^A-Za-z0-9._-]+', '_', mapName)
    canv.writePDFfile(mapsPath + '/' + 'SE.' + safe_map_name + '.pdf')
    ## write to PDF

    try:
        canv.writeGSfile(mapsPath + '/jj.png', "png16m", resolution=100)
        os.rename(mapsPath + '/jj.png', mapsPath + '/' + safe_map_name +
                  '.png')
    except FileNotFoundError:
        logging.debug("Ghostscript not found; skipping PNG export for %s", mapName)


def plotMotifs(
    d, mc, mapName, pdic_up, pdic_dn, motifs,
    cdist):  ## plotting motifs (count, name, pvalue.up.vs.bg, pvalue.dn.vs.bg)

    global uNum, dNum, bNum
    upc = {}
    dnc = {}
    bgc = {}
    #print len(motifs)

    #upc=getCombinedCounts(mc['up']);
    #dnc=getCombinedCounts(mc['dn']);
    #bgc=getCombinedCounts(mc['bg']);
    upc = ccc(cdist['up'])
    dnc = ccc(cdist['dn'])
    bgc = ccc(cdist['bg'])

    ## setting up a canvas
    ## set canvas
    drawPoints = []
    negPvalPoints = []
    maxPointValue = 0.0000000000000000000000000000000000000000000000000001
    maxNegPval = 0.00000000000000000000000000000000000000000000000000000000000000001
    for zz in range(8):  ## for 8 regions
        drawPoints.append([])
        negPvalPoints.append([])
        for ind in range(
                len(upc[zz])
        ):  ## everything has the same items. It's okay to use this.
            drawPoints[zz].append([
                min(1.0,
                    float(upc[zz][ind]) / float(uNum * len(motifs))),
                min(1.0,
                    float(dnc[zz][ind]) / float(dNum * len(motifs))),
                min(1.0,
                    float(bgc[zz][ind]) / float(bNum * len(motifs)))
            ])
            negPvalPoints[zz].append([
                -numpy.log10(pdic_up[zz][ind]), -numpy.log10(pdic_dn[zz][ind])
            ])
            for yy in drawPoints[zz][ind]:  ## examine three values
                if yy > maxPointValue:  ## new max here
                    maxPointValue = yy
            for yy in negPvalPoints[zz][ind]:  ## examine two values
                if yy > maxNegPval:  ## new max here
                    maxNegPval = yy
    logging.debug("Max average height is: %f" % maxPointValue)
    logging.debug("Max negPval is: %f" % maxNegPval)
    val = [mapName, maxPointValue, maxNegPval, drawPoints, negPvalPoints]
    d.append(val)


def printCounts(
    counts, dFile, eNum
):  ## print counts on the destination file, eNum is the nubmer of exons in the group
    global rName
    dFile.write('Region\tPosition\tRawPeakCounts\tAverageRawPeakCount\n')
    if eNum == 0:  ## no exons in the group
        return
    for zz in range(8):  ## for 8 regions
        for ind in range(len(counts[zz])):
            dFile.write('\t'.join([
                rName[zz],
                str(ind),
                str(counts[zz][ind]),
                str(float(counts[zz][ind]) / float(eNum))
            ]) + '\n')
        #dFile.write('\n');


def printCombinedCounts(upc, dnc, bgc, destFile, uNum, dNum,
                        bNum):  ## print combined counts
    global rName
    destFile.write(
        'Region\tPosition\tRawPeakCounts_up\tRawPeakCounts_down\tRawPeakCounts_background\tAverageRawPeakCount_up\tAverageRawPeakCount_down\tAverageRawPeakCount_background\n'
    )
    if uNum * dNum * bNum == 0:  ## no exons in the group
        logging.debug(
            "No exons in one of the upregulated, downregulated, or background exon groups. No combined result available."
        )
        return
    for zz in range(8):  ## for 8 regions
        for ind in range(
                len(upc[zz])
        ):  ## everything has the same items. It's okay to use this.
            outString = '\t'.join([
                rName[zz],
                str(ind),
                str(upc[zz][ind]),
                str(dnc[zz][ind]),
                str(bgc[zz][ind])
            ]) + '\t'
            outString += '\t'.join([
                str(float(upc[zz][ind]) / float(uNum)),
                str(float(dnc[zz][ind]) / float(dNum)),
                str(float(bgc[zz][ind]) / float(bNum))
            ])
            destFile.write(outString + '\n')
        #dFile.write('\n');


def printCountDist(cdist, tName, tMotif):
    rName = {
        0: 'UpstreamExon_3prime',
        1: 'UpstreamExonIntron',
        2: 'UpstreamIntron',
        3: 'TargetExon_5prime',
        4: 'TargetExon-3prime',
        5: 'DownstreamIntron',
        6: 'DownstreamExonIntron',
        7: 'DownstreamExon_5prime'
    }
    for fff in ['up', 'dn', 'bg']:
        desFile = open(
            tempPath + '/' + tName + '.' + tMotif + '.countDist.' + fff +
            '.txt', 'w')
        desFile.write('Region\tposition\tsum\tvalues\n')
        #if exNum==0: ## no exons in the group
        #  return;
        for zz in range(8):  ## for 8 regions
            for locus in range(len(cdist[fff][zz])):
                desFile.write(rName[zz] + '\t' + str(locus) + '\t' +
                              str(sum(cdist[fff][zz][locus])) + '\t[')
                desFile.write(
                    ','.join([str(ccc)
                              for ccc in cdist[fff][zz][locus]]) + ']\n')
        desFile.close()


def computeWilcoxonP(cdist_one, cdist_two,
                     test_p):  ## count p value for one vs. two
    rName = {
        0: 'UpstreamExon_3prime',
        1: 'UpstreamExonIntron',
        2: 'UpstreamIntron',
        3: 'TargetExon_5prime',
        4: 'TargetExon-3prime',
        5: 'DownstreamIntron',
        6: 'DownstreamExonIntron',
        7: 'DownstreamExon_5prime'
    }
    for zz in range(8):  ## for 8 regions
        test_p[zz] = {}
        for locus in range(len(cdist_one[zz])):
            #test_p[zz][locus]=stats.ranksums(cdist_one[zz][locus] ,cdist_two[zz][locus])[1]; ## store p-value
            test_p[zz][locus] = stats.fisher_exact(
                [[
                    sum(cdist_one[zz][locus]),
                    max(0,
                        len(cdist_one[zz][locus]) - sum(cdist_one[zz][locus]))
                ],
                 [
                     sum(cdist_two[zz][locus]),
                     max(0,
                         len(cdist_two[zz][locus]) - sum(cdist_two[zz][locus]))
                 ]], 'greater')[1]
            ## store p-value
    return test_p


def printPval(pdic, dFile, eNum):  ## print p values per position
    rName = {
        0: 'UpstreamExon_3prime',
        1: 'UpstreamExonIntron',
        2: 'UpstreamIntron',
        3: 'TargetExon_5prime',
        4: 'TargetExon-3prime',
        5: 'DownstreamIntron',
        6: 'DownstreamExonIntron',
        7: 'DownstreamExon_5prime'
    }
    dFile.write('Region\tposition\twilcoxon.ranksum.pVal\n')
    if eNum == 0:  ## no exons in the group
        return
    for zz in range(8):  ## for 8 regions
        for locus in range(len(pdic[zz])):
            dFile.write(rName[zz] + '\t' + str(locus) + '\t' +
                        str(pdic[zz][locus]) + '\n')


def makeIndividualMaps(d, line):
    global mCount, totalExonCount

    test_up1 = {}
    test_dn1 = {}
    ## p-value storage
    ele = line.strip().split('\t')
    tName = ele[0]
    tMotif = ele[1]
    ## regular expression

    tmpFile = open(tempPath + '/' + tName + '.' + tMotif + '.txt', 'w')
    tmpFile.write(tHeader + '\n')
    tmpFile.write(tName + '\t' + tMotif + '\n')
    tmpFile.close()
    tmpFile = open(tempPath + '/' + tName + '.' + tMotif + '.txt')
    mCount = {}
    motifs = initMotif(tmpFile, mCount)
    cdist = countMotif(mCount, tName, tMotif, motifs)
    printCountDist(cdist, tName, tMotif)
    test_up = computeWilcoxonP(cdist['up'], cdist['bg'], test_up1)
    test_dn = computeWilcoxonP(cdist['dn'], cdist['bg'], test_dn1)
    upbgPFile = open(
        tempPath + '/' + tName + '.' + tMotif + '.pVal.up.vs.bg.txt', 'w')
    dnbgPFile = open(
        tempPath + '/' + tName + '.' + tMotif + '.pVal.dn.vs.bg.txt', 'w')
    printPval(test_up, upbgPFile, totalExonCount['up'])
    printPval(test_dn, dnbgPFile, totalExonCount['dn'])
    plotMotifs(d, mCount, tName + '-' + tMotif, test_up, test_dn, motifs,
               cdist)
    upbgPFile.close()
    dnbgPFile.close()
    tmpFile.close()

    ## end for line in mf


##
### Count line in a file ###
##
def wccount(filename):
    count = 0
    with open(filename, "rb") as f:
        for _ in f:
            count += 1
    return count


## end of wccount()


##
### File out of minimum pValue list
##
def minPvalueOut(exonType):
    global outPath
    global tempPath

    fileList = ""
    resultFile = ""
    if exonType == "up":
        fileList = glob.glob(tempPath + "/" + "*.pVal.up.vs.bg.txt")
        resultFile = "pVal.up.vs.bg.RNAmap.txt"
    elif exonType == "down":
        fileList = glob.glob(tempPath + "/" + "*.pVal.dn.vs.bg.txt")
        resultFile = "pVal.dn.vs.bg.RNAmap.txt"
    else:
        logging.debug("Exon type is not correct in Minimum p-Value list...")
        return
    pValList = []
    for fN in fileList:
        pValData = open(fN, "r")
        pValUpExon3pChunk = []
        pValUpExonIntronChunk = []
        pValUpIntronChunk = []
        pValTg5pChunk = []
        pValTg3pChunk = []
        pValDnIntronChunk = []
        pValDnExonIntronChunk = []
        pValDnExon5pChunk = []
        for eachLine in pValData:
            tmpStr = eachLine.strip().split('\t')
            if tmpStr[0] == "DownstreamExon_5prime":
                pValDnExon5pChunk.append(float(tmpStr[-1]))
            elif tmpStr[0] == "DownstreamExonIntron":
                pValDnExonIntronChunk.append(float(tmpStr[-1]))
            elif tmpStr[0] == "DownstreamIntron":
                pValDnIntronChunk.append(float(tmpStr[-1]))
            elif tmpStr[0] == "TargetExon-3prime":
                pValTg3pChunk.append(float(tmpStr[-1]))
            elif tmpStr[0] == "TargetExon_5prime":
                pValTg5pChunk.append(float(tmpStr[-1]))
            elif tmpStr[0] == "UpstreamIntron":
                pValUpIntronChunk.append(float(tmpStr[-1]))
            elif tmpStr[0] == "UpstreamExonIntron":
                pValUpExonIntronChunk.append(float(tmpStr[-1]))
            elif tmpStr[0] == "UpstreamExon_3prime":
                pValUpExon3pChunk.append(float(tmpStr[-1]))
            else:
                pass
        pValUpExon3pChunk.sort()
        pValUpExonIntronChunk.sort()
        pValUpIntronChunk.sort()
        pValTg5pChunk.sort()
        pValTg3pChunk.sort()
        pValDnIntronChunk.sort()
        pValDnExonIntronChunk.sort()
        pValDnExon5pChunk.sort()
        fileName = fN.strip().split('/')[-1]
        pValList.append([
            fileName[:-18], pValUpExon3pChunk[0], pValUpExonIntronChunk[0],
            pValUpIntronChunk[0], pValTg5pChunk[0], pValTg3pChunk[0],
            pValDnIntronChunk[0], pValDnExonIntronChunk[0],
            pValDnExon5pChunk[0]
        ])
    pValList.sort(key=operator.itemgetter(2))
    fileOpen = open(outPath + "/" + resultFile, "w")
    fileOpen.write(
        "RBP\tsmallest_p_in_upstreamExon-3prime\tsmallest_p_in_upstreamExonIntron\tsmallest_p_in_upstreamIntron\tsmallest_p_in_targetExon-5prime\tsmallest_p_in_targetExon-3prime\tsmallest_p_in_downstreamIntron\tsmallest_p_in_downstreamExonIntron\tsmallest_p_in_downstreamExon-5prime\n"
    )
    for sortedPvalue in pValList:
        fileOpen.write(sortedPvalue[0] + '\t' + str(sortedPvalue[1]) + '\t' +
                       str(sortedPvalue[2]) + '\t' + str(sortedPvalue[3]) +
                       '\t' + str(sortedPvalue[4]) + '\t' +
                       str(sortedPvalue[5]) + '\t' + str(sortedPvalue[6]) +
                       '\t' + str(sortedPvalue[7]) + '\t' +
                       str(sortedPvalue[8]) + '\n')
    fileOpen.close()
    ## End of minPvalueOut()


## end of function

########## 0. preprocessing. Getting PYGR object ####################
#
#
logging.debug("================================")
logging.debug("GETTING GENOME FASTA OBJECT")
#genomeToPygr = {'hg38': seqdb.BlastDB(pygr_station_PATH+'hg38/hg38'), 'hg19': seqdb.BlastDB(pygr_station_PATH+'hg19/hg19'), 'dm3': seqdb.BlastDB(pygr_station_PATH+'dm3/dm3'), 'dm6': seqdb.BlastDB(pygr_station_PATH+'dm6/dm6'), 'rn6': seqdb.BlastDB(pygr_station_PATH+'rn6/rn6'), 'ce11': seqdb.BlastDB(pygr_station_PATH+'ce11/ce11'), 'danRer10': seqdb.BlastDB(pygr_station_PATH+'danRer10/danRer10'), 'mm10': seqdb.BlastDB(pygr_station_PATH+'mm10/mm10'), 'araTha1': seqdb.BlastDB(pygr_station_PATH+'araTha1/araTha1'), 'danRer11': seqdb.BlastDB(pygr_station_PATH+'danRer11/danRer11'), 'oSa7': seqdb.BlastDB(pygr_station_PATH+'oSa7/oSa7'), 'xenLae2': seqdb.BlastDB(pygr_station_PATH+'xenLae2/xenLae2'), 'xenTro9': seqdb.BlastDB(pygr_station_PATH+'xenTro9/xenTro9'), 'susScr11': seqdb.BlastDB(pygr_station_PATH+'susScr11/susScr11'), 'galGal5': seqdb.BlastDB(pygr_station_PATH+'galGal5/galGal5'), 'galGal6': seqdb.BlastDB(pygr_station_PATH+'galGal6/galGal6')};
#tempDB = seqdb.SequenceFileDB(pygr_station_PATH+'mm10/mm10'); ## to make proper seqlen and pureseq file
#genomeToPygr = {'mm10': seqdb.BlastDB(pygr_station_PATH+'mm10/mm10')};
#genomeToPygr = {'hg38': seqdb.BlastDB(pygr_station_PATH+'hg38/hg38')};
#genomeToPygr = {'hg19': seqdb.BlastDB(pygr_station_PATH+'hg19')};
genome = load_genome(args.genome, pygr_station_PATH)
logging.debug("DONE GETTING PYGR OBJECT")
logging.debug("================================")

########## 1. Making input files for all three categories #################
#
#
logging.debug("================================")
logging.debug("MAKING INPUT FILES FROM rMATS FILE")
try:
    #if rMATS=="NA":
    #  logging.debug("rMATS file was not provided. Skip making up, dn, and bg files..");
    #else:
    #  pass;
    makeInputFiles(rMATS)
except:
    logging.debug("There is an exception in making input from rMATS output")
    logging.debug("Exception: %s" % sys.exc_info()[0])
    logging.debug("Detail: %s" % sys.exc_info()[1])
    sys.exit(-1)
logging.debug("DONE MAKING INPUT FILES FROM rMATS FILE")
logging.debug("================================")

#
########## 2. Making fasta files  #################
#
logging.debug("================================")
logging.debug("MAKING FASTA FILES")
try:
    getFasta('up')
    getFasta('dn')
    getFasta('bg')
    for jj in ['up', 'dn', 'bg']:
        rg = region[0]
        fFile = open(fastaPath + '/' + jj + '.' + rg + '.fasta')
        c = 0
        for dummy in fFile:
            noUse = next(fFile).strip()
            c += 1

        if jj == 'up':
            uNum = c
        elif jj == 'dn':
            dNum = c
        elif jj == 'bg':
            bNum = c
            fFile.close()

    logging.debug(
        "Number of events for upregulated, downregulated, and background: %d, %d, %d"
        % (uNum, dNum, bNum))
    pass

except:
    logging.debug("There is an exception in making fasta files")
    logging.debug("Exception: %s" % sys.exc_info()[0])
    logging.debug("Detail: %s" % sys.exc_info()[1])
    sys.exit(-2)
logging.debug("DONE MAKING FASTA FILES")
logging.debug("================================")

#sys.exit(0); ## check to see if it make fasta files okay

#
########## 3. Counting motifs for combined #################
#
#logging.debug("================================");
#logging.debug("COUNTING MOTIFS");
#try:
#  initMotif(mFile,mCount);
#  countMotif(mCount);
#  pass;
#
#except:
#  logging.debug("There is an exception in counting motifs");
#  logging.debug("Exception: %s" % sys.exc_info()[0]);
#  logging.debug("Detail: %s" % sys.exc_info()[1]);
#  sys.exit(-3);
#logging.debug("DONE COUNTING MOTIFS");
#logging.debug("================================");

#
########## 4. Plotting motifs for combined #################
#
#logging.debug("================================");
#logging.debug("PLOTTING MOTIFS");
#try:
#  #plotMotifs(mCount, 'ESRP-like');
#  plotMotifs(mCount, motifLabel);
#  pass;
#
#except:
#  logging.debug("There is an exception in counting motifs");
#  logging.debug("Exception: %s" % sys.exc_info()[0]);
#  logging.debug("Detail: %s" % sys.exc_info()[1]);
#  sys.exit(-4);
#logging.debug("DONE PLOTTING MOTIFS");
#logging.debug("================================");

#
########## 5. MAP for individual motifs  ##############
#
logging.debug("================================")
logging.debug("MAKING INDIVIDUAL MAPS")
#makeIndividualMaps(kFile);
try:
    kFile.seek(0)
    ## rewind file
    tHeader = kFile.readline().strip()
    ## header

    d = []
    for line in kFile:
        makeIndividualMaps(d, line)

    for i in range(len(d)):

        [mapname, maxPointValue, maxNegPval, drawPoints, negPvalPoints] = d[i]
        plotMotifs_finale(mapname, maxPointValue, maxNegPval, drawPoints,
                          negPvalPoints, args.separate)

    stop = timeit.default_timer()
    print(stop - start)
    if args.motif != 'NA':
        mFile.seek(0)
        ## rewind file
        tHeader = mFile.readline().strip()
        ## header

        d2 = []
        for line2 in mFile:
            makeIndividualMaps(d2, line2)

        for i in range(len(d2)):

            [mapname, maxPointValue, maxNegPval, drawPoints,
             negPvalPoints] = d2[i]
            plotMotifs_finale(mapname, maxPointValue, maxNegPval, drawPoints,
                              negPvalPoints)

    pass

except:
    logging.debug("There is an exception in making individual maps")
    logging.debug("Exception: %s" % sys.exc_info()[0])
    logging.debug("Detail: %s" % sys.exc_info()[1])
    sys.exit(-3)
logging.debug("DONE MAKING INDIVIDUAL MAPS")
logging.debug("================================")

#
########## 5. Printing out motif counts  #################
#
###logging.debug("================================");
###logging.debug("PRINTING OUT MOTIF COUNTS");
###
###upFile = open(outPath+'/'+'upregulated.motif.txt', 'w');
###dnFile = open(outPath+'/'+'downregulated.motif.txt', 'w');
###bgFile = open(outPath+'/'+'background.motif.txt', 'w');
###cbFile = open(outPath+'/'+'combined.motif.txt', 'w');
###
###try:
###  printCounts(upc, upFile, uNum); ## print counts
###  printCounts(dnc, dnFile, dNum); ## print counts
###  printCounts(bgc, bgFile, bNum); ## print counts
###  printCombinedCounts(upc,dnc,bgc,cbFile,uNum,dNum,bNum); ## print combined counts
###except:
###  logging.debug("There is an exception in printCounts function");
###  logging.debug("Exception: %s" % sys.exc_info()[0]);
###  logging.debug("Detail: %s" % sys.exc_info()[1]);
###  print "There is an exception. Please check your log file";
###  sys.exit(-5);
###logging.debug("DONE PRINTING OUT MOTIF COUNTS");
###logging.debug("================================");
###
###upFile.close(); dnFile.close(); bgFile.close(); cbFile.close();
#
########## 5. Sorting RBPs by mininum p-value  #################
#
###
logging.debug("================================")
logging.debug("SORTING BBPs by minimum P-values")
#
### upFileList
#
logging.debug("processing up regulated exons..")
minPvalueOut("up")
logging.debug("Done processing up regulated exons..")
#
### dnFileList ###
#
logging.debug("processing dn regulated exons..")
minPvalueOut("down")
logging.debug("Done processing dn regulated exons..")
#
logging.debug("DONE SORTING BBPs by minimum P-values")
logging.debug("================================")
###
#

if args.motif != "NA":
    mFile.close()

#############
## calculate total running time
#############
logging.debug("Program ended")
currentTime = time.time()
runningTime = currentTime - startTime
## in seconds
logging.debug("Program ran %.2d:%.2d:%.2d" %
              (runningTime / 3600,
               (runningTime % 3600) / 60, runningTime % 60))

sys.exit(0)


