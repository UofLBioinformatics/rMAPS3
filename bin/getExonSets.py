#
## this program gets the rMATS + cuffdiff output to make upregulated, downregulated, and background SE events
#

### import necessary libraries
import re,os,sys,logging,time,datetime;

### checking out the number of arguments
if (len(sys.argv)<5): 
  print('Not enough arguments!!');
  print ('It takes at least 4 arguments.');
  print ('Usage:\n\tProgramName rMATS.cuffdiff sample_1_name sample_2_name outFolder');
  print ('Example\n\tProgramName /u/project/yxing/sstein93/Epidermis.Mouse.UPenn/MATS.DKO.WT.0001/MATS_output/MATS.cuffdiff/common.txt DKO WT DKO.WT');
  sys.exit();

def listToString(x):
  rVal = '';
  for a in x:
    rVal += a+' ';
  return rVal;


### setting up the logging format 
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(message)s',
                    filename='log.getExonSets.'+  str(datetime.datetime.now())+'.txt',
                    filemode='w')

##### Getting Start Time ######
logging.debug('Start the program with [%s]\n', listToString(sys.argv));
startTime = time.time();

###
iFile = open(sys.argv[1]); ## input file
S1 = sys.argv[2]; ## sample 1 name
S2 = sys.argv[3]; ## sample 2 name
outDir = sys.argv[4]; ## output folder
os.system('mkdir -p ' + outDir);
#
outPath = os.path.abspath(outDir); ## absolute output path
#
upFile=open(outPath+'/highIn'+S1+'.txt','w');
dnFile=open(outPath+'/highIn'+S2+'.txt','w');
bgFile=open(outPath+'/background.txt','w');
#
################### global variables #######################
#
#nUP=0; nDN=0; nBG=0; 
sigFDR=0.05;  ## fdr must be smaller than this to be significant
sigDeltaPSI=0.05; ## abs(deltaPSI) must be greater than or equal to this to be significant
bgFDR=0.5;  ## fdr must be greater than this to be background
bgFPKM=5.0; ## FPKM must be greater than this to be background
minPSI=0.85; maxPSI=0.15;
u={}; d={}; b={};
#
#
############################################################
#
hh=iFile.readline().strip().split('\t');
header = '\t'.join(hh[1:11]);
upFile.write(header+'\n');
dnFile.write(header+'\n');
bgFile.write(header+'\n');

for line in iFile: ## for each line
  ele = line.strip().split('\t');
  key = ':'.join(ele[3:7]);###print key; sys.exit();
  value = [1, '\t'.join(ele[1:11])];
  fdr = float(ele[19]);
  if fdr<sigFDR: ## it could be significant 
    deltaPSI=float(ele[22]);
    if deltaPSI>=sigDeltaPSI: ## it's upregulated. high in sample 1
      u[key]=value;
    elif deltaPSI<=-sigDeltaPSI: ## it's downregulated. high in sample 2
      d[key]=value;
  elif fdr>bgFDR and min(float(ele[-6]),float(ele[-5]))>bgFPKM: ## it could be background
    psi = ele[20].replace('"','').split(',');
    t=0;sum=0.0;
    PSI_1=0.0;
    for p in psi: ##
      if p!="NA": ## number here
        t+=1;
        sum+=float(p);
    if t>0:
      PSI_1 = sum/t;
  
    psi = ele[21].replace('"','').split(',');
    t=0;sum=0.0;
    PSI_2=0.0;
    for p in psi: ##
      if p!="NA": ## number here
        t+=1;
        sum+=float(p);
    if t>0:    
      PSI_2 = sum/t;

    if min(PSI_1,PSI_2)<minPSI and max(PSI_1,PSI_2)>maxPSI: ## it is a background event
      b[key]=value;

logging.debug("Done populating initial dictionaries with possible duplicates");
logging.debug("Number of up, down, and background exons are: %d, %d, %d" % (len(u),len(d),len(b)));

logging.debug("Removing exons included in more than one dictionaries..");

for key in u: ## going through u
  if key in d: ## same exon in d
    d[key][0]+=1;
    u[key][0]+=1;
  if key in b: ## same exon in b
    b[key][0]+=1;
    u[key][0]+=1;

for key in d: ## going through d
  if key in b: ## same exon in b
    b[key][0]+=1;
    d[key][0]+=1;

### now write to up,down,background only if count value is 1
nu=0;nd=0;nb=0; ## number of up,down,background
for key in u:
  if u[key][0]==1: ## it is unique
    nu+=1;
    upFile.write(u[key][1] + '\n');
for key in d:
  if d[key][0]==1: ## it is unique
    nd+=1;
    dnFile.write(d[key][1] + '\n');
for key in b:
  if b[key][0]==1: ## it is unique
    nb+=1;
    bgFile.write(b[key][1] + '\n');

logging.debug("Number of upregulated (high in sample_1) exons: %d" % nu);
logging.debug("Number of downregulated (high in sample_2) exons: %d" % nd);
logging.debug("Number of background exons: %d" % nb);

iFile.close();
upFile.close();
dnFile.close();
bgFile.close();

#############
## calculate total running time
#############
logging.debug("Program ended");
currentTime = time.time();
runningTime = currentTime-startTime; ## in seconds
logging.debug("Program ran %.2d:%.2d:%.2d" % (runningTime/3600, (runningTime%3600)/60, runningTime%60));

sys.exit(0);
