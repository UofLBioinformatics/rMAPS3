#
### script from UofL
#
## rMATS
#
python ../motifMap.py -g hg19 -o PC3E.GS689.Li -r /media/bio/jwpark/work/motifTool/underDevelopment/testData/SE.MATS.ReadsOnTargetAndJunctionCounts.txt -mi NA -u NA -d NA -b NA --intron 250 --exon 50 --window 50 --step 1 --label "PC3E.vs.GS689.Li ESRP-like" -p /media/bio/data/pygr/ -m ../data/ESRP.like.motif.txt -k ../data/knownMotifs.human.mouse.txt 
#
### no optional motif
#
python ../motifMap.py -g hg19 -o without.opt.motif.PC3E.GS689.Li -r /media/bio/jwpark/work/motifTool/underDevelopment/testData/SE.MATS.ReadsOnTargetAndJunctionCounts.txt -mi NA -u NA -d NA -b NA --intron 250 --exon 50 --window 50 --step 1 --label "PC3E.vs.GS689.Li ESRP-like" -p /media/bio/data/pygr/ -m NA -k ../data/knownMotifs.human.mouse.txt 
#
### running with MISO
#
## hg19
#
python ../motifMap.py -g hg19 -o hg19.miso.test -r NA -mi /media/bio/jwpark/work/motifTool/underDevelopment/testData/ESRP.OE.miso_bf -u NA -d NA -b NA --intron 250 --exon 50 --window 50 --step 1 --label "hg19.miso ESRP-like" -p /media/bio/data/pygr/ -m ../data/ESRP.like.motif.txt -k ../data/knownMotifs.human.mouse.txt 
#
## dm3
python ../motifMap.py -g dm3 -o dm3.miso.test -r NA -mi /media/bio/jwpark/work/motifTool/underDevelopment/testData/example.dm3.miso_bf -u NA -d NA -b NA --intron 250 --exon 50 --window 50 --step 1 --label "dm3.miso ESRP-like" -p /media/bio/data/pygr/ -m ../data/ESRP.like.motif.txt -k ../data/knownMotifs.human.mouse.txt 
#
#
#### ignore this part for now
#
#python ../motifMap.py -g mm10 -o r.motifMap.DKO.WT -r /u/project/yxing/sstein93/Epidermis.Mouse.UPenn/MATS.DKO.WT.0001/MATS_output/SE.MATS.ReadsOnTargetAndJunctionCounts.txt -u NA -d NA -b NA --intron 250 --exon 50 --window 50 --step 1 --label "DKO.vs.WT ESRP-like" -p /u/nobackup/yxing/jwpark/data/pygr/ -m ../data/ESRP.like.motif.txt -k ../data/knownMotifs.human.mouse.txt
#python ../motifMap.py -g mm10 -o c.motifMap.DKO.WT -r NA -u upregulated.txt -d downregulated.txt -b background.txt --intron 250 --exon 50 --window 50 --step 1 --label "DKO.vs.WT ESRP-like" -p /u/nobackup/yxing/jwpark/data/pygr/ -m ../data/ESRP.like.motif.txt -k ../data/knownMotifs.human.mouse.txt
#python ../motifMap.py -g mm10 -o r.motifMap.iPS0.iPS7 -r NA -u /media/bio/jwpark/work/motifTool/underDevelopment/testData/rMAPS.highIniPS0.unique.txt -d /media/bio/jwpark/work/motifTool/underDevelopment/testData/rMAPS.highIniPS7.unique.txt -b /media/bio/jwpark/work/motifTool/underDevelopment/testData/rMAPS.background.unique.txt --intron 250 --exon 50 --window 50 --step 1 --label "DKO.vs.WT ESRP-like" -p /media/bio/data/pygr/ -m ../data/ESRP.like.motif.txt -k ../data/knownMotifs.human.mouse.txt
