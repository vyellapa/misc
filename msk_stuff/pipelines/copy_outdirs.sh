PJ=$1
lk get_outdirs -fi name CAVEMAN -fi projects__pk ${PJ} |awk '{print "rsync -aP "$1"/*.caveman.muts.annot.vcf.gz ./"}'
