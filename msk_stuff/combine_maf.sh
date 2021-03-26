rm -f all.maf
cat *.maf | egrep "^Hugo" |head -1 > temp
cat *.maf | egrep -v "^#|^Hugo" |sort -u  >> temp

mv temp all.maf
