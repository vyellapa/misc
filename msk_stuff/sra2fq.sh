SRA=$1

fastq-dump -I --split-files ${SRA}

if [ "$?" == "0" ]
then

echo "Delete ${SRA}"

fi

