#!/bin/sh


pathlist=$1
LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p $pathlist)

echo "SCRIPT START FOR $LINE ----------------------------"
printf "Path to genomeark file:  $LINE \n\n"

IFS='/' read -r -a array <<< $LINE 
fastaname=${array[4]} ##make more robust - basename 
echo "Local filename: $LINE"

printf "Path to genomeark file:  $LINE \n\n"

statsname=$LINE.gfastats

aws s3 cp s3://genomeark/${LINE} ./${fastaname}

ziptime=$(wc -c $fastaname)
echo "Time on compressed $fastaname"
TIMEFORMAT=%R
ziptime=$(time (gfastats $fastaname > ${fastaname}_temp_out.txt) 2>&1)

ziplength=$(grep "Total scaffold length" ${fastaname}_temp_out.txt | grep -Eo "[0-9]+") 
echo "Compressed length: $ziplength"
printf "$LINE\t $ziplength\t $ziptime\t gzip\n"  >> gfastats_stats_out.txt

uncomp=$(echo $fastaname | sed 's/.gz//g') 
echo "Decompressing $fastaname \n\n"
gunzip $fastaname 

echo "Time on uncompressed fasta: $uncomp"
unziptime=$(time (gfastats $uncomp > ${fastaname}_temp_out_2.txt) 2>&1)

unziplength=$(grep "Total scaffold length" ${fastaname}_temp_out_2.txt | grep -Eo "[0-9]+")
echo "Decompressed length: $unziplength"

printf "$LINE\t $unziplength\t $unziptime\t plain text\n" >> gfastats_stats_out.txt

rm $uncomp
rm ${fastaname}_temp_out.txt 
rm ${fastaname}_temp_out_2.txt 

echo "SCRIPT COMPLETE FOR $LINE --------------------"

