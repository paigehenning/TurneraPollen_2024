#!/bin/bash
echo "Please provide name for output files. i.e. name of species or project"
read varname
output_DownloadSRA=$varname"downloadSRA.sh"
output_Map=$varname"map.sh"
output_quantification=$varname"quantificationRNAseq.sh"
echo ""
echo ""
echo "Please provide the name for your genome including file extension (i.e. fasta)"
read genome

while true; do
    read -p "Your genome file is named $genome, is this correct?" yn
    case $yn in
       [Yy]* ) echo "Great!";break;;
        [Nn]* ) echo "restarting";bash process_SRA.sh;exit;;
        * ) echo "Please answer yes or no.";;
    esac
done

echo "Please provide the name for your gff including file extension (i.e. fasta)"
read gff

while true; do
    read -p "Your gff file is named $gff, is this correct?" yn
    case $yn in
       [Yy]* ) echo "Great!";break;;
        [Nn]* ) echo "restarting";bash process_SRA.sh;exit;;
        * ) echo "Please answer yes or no.";;
    esac
done
cp $gff gff.gff

echo ""
echo "Your bash script for downloading SRA files from genbank will be named $output_DownloadSTA"
echo ""
echo "Your bash scrip for mapping will be named $output_Map"
echo ""
echo "Your bash script for quantification of RNAseq will be named $output_quantification"

cat SraRunTable.txt| awk -F '\t' '
										{print "fasterq-dump",$1}}' > temp

mkdir map
awk -F '\t' '{print "hisat2 -p 12 --dta -x index/index -1",$1"_1.fastq -2",$1"_2.fastq -S sam/"$1"\nwait\nsamtools sort -@ 8 -o map/"$1".bam sam/"$1,"\nwait\nrm sam/"$1}' SraRunTable.txt >  temp2
mkdir RNAseq

echo "#!/usr/bin/"> $output_DownloadSRA
echo "echo ' Now  downloading rnaseq files'" >> $output_DownloadSRA
paste temp >> $output_DownloadSRA
echo "bash $output_Map">> $output_DownloadSRA


echo "#!/usr/bin/">$output_Map
echo "mkdir sam">$output_Map
echo "echo 'Now mapping RNA'" >> $output_Map
paste temp2 >>$output_Map
echo "rm -rf sam/" >> $output_Map
echo "bash $output_quantification">>$output_Map
rm temp*

mkdir expression
mkdir ballgown

awk -F '\t' '{if ($1=="") {} else if ($1=="Run"){} else {print "stringtie map/"$1".bam -G gff.gff -p 20 -e -A expression/"$1"_gene_abund.tab -b ballgown/rep"$1"_gen_abundance"}}' SraRunTable.txt > temp

echo "#!/usr/bin/"> $output_quantification
echo "echo 'Now quantifying gene expression'" >> $output_quantification
paste temp >> $output_quantification
echo  "bash cleanup.sh">> $output_quantification
rm temp*

echo "#!/usr/bin/" >cleanup.sh
echo "rm temp*" >>cleanup.sh
echo "rm gff.gff">>cleanup.sh
echo "mkdir shellscripts">>cleanup.sh
echo "mv *fastq RNAseq" >> cleanup.sh
echo "mv *.sh shellscripts/">> cleanup.sh
echo "ls expression/ > temp">> cleanup.sh
echo "cat temp| awk '{print "sort -o expression/"$1,"expression/"$1}'">> cleanup.sh
echo "awk -i inplace -v ORS='\r\n' 'FNR==1{print FILENAME}1' expression/*">> cleanup.sh
echo "echo 'analysis complete!'" >> cleanup.sh

#bash $output_DownloadSRA
