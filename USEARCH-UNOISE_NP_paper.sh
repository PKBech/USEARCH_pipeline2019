#wget http://www.drive5.com/downloads/usearch11.0.667_i86linux32.gz
#collect each F and R files from qiime artifact
qiime tools extract \
    --input-path demultiplexed-seqs_NP_paper.qza \
  --output-path demultiplexed-seq_NP_paper-files
#move data to Input dir in Usearch-unoise dir

#Removine primers from F and R
for file in Input/*.fastq
  do
  # Create variable filename with path of current file
  filename="$file"
  # If current path / filename contains "_1" ...
  if [[ $filename = *"R1"* ]]
  then
  # Generate the new file name with appended TP to indicate trimmed primer file
  newfilename="$(echo "$filename" | sed s/R1_001.fastq/R1_TP.fq/)"
  # Write to console what file is being processed...
  echo "Processing $newfilename"
  # Trim primer in file using usearch...
  ./usearch -fastx_truncate "$filename" -stripleft 18 -fastqout "$newfilename"
  # Same process as above for files containing *R2*
  elif [[ $filename = *"R2"* ]]
  then
  newfilename="$(echo "$filename" | sed s/R2_001.fastq/R2_TP.fq/)"
  echo "Processing $newfilename"
  ./usearch -fastx_truncate "$filename" -stripleft 20 -fastqout "$newfilename"
  mkdir ../Primer_trimmed-files
  mv `ls "$file" | grep "TP.fq"` ../Primer_trimmed-files
  fi
  done

mkdir Primer_trimmed-files
mv Input/*TP.fq Primer_trimmed-files

#Trancates F and R 

for file in Primer_trimmed-files/*
  do
  # Create variable filename with path of current file
  filename="$file"
  
  # If current path / filename contains "_1" ...
  if [[ $filename = *"R1"* ]]
  then
  
  # Generate the new file name with appended T to indicate truncated file
  newfilename="$(echo "$filename" | sed s/R1_TP.fq/R1.fq/)"
  
  # Write to console what file is being processed...
  echo "Processing $newfilename"
  
  # Truncate file using usearch...
  ./usearch -fastx_truncate "$filename" -trunclen 240 -fastqout "$newfilename"
  
  # Same process as above for files containing _2
  elif [[ $filename = *"R2"* ]]
  then
  newfilename="$(echo "$filename" | sed s/R2_TP.fq/R2.fq/)"
  echo "Processing $newfilename"
  ./usearch -fastx_truncate "$filename" -trunclen 220 -fastqout "$newfilename"
  fi
 done

mkdir Truncated-files_new
mv Primer_trimmed-files/*T.fq Truncated-files_new


#All forward and reverse reads are sorted according to their respective pairs with the script resync.pl downloaded from gibhub

for file in Truncated-files_new/*
do	
# Create variable filename with path of current file
	filename="$file"
	if [[ $filename = *"R1"* ]]
	then
	R1="$(echo "$filename" | sed s/R1.fq//)"
	elif [[ $filename = *"R2"* ]]
        then
        R2="$(echo "$filename" | sed s/R2.fq//)"
	#echo "$R1"
	fi
	if [ "$R1" == "$R2" ]
	then
	#newfilename="$(echo "$file" | sed s/.fq/_Sort.fq)"
	#echo "Processing $newfilename"
	perl resync.pl --fwd="$R1"R1.fq -rev="$R2"R2.fq 
	fi
done

mkdir sync_files
mv Truncated-files_new/*.syns.fastq sync_files

#Merge files
for file in sync_files/*
	newfilename="$(echo "$file" | sed s/_R1.sync.fastq/_merged.fq/)"
	echo "Processing $newfilename"
	./usearch -fastq_mergepairs *R1.sync.fastq -sample @ -fastqout "$newfilename"
	done

for file in sync_files/*R1*
do	
	newfilename="$(echo "$file"| sed s/1.sync.fastq/_merged.fq/)"
	echo "Processing $newfilename"
	./usearch -fastq_mergepairs "$file" -sample @ -fastqout "$newfilename" 
done

mkdir merged-files
mv sync_files/*merged.fq merged-files

#Quality Filtering
for file in merged-files/*
do	
	newfilename="$(echo "$file"| sed s/_merged.fq/_Qfiltered.fa/)"
	echo "Processing $newfilename"
	./usearch -fastq_filter "$file" -fastq_maxee 1.0 -fastaout "$newfilename" 
done

mkdir Qfiltered-seqs
mv merged-files/*_Qfiltered.fa Qfiltered-seqs

#Deduplicating
for file in Qfiltered-seqs/*
do	
	newfilename="$(echo "$file"| sed s/_Qfiltered.fa/_Uniques.fa/)"
	echo "Processing $newfilename"
	./usearch -fastx_uniques "$file" -fastaout "$newfilename" -sizeout -relabel Uniq
done

mkdir uniq-seqs
mv Qfiltered-seqs/*_Uniques.fa uniq-seqs

#Pool all uniques to give the highest possible abundances for correct reads.

 

#UNOISE3 denoised sequences (ZOTUs)
for file in uniq-seqs/*
do	
	newfilename="$(echo "$file"| sed s/_Uniques.fa/_unoise3/)"
	echo "Processing $newfilename"
	./usearch -unoise3 "$file" -zotus "$newfilename"_zotus.fa -tabbedout "$newfilename".txt
done

mkdir denoised-seqs
mv uniq-seqs/*_unoise3* denoised-seqs

#Pool all denoised ZOTUs 
cat ./6_AD_denovo_nonchimera/*_denovo.nonchimeras.fa > ./6_AD_denovo_nonchimera/all.denovo.nonchimeras.fa

#MAke ASV tables

for file in ./6_AD_denovo_nonchimera/*_denovo.nonchimeras.fa
do	
	newfilename="$(echo "$file" | sed s/_denovo.nonchimeras.fa/_OTUtab.txt/)"
	echo "Processing $newfilename"
	./usearch -otutab  "$file" -otus db_all_AD_95_filtered.fa -otutabout "$newfilename" 
done

for file in ./6_AD_denovo_nonchimera/*_denovo.nonchimeras.fa
do	
	newfilename="$(echo "$file" | sed s/_denovo.nonchimeras.fa/_OTUtab.txt/)"
	echo "Processing $newfilename"
	./usearch -otutab  "$file" -otus db_all_AD_95_filtered.fa -otutabout "$newfilename" 
done


usearch -otutab reads.fq -zotus zotus.fa -otutabout zotutab.txt -mapout zmap.txt