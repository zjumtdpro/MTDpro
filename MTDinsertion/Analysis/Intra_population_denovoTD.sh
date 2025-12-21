#!/bin/bash

# Initialize count variable
inter_count=0

# Add TRF to the PATH (ensure this path is correctly specified)
export PATH=/data/home/helab/weixianfang/bio-tools/trfinder/TRF-4.09.1/build/src:$PATH

# Input species names (Please enter species names here)
strname=#Please enter species names here (e.g., "species1 species2")

# Reference genome file path (Please enter the reference file path here)
ref=#Please enter the reference genome file path here (e.g., /path/to/ref.fna)

# Loop over species names
for a1 in $strname;
do
  species_name="$a1"
  
  # Input start file path (Please enter the start file path here)
  start_file=#Please enter the start file path here (e.g., /path/to/start_file.txt)

  # Output directory for MTD-related files (Please enter the output file directory here)
  mtdfile=#Please enter the output file directory here (e.g., /path/to/output/)

  # Define input and output file paths
  indelfile="$mtdfile"
  indelfile_path=#Please enter the directory path where data is stored here (e.g., /path/to/data/)

  # Filter and process insertion data from start file
  awk -F'\t' 'BEGIN {OFS="\t"; count=0} $19!="VarAllele" && length($19)>4 {print $1, $2, $3, $19}' "$start_file" \
    | awk '!arr[$0]++' \
    | awk -F'\t' -v OFS='\t' '{ if ($4 ~ /^\+/) sub(/^\+/, "", $4); else next } { print $1, $2, $3, $4 }' \
    > "$mtdfile"indel-4bp-5.txt

  # Add length of insertion to file
  awk -F'\t' 'BEGIN {OFS="\t"} {print $0, length($4)}' indel-4bp-5.txt > "$mtdfile""$a1".indel-4bp.txt

  # Calculate the number of insertions
  mtd_start_file="$mtdfile""$a1".indel-4bp.txt
  awk 'END{print "This species has " NR " insertions greater than 4bp."}'  "$mtd_start_file"

  # Clean up intermediate files
  rm "$mtdfile"indel-4bp-5.txt

  # Prepare files for further processing (converting insertions into bed format)
  awk -F'\t' 'BEGIN {OFS="\t"} {print $1, $2, $2+$5}' "$mtd_start_file" > indel.bed

  # Extract sequences from the specified positions in the reference genome using bedtools
  while IFS=$'\t' read -r chr start end; do
    echo -e "$chr\t$start\t$end" > position_per.bed
    bedtools getfasta -fi "$ref" -bed position_per.bed -fo indel-seq.fa
    rm position_per.bed

    if [[ -s indel-seq.fa ]]; then
      awk '/^>/{getline; print}' indel-seq.fa >> indel-seq-2.txt
    else
      echo " " >> indel-seq-2.txt
    fi
  done < indel.bed

  # Combine sequences with insertion data
  paste -d $'\t' "$mtd_start_file" indel-seq-2.txt > indel-4bp-7.txt
  rm indel-seq-2.txt

  # Further processing of second segment
  awk -F'\t' 'BEGIN {OFS="\t"} {print $1, $2+$5, $2+$5+$5}' "$mtd_start_file" > indel.bed
  while IFS=$'\t' read -r chr start end; do
    echo -e "$chr\t$start\t$end" > position_per.bed
    bedtools getfasta -fi "$ref" -bed position_per.bed -fo indel-seq.fa
    rm position_per.bed

    if [[ -s indel-seq.fa ]]; then
      awk '/^>/{getline; print}' indel-seq.fa >> indel-seq-2.txt
    else
      echo " " >> indel-seq-2.txt
    fi
  done < indel.bed

  # Combine the second segment with the first
  paste -d $'\t' indel-4bp-7.txt indel-seq-2.txt > indel-4bp-8.txt
  rm indel-seq-2.txt

  # Make all sequences uppercase and filter out insertions with 'N'
  awk -F'\t' 'BEGIN {OFS="\t"} {print $1, $2, $3, toupper($4), $5, toupper($6), toupper($7)}' indel-4bp-8.txt > tmp-8-1.txt
  awk -F'\t' 'BEGIN {OFS="\t"} $4 !~ /N/ {print}' tmp-8-1.txt > tmp-8-2.txt

  # Filter insertions that match between two columns
  awk -F'\t' 'BEGIN {OFS="\t"} $4 == $6 && $6 != $7 {print}' tmp-8-2.txt > tmp1_final.txt
  awk '{if ($6 != $7) {for (i=1; i<=length($4) && i<=length($7); i++) {if(substr($4, i, 1) != substr($7, i, 1)) {print $0,"\t",i-1; break}}}}' tmp1_final.txt > "${mtdfile}${species_name}_MTD.tsv"

  # Clean up intermediate files
  rm tmp-8-1.txt tmp-8-2.txt tmp1_final.txt

  # Further processing for filtering insertions based on length
  rm indel-4bp-7.txt indel-4bp-8.txt indel.bed indel-seq.fa

  # Step for TRF filtering: Filter insertions of at least 8bp
  awk -F'\t' 'BEGIN {OFS="\t"} $5 >= 8 {print}' "${mtdfile}${species_name}_MTD.tsv" > "${mtdfile}before_trf_file.tsv"

  # TRF filtering using Tandem Repeat Finder
  tmp_filetrf="/data/home/helab/weixianfang/trfinder-TRfinder/add_trf_analy_mha/"
  output_filetrf="${mtdfile}${species_name}_MTD_trf.tsv"
  trf_format_output="${mtdfile}${species_name}_trf_filter_trf_format.txt"

  # TRF filtering loop
  while IFS=$'\t' read -r chr insert_start insert_alt insert_seq insert_len ref_unit1_seq ref_unit2_seq mha_lengh; do
    seqname=">corseq"
    echo -e "$seqname\n$insert_seq$ref_unit1_seq$ref_unit2_seq" > "$tmp_filetrf"corseq.fa
    trf "$tmp_filetrf"corseq.fa 2 7 7 80 10 30 1000 -f -d -m
    python TRFconvert_script.py
    tr_file="corseq.fa.2.7.7.80.10.30.1000-tr.dat"

    # Apply TRF filter criteria
    awk -v value1="$insert_len" -v value2="$insert_seq" '$4 == value1 && $15 == value2 && $4 >= 8 { print $0 > "cor_filter_trf.txt" }' "$tr_file"
    cat cor_filter_trf.txt >> "$trf_format_output"
    
    if [[ -s cor_filter_trf.txt ]]; then
      echo -e "$chr\t$insert_start\t$insert_alt\t$insert_seq\t$insert_len\t$ref_unit1_seq\t$ref_unit2_seq\t$mha_lengh" >> "$output_filetrf"
    fi

    rm ./corseq*
    > cor_filter_trf.txt
  done < "${mtdfile}before_trf_file.tsv"

  # Clean up intermediate files
  rm "${mtdfile}before_trf_file.tsv" cor_filter_trf.txt

  # Step 7: Output statistical results
  cd "$indelfile_path"
  awk -F'\t' -v species="$species_name" -v inter_count="$inter_count" 'BEGIN {OFS="\t"; count1=0; count2=0; count3=0; count4=0; count5=0}
    $6 != $7 {count1++}
    $6 != $7 && substr($6, 1, 1) != substr($7, 1, 1) {count2++}
    $6 != $7 && substr($6, 1, 1) == substr($7, 1, 1) {count3++}
    $6 != $7 && substr($6, 1, 2) == substr($7, 1, 2) {count4++}
    $6 != $7 && substr($6, 1, 3) == substr($7, 1, 3) {count5++}
    END {
      print species, inter_count, count1, count2, count3, count4, count5 >> "'"${mtdfile}"'fungi-insert-mtd-8bp-number_trf.tsv"
    }' "${mtdfile}${species_name}_MTD_trf.tsv"
done
