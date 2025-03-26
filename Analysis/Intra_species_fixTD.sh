#!/bin/bash

# Define the reference genome ID and paths
ref_id="" ##Please enter the file address here
ref_fna="" ##Please enter the file address here
species_name="" ##Please enter the file address here
indelfile_path="" ##Please enter the file address here

# Step 1: Filter the fungal genomes
while IFS= read -r a; do
    raw_fna=$(find "" -name "$a*genomic.fna.gz") ##Please enter the file address here
    if [[ -n "$raw_fna" && "$a" != "$ref_id" ]]; then
        echo -e "$a" >> Saccharomyces-cerevisiae-filter.txt
    fi
done < ##Please enter the file address here

# Step 2: Run Mugsy for whole genome alignment and insertions detection
#SBATCH --job-name=""  ##Please enter here
#SBATCH --partition=cpu
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --time=336:00:00
#SBATCH --mem-per-cpu=70gb
#SBATCH --cpus-per-task=1

# Directory for storing files
cd $indelfile_path

# Counter for the processed species
inter_count=0

# Process each genome
while IFS= read -r a; do
    raw_fna=$(find /data/home/helab/weixianfang/trfinder-TRfinder/rsgb/fungi/ -name "$a*genomic.fna.gz")
    cp "$raw_fna" "$indelfile_path""$a".fna.gz
    gzip -d -c "$indelfile_path""$a".fna.gz > "$indelfile_path""$a".fna
    ((inter_count++))

    # Generate MUGSY alignment and identify insertions
    a_name=$(echo "$a" | sed 's/\./_/g')
    mugsy --directory "$indelfile_path" --prefix "$a_name" "$ref_fna" "$indelfile_path""$a".fna
    python identy_insert-fwd-rev-done.py "$a_name".maf "$ref_id" "$a"_insert.tsv

    # Append insertions to the cumulative file
    cat "$indelfile_path""$a"_insert.tsv >> "${indelfile_path}${species_name}_insert.tsv"

    # Clean up temporary files
    rm "$indelfile_path""$a".fna "$indelfile_path""$a".fna.gz "$a_name".maf "$a_name" "$a_name".mugsy.log
done < /data/home/helab/weixianfang/trfinder-TRfinder/rsgb/rsgb_assembly_summary/fungi_ni/sc-zhongpao/mugsy/Saccharomyces-cerevisiae-filter.txt

# Remove duplicate entries
awk '!seen[$1" "$2" "$4]++' "${indelfile_path}${species_name}_insert.tsv" > "${indelfile_path}${species_name}_insert_rmdup.tsv"

# Get the number of insertions
insert_num=$(wc -l < "${indelfile_path}${species_name}_insert_rmdup.tsv")

# Step 3: Process insertion sequences
awk -F'\t' 'BEGIN {OFS="\t"} {print $1, $2, $2+$5}' "$indelfile_path${species_name}_insert_rmdup.tsv" > ecoli-indel.bed
while IFS=$'\t' read -r chr start end; do
    echo -e "$chr\t$start\t$end" > position_per.bed
    bedtools getfasta -fi "$ref_fna" -bed position_per.bed -fo ecoli-indel-seq.fa
    rm position_per.bed

    if [[ -s ecoli-indel-seq.fa ]]; then
        awk '/^>/{getline; print}' ecoli-indel-seq.fa >> ecoli-indel-seq-2.txt
    else
        echo " " >> ecoli-indel-seq-2.txt
    fi
done < ecoli-indel.bed

# Combine insertion data with sequences
paste -d $'\t' "$indelfile_path${species_name}_insert_rmdup.tsv" ecoli-indel-seq-2.txt > ecoli-indel-4bp-7.txt
rm ecoli-indel-seq-2.txt

# More processing
awk -F'\t' 'BEGIN {OFS="\t"} {print $1, $2+$5, $2+$5+$5}' "$indelfile_path${species_name}_insert_rmdup.tsv" > ecoli-indel.bed
while IFS=$'\t' read -r chr start end; do
    echo -e "$chr\t$start\t$end" > position_per.bed
    bedtools getfasta -fi "$ref_fna" -bed position_per.bed -fo ecoli-indel-seq.fa
    rm position_per.bed

    if [[ -s ecoli-indel-seq.fa ]]; then
        awk '/^>/{getline; print}' ecoli-indel-seq.fa >> ecoli-indel-seq-2.txt
    else
        echo " " >> ecoli-indel-seq-2.txt
    fi
done < ecoli-indel.bed

# Final combination of sequences
paste -d $'\t' ecoli-indel-4bp-7.txt ecoli-indel-seq-2.txt > ecoli-indel-4bp-8.txt
rm ecoli-indel-seq-2.txt

# Filter sequences with no 'N' characters
awk -F'\t' 'BEGIN {OFS="\t"} $4 !~ /N/ {print}' ecoli-indel-4bp-8.txt > tmp-8-2.txt

# Output final insertion data
awk -F'\t' 'BEGIN {OFS="\t"} $4 == $6 && $6 != $7 {print}' tmp-8-2.txt > tmp1_final.txt
awk '{if ($6 != $7) {for (i=1; i<=length($4) && i<=length($7); i++) {if(substr($4, i, 1) != substr($7, i, 1)) {print $0,"\t",i-1; break}}}}' tmp1_final.txt > "${indelfile_path}${species_name}_MTD.tsv"

# Cleanup temporary files
rm tmp-8-1.txt tmp-8-2.txt tmp1_final.txt

# Step 4: TRF Filtering
awk -F'\t' 'BEGIN {OFS="\t"} $5 >= 8 {print}' "${indelfile_path}${species_name}_MTD.tsv" > "${indelfile_path}before_trf_file.tsv"

# TRF filtering process
cd /data/home/helab/weixianfang/trfinder-TRfinder/add_trf_analy_mha

while IFS=$'\t' read -r chr insert_start insert_alt insert_seq insert_len ref_unit1_seq ref_unit2_seq mha_lengh; do
    seqname=">corseq"
    echo -e "$seqname\n$insert_seq$ref_unit1_seq$ref_unit2_seq" > "$tmp_filetrf"corseq.fa
    trf corseq.fa 2 7 7 80 10 30 1000 -f -d -m
    python /data/home/helab/weixianfang/trfinder-TRfinder/add_trf_analy_mha/convert_script2.py
    tr_file="corseq.fa.2.7.7.80.10.30.1000-tr.dat"
    awk -v value1="$insert_len" -v value2="$insert_seq" '$4 == value1 && $15 == value2 && $4 >= 8 { print $0 > "cor_filter_trf.txt" }' "$tr_file"
    cat cor_filter_trf.txt >> "$trf_format_output"
    if [[ -s cor_filter_trf.txt ]]; then
        echo -e "$chr\t$insert_start\t$insert_alt\t$insert_seq\t$insert_len\t$ref_unit1_seq\t$ref_unit2_seq\t$mha_lengh" >> "$output_filetrf"
    fi
    rm ./corseq*
    > cor_filter_trf.txt
done < "${indelfile_path}before_trf_file.tsv"

rm "${indelfile_path}before_trf_file.tsv" cor_filter_trf.txt


# Step 5: Output statistical results
awk -F'\t' -v species="$species_name" -v inter_count="$inter_count" '
BEGIN {
    OFS="\t"; 
    count1=0; 
    count2=0; 
    count3=0; 
    count4=0; 
    count5=0
}
# Count different conditions based on the comparison of columns 6 and 7
$6 != $7 {count1++}
$6 != $7 && substr($6, 1, 1) != substr($7, 1, 1) {count2++}
$6 != $7 && substr($6, 1, 1) == substr($7, 1, 1) {count3++}
$6 != $7 && substr($6, 1, 2) == substr($7, 1, 2) {count4++}
$6 != $7 && substr($6, 1, 3) == substr($7, 1, 3) {count5++}
END {
    # Output the results into a specific file
    print species, inter_count, count1, count2, count3, count4, count5 >> "'"${mtdfile}"'fungi-insert-mtd-8bp-number_trf.tsv"
}' "${mtdfile}${species_name}_MTD_trf.tsv"


# Step 6: Count occurrences of 2x and 3x TD loci in isolates and add to final output file
output_file2x_final="${mtdfile}${species_name}_MTD_trf_final.tsv"

while IFS=$'\t' read -r chr insert_start insert_alt insert_seq insert_len ref_unit1_seq ref_unit2_seq mha_lengh; do
    # Count occurrences of 2x TD loci
    num2x=$(awk -v chr="$chr" -v pos="$insert_start" -v mut_alt="$insert_seq" -v length_bp="$insert_len" '
    BEGIN {count=0}
    $1==chr && $2==pos && $4==mut_alt {count++}
    END {print count}
    ' "${indelfile}${species_name}_insert.tsv")

    # Count occurrences of 3x TD loci
    num3x=$(awk -v chr="$chr" -v pos="$insert_start" -v mut_alt="$insert_seq" -v length_bp="$insert_len" '
    BEGIN {count=0}
    $1==chr && $2==pos && $5==length_bp*2 && substr($4, 1, length_bp)==mut_alt && substr($4, length_bp+1, length_bp)==mut_alt {count++}
    END {print count}
    ' "${indelfile}${species_name}_insert.tsv")

    # Append the result to the final output file
    echo -e "$chr\t$insert_start\t$insert_alt\t$insert_seq\t$insert_len\t$ref_unit1_seq\t$ref_unit2_seq\t$mha_lengh\t$num2x\t$num3x" >> "$output_file2x_final"
done < "${mtdfile}${species_name}_MTD_trf.tsv"

# Remove the intermediate file after processing
rm "${mtdfile}${species_name}_MTD_trf.tsv"

# The final output file "${mtdfile}${species_name}_MTD_trf_final.tsv" contains columns:
# 1. Chromosome
# 2. Insertion start position
# 3. Insertion alternate sequence
# 4. Insertion sequence
# 5. Insertion length
# 6. Reference sequence unit 1
# 7. Reference sequence unit 2
# 8. MHA length
# 9. Count of 2x occurrences in isolates
# 10. Count of 3x occurrences in isolates
