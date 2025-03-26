#!/bin/bash

# Define the reference genome ID and paths
ref_id="GCF_000146045.2"
ref_fna="/share/home/helab/weixianfang/MTDtool/ref/GCF_000146045.2_R64_genomic.fna"
species_name="Saccharomyces_cerevisiae"
indelfile_path="/data/home/helab/weixianfang/trfinder-TRfinder/rsgb/rsgb_assembly_summary/fungi_ni/sc-zhongpao/mugsy/Saccharomyces_cerevisiae"

# Step 1: Filter the fungal genomes
while IFS= read -r a; do
    raw_fna=$(find /data/home/helab/weixianfang/trfinder-TRfinder/rsgb/fungi/ -name "$a*genomic.fna.gz")
    if [[ -n "$raw_fna" && "$a" != "$ref_id" ]]; then
        echo -e "$a" >> Saccharomyces-cerevisiae-filter.txt
    fi
done < /data/home/helab/weixianfang/trfinder-TRfinder/rsgb/rsgb_assembly_summary/fungi_ni/sc-zhongpao/mugsy/Saccharomyces-cerevisiae.txt

# Step 2: Run Mugsy for whole genome alignment and insertions detection
#SBATCH --job-name=sc3
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
    python /data/home/helab/weixianfang/trfinder-TRfinder/mugsy-multiple-whole-genome-alignment/python_concordance/identy_insert-fwd-rev-done.py "$a_name".maf "$ref_id" "$a"_insert.tsv

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
