# MTDfinder-pro
MTDfinder-Pro provides a unified pipeline to detect Microhomology-Mediated Tandem Duplications (MTDs) at three evolutionary levels: intra-population, intra-species, and species-level. 
# Introduction
MTDfinder-Pro is a computational framework designed to identify Microhomology-Mediated Tandem Duplications (MTDs) using Illumina short reads, PacBio HiFi reads, multi-genome alignments, and TRF-based reference genome analysis.

The pipeline performs:
Intra-population (read-based) detection of de novo MTDs
Intra-species (comparative genomic) detection of fixed/polymorphic MTDs
Species-level (reference genome) scanning of tandem duplications

Algorithmic details correspond to the revised manuscript and workflow.
<img width="865" height="479" alt="image" src="https://github.com/user-attachments/assets/c3c48904-eaf8-44bc-86a6-3772a5d90327" />

## Requirements
| Tool                    | Version |
| ----------------------- | ------- |
| python3                 | ≥3.8    |
| bwa                     | ≥0.7    |
| pbmm2 (for PacBio HiFi) | latest  |
| bedtools                | ≥2.30   |
| samtools                | ≥1.10   |
| VarScan2                | 2.3.9   |
| TRF                     | 4.09    |
| MUGSY                   | 1.2.3   |
| awk                     | system  |
| fastp                   | ≥0.22   |
| biopython               | ≥1.79   |

Install using conda or system package managers as needed.

---

**The detailed protocols were as follows:**  


### I. Intra-population Analysis (_de novo_ MTD detection)  
Read-based detection of MTDs in a single-population deep-sequencing dataset  

### 1. Preprocessing reads

For Illumina or T7 NGS data, we recommend filtering with fastp (or equivalent tools) to remove low-quality bases and short reads.

```bash
fastp -i sample_R1.fq.gz -I sample_R2.fq.gz \
      -o sample.clean.R1.fq.gz -O sample.clean.R2.fq.gz \
      -q 20 -l 50
```
For PacBio HiFi data, basic QC is usually performed upstream and the pipeline can start directly from the alignment step.

### 2. Read alignment

For Illumina/T7 short reads, use bwa mem to align reads to the reference genome, enabling -Y for soft clipping to retain structural information.
```
bwa mem -Y ref.fa sample.clean.R1.fq.gz sample.clean.R2.fq.gz \
    | samtools view -b - \
    | samtools sort -o sample.sort.bam

samtools index sample.sort.bam
```

For PacBio HiFi data, we recommend pbmm2 with the --preset HiFi option.
```
pbmm2 align ref.fa input.bam output.aln.bam \
    --preset HiFi --sort -j $(nproc)

samtools index output.aln.bam
```
### 3. INDEL calling

Use samtools mpileup to generate pileup files and VarScan2 pileup2indel to call INDELs. Parameters are kept permissive to retain low-frequency de novo events.
```
samtools mpileup -d 10000 -f ref.fa sample.sort.bam > sample.mpileup.txt

java -Xmx20g -jar VarScan.v2.3.9.jar pileup2indel \
    --min-var-freq 0 \
    --p-value 0.99 \
    --min-avg-qual 0 \
    --min-coverage 1 \
    --min-reads2 1 \
    --output-indel sample.indel.txt < sample.mpileup.txt
```

### 4. MTD detection using MTDfinder-Pro (TD/MTD identification)

This step identifies true tandem-duplication events (1×→2×) and determines whether each event is microhomology-mediated (MTD) or non-MTD (NMTD).  
Our script first scans all INDELs, extracts candidate insertion events, validates tandem-duplication structure using TRF, and then performs microhomology detection using a dynamic programming algorithm.

**Important:**  
`Analysis/TRFconvert_script.py` must be in the **same directory** as `Intra_population_denovoTD.sh`.

Run the script as follows:

```bash
bash Analysis/Intra_population_denovoTD.sh ref.fa sample.indel.txt output_prefix
```
ref.fa — reference genome

sample.indel.txt — VarScan2 INDEL file

output_prefix — prefix of the output file

Final output file:
<species_name>_MTD_trf_final.tsv

Output columns:

Chromosome

Insertion start position

Insertion alternate sequence

Insertion sequence

Insertion length

Reference unit 1

Reference unit 2

MHA length (bp)

### II. Intra-species Analysis (fixed/polymorphic MTD detection)

This module detects tandem-duplication events fixed in some isolates but absent in others, enabling the study of MTD polymorphism within a species.

### 1. Whole-genome alignment using MUGSY

Use MUGSY to align isolate genomes to a reference genome and generate a multi-genome MAF file:
```
mugsy --directory <output_dir> \
      --prefix <species_name> \
      <ref_fna> <isolate_genome.fna>
```

<ref_fna> — species reference genome (FASTA)

<isolate_genome.fna> — isolate assembly

<species_name>.maf will be generated
### 2. Extract species-specific insertions (candidate TDs)

Run our extraction script to identify insertion mutations (≥8 bp) from the MAF file:
```
python identify_insertion.py <species_name>.maf <ref_id> <output_insert.tsv>
```

<ref_id> — reference genome ID within the MAF

<output_insert.tsv> — list of all insertions in this isolate

### 3. Detect fixed TD/MTD events among isolates

Use Analysis/Intra_species_fixTD.sh to annotate tandem-duplication structure, evaluate microhomology, and quantify allele counts across isolates.

Before running, edit the following variables inside the script:
```
ref_id=""         ## reference genome ID
ref_fna=""        ## reference genome fasta file
species_name=""   ## species name
indelfile_path="" ## path to insertion/maf files
```

Then run:
```
bash Analysis/Intra_species_fixTD.sh
```

Final output file:
<species_name>_MTD_trf_final.tsv

Output columns:
1–8: same as intra-population module
9. Number of isolates carrying a 2× allele
10. Number of isolates carrying a 3× allele

### III. Species-level Analysis (reference genome TRF-based scanning)

This module characterizes tandem-duplication structures across species using reference genomes.
It reflects long-term evolutionary retention and decay of microhomology-mediated events.

### 1. Run TRF on the reference genome

Use TRF to detect tandem repeats:
```
trf ref.fa 2 7 7 80 10 30 1000 -f -d -m
```
This command outputs a .dat file containing all detected tandem repeats.

### 2. Convert TRF .dat output to TSV format

Use our conversion script:
```
python Analysis/TRFconvert_script.py ref.fa.dat > ref.trf.tsv
```

The TSV file includes repeat coordinates and repeat-unit structures.

### 3. Identify tandem-duplication candidates (length & copy-number filtering)

Example: Keep only repeat units ≥8 bp and unit match% == 100% :
```
awk 'BEGIN { count=0 } $4>7 && $7>=100 && $18 >= 2 && $18 < 3 { count++ } END {  print count }' "$tr2_file"
```
---
These steps constitute the core usage of MTDfinder-Pro at intra-population, intra-species, and species-level scales, with parameters adaptable to different species and sequencing platforms.

## Citing this work

