# MTDfinder-pro
MTDfinder-pro provides a powerful approach to detect MicroHomologyMediatedTandemDuplications(MTD) in multiple species.  
This pipeline is used to analyze DEtail-seq data.
## Requirements
|[python]|[bwa]|[bedtools]|[samtools]|[Tandem Repeat Finder]|
|---|---|---|---|---|
|[**VarScan2**]|[**MUGSY**]|[**awk**]|[**fastp**]|[**biopython**]|
## Documentation
First, duplicated reads which had same sequences for both forward and reverse reads were removed. And our own scripts RemoveSamReads.py were used to remove duplicated reads. Then, reads were aligned to the reference genome with Bowtie2 using --local settings. For visualization and spliting strand, the aligned reads files (BAM) were converted to single base bigWig file with 1 bp bins using bamCoverage from deepTools. Finally, poisson distribution is used for Hotspots calling, which is implemented by our own script Hotspotcalling.py.  
  
**The detailed protocols were as follows:**
🧬 MTDfinder-Pro

English version is followed by the Chinese version.

MTDfinder-Pro is a comprehensive computational pipeline for detecting Microhomology-Mediated Tandem Duplications (MTDs) across three evolutionary levels—intra-population, intra-species, and species-level.
The pipeline supports both Illumina short reads and PacBio HiFi long reads, and integrates INDEL calling, 1×→2× tandem-duplication confirmation, and microhomology detection.

A full description of the MTDfinder-Pro algorithm is provided in the manuscript (Revised Methods, lines xx–xx), and demonstrated in Author Response Figure 7.

📦 Requirements 
Required tools	Version (recommended)
python3	≥3.8
bwa	≥0.7
pbmm2 (for PacBio HiFi)	latest
bedtools	≥2.30
samtools	≥1.10
VarScan2	2.3.9
TRF (Tandem Repeat Finder)	4.09
MUGSY (for multi-genome alignment)	1.2.3
awk	system
fastp	≥0.22
biopython	≥1.79




### 1.Removing duplication reads
```
$ RemoveSamReads.py test_R1.fastq.gz test_R2.fastq.gz test
```

### 2.Alignment
```
$ bowtie2-build --threads 20 RefGenome.fasta RefGenome
$ bowtie2 --local --phred33 -p 20 -t -x RefGenome -1 test_dup_R1.fq.gz -2 test_dup_R2.fq.gz\
2>test_align.info|samtools view -bS -1 |samtools sort -@ 10 -m 5G -l 9 -o test.sort.bam
```
### 3.Converting BAM to single base bigWig and splitting strand
--minMappingQuality 30 may sometimes be needed
```
$ samtools index test.sort.bam
$ bamCoverage -v -p 60 -b test.sort.bam -o test_Crick.bw --binSize 1\
--Offset 1 --samFlagInclude 128 --filterRNAstrand forward
$ bamCoverage -v -p 60 -b test.sort.bam -o test_Watson.bw --binSize 1\
--Offset 1 --samFlagInclude 128 --filterRNAstrand reverse
```
### 4.Hotspot calling
```
usage: Hotspotcalling.py [-h] [--bw BIGWIG] [--strand {fwd,rev}]
                         [--qvalue QVALUE] [--res REST] [--prefix PREFIX]

optional arguments:
  -h, --help          show this help message and exit
  --bw BIGWIG
  --strand {fwd,rev}
  --qvalue QVALUE
  --filter FILTER     A list of comma-delimited chromosome names, containing those chromosomes that 
                      should be excluded for computing the hotspotcalling and normalization.
  --res REST          Restriction sites bed file
  --norm              normalizing bw file
  --debug
  --prefix PREFIX

$ Hotspotcalling.py --bw test_Watson.bw --strand fwd --qvalue 0.01 --prefix test_Watson --norm
$ Hotspotcalling.py --bw test_Crick.bw --strand rev --qvalue 0.01 --prefix test_Crick --norm

```
If restriction enzymes are used, you may need to use the --res parameter to add restriction site information

test_Crick_basepeaks.bed and test_Watson_basepeaks.bed are the hotspot results in bed format.
test_Crick_basepeaks.xls and test_Watson_basepeaks.xls are the detailed hotspot results.

**If you want to process multiple samples in batches, please refer to DEtailPip.sh**

## Citing this work
Xu, W., Liu, C., Zhang, Z., Sun, C., Li, Q., Li, K., Jiang, H., Li, W., & Sun, Q. (2023). [DEtail-seq is an ultra-efficient and convenient method for meiotic DNA break profiling in multiple organisms.](https://pubmed.ncbi.nlm.nih.gov/36723795/) Science China. Life sciences, 66(6), 1392–1407.
