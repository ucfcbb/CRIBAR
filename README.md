# CRIBAR

CRIBAR is a fast and flexible sgRNA design tool for CRISPR-based imaging. For enhanced brightness, CRISPR-based imaging systems adapt multiple CRISPR
single guide RNA (sgRNA) scaffolds to mark binding sites for fluorescent proteins. Therefore, CRISPR-based imaging systems require clustered binding sites in the target regions. In regions outside the target regions, the sgRNAs should not form cluster off-target copies that may generate
excessive noise at undesired loci. CRIBAR can generate an optimal set of sgRNA to reduce the cost and increase the efficiency of the CRISPR-based imaging experiment.

## Formulation

On-target activity score is the weighted sum of the number of binding sites. For each binding site, the score is between 0 and 1. Most of the perfect binding sites have activity score of 1. The score may vary based on the number, type, and position of the mismatches/insertions/deletions.  
* Formulation 1: Given the minimum on-target activity score (constraint), output the minimal set of sgRNAs (objective)   
 <!-- In a sub-optimal case, we can use 20 sgRNAs to have 20 on-target activity score. CRIBAR may achieve the same on-target activity score by using fewer sgRNAs. This can help to save the cost of the experiment. -->
* Formulation 2: Given the maximum number of sgRNAs (constraint), output the maximum on-target activity score (objective)   
 <!-- In a sub-optimal case, by using 20 sgRNAs, we can have 20 on-target activity score. CRIBAR may generate a set of 20 sgRNAs to achieve the a higher on-target activity score. This can help to maximize the brightness of the target region in CRISPR-based imaging experiment. -->
  
## Dependencies

Python 3.8, pysam 0.21.0, PuLP 2.7.0: ILP solver dependencies\
bedtools: To get fasta for the target regions\
bowtie: To search for off-target sites\
CRISPRitz(Optional): To exhaustively search for off-target sites

**Install prelimineries in Ubuntu**

For Ubuntu, the dependencies can be installed by using the following commands:
<!-- If you are using Ubuntu, you can install dependentcies by using the following commands: -->

```
sudo apt install bedtools bowtie glpk-utils
```

**Install required Python libraries**

```
pip install -r requirements.txt
```

**Collect bowtie indexes**

Download hg38 bowtie index from ftp://ftp.ccb.jhu.edu/pub/data/bowtie_indexes/GRCh38_no_alt.zip. Put the *GRCh38_noalt_as folder* in the root directory of CRIBAR. Make sure you have the following files in CRIBAR/GRCh38_noalt_as/
```
GRCh38_noalt_as.1.ebwt
GRCh38_noalt_as.2.ebwt
GRCh38_noalt_as.3.ebwt
GRCh38_noalt_as.4.ebwt
GRCh38_noalt_as.fa
GRCh38_noalt_as.rev.1.ebwt
GRCh38_noalt_as.rev.2.ebwt
```

If you only have the bowtie index files, you can generate the .fa file by running the following command in CRIBAR/GRCh38_noalt_as/ directory:
```
bowtie-inspect GRCh38_noalt_as > GRCh38_noalt_as.fa
```
Make sure the Genome prefix provided with *--genome_prefix* parameter exactly matches with the bowtie index file prefixes.

**Collect CRISPRitz and provide CRISPRitz directory location as parameter**

Download CRISPRitz from https://anaconda.org/bioconda/crispritz/files and put the extracted files inside lib directory under the root directory of CRIBAR. For example, if linux-64_crispritz-2.6.6-py39h68928f9_1.tar.bz2 file is downloaded from the mentioned source, the extracted files and folder should be placed inside CRIBAR/lib/linux-64_crispritz-2.6.6-py39h68928f9_1/ directory. The CRISPRitz version used here is 2.6.6 and the default CRISPRitz directory is set to the one mentioned as example here. The directory can be provided with *--crispritz_dir* parameter.

## Usage steps

**Example command**

Open terminal from the folder that contains cribar.py and run the following command:

```
python3 cribar.py --chr chr20 --start 20000001 --end 20005001 --work_dir ../work/ --tar_seq_dir ../tar_seq/ --genome_prefix ../GRCh38_noalt_as/GRCh38_noalt_as --src_dir ./ --len 20 --mismatch 3 --formulation 1 --constraint 5 --pam_seq NGG --off_target_window 0 --off_target_ratio 0.2 --min_gap 30 --excluded_substring AAAA,TTTT
```

## Parameter description

**genome:**	The reference genome used in CRISPR imaging experiment design.  
**chr**:	The chromosome of the target region.  
**start**:	The starting position of target region.  
**end**:	The ending position of target region.  
**len**:	The length of gRNA (not include pam).  
**formulation**:	1 - Given the minimum on-target activity score, output the minimal set of sgRNAs.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2 - Given the maximum number of sgRNAs, output a set of sgRNAs with the maximum on-target activity score.  
**constraint:**	When formulation=1, constrain is the minimum on-target activity score. When formulation=2, constrain is the maximum gRNA number.  
**pam_seq:**	Protospacer adjacent motif. Default is NGG.  
**mismatch:**	The number of mismatches allowed in each gRNA.  
**off_target_window:**	The window size that CRIBAR uses to check the off-targets. Default value is the same length of the target region.  
**off_target_ratio:**	The threshold of (off-target binding site density) / (on-target binding site density).  
**min_gap:**	The minimum distance between two binding sites.  
**excluded_substring:**	Excluded substrings separated by a comma. Example: AAAA,TTTT.  

## CONTACTS

For bug reports or comments please contact xiaolichen.cs@gmail.com or mahfuz@ucf.edu.
