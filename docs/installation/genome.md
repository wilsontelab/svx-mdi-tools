---
title: Download genomes
parent: Installation and usage
has_children: false
nav_order: 30
published: true
---

## {{page.title}}

The svx-mdi-tools code suite uses 
[Illumina iGenomes](https://support.illumina.com/sequencing/sequencing_software/igenome.html) 
for consistent formatting of reference genome FASTA and aligner index files,
plus additional metadata files from the UCSC genome browser and the ENCODE project.

Specifically, when called as follows:

```sh
mdi <pipeline> <action> --genomes-dir /path/to/genomes --genome hg38
```

pipelines in the svx-mdi-tools suite expect to find folder structure:

```sh
tree -L 4 /path/to/genomes/iGenomes/Homo_sapiens
/path/to/genomes/iGenomes/Homo_sapiens
└── UCSC
    └── hg38
        ├── Annotation
        │   ├── Archives
        │   ├── Genes -> Archives/archive-2015-08-14-08-18-15/Genes
        │   ├── Genes.gencode -> Archives/archive-2015-08-14-08-18-15/Genes.gencode
        │   └── SmallRNA -> Archives/archive-2015-08-14-08-18-15/SmallRNA
        └── Sequence
            ├── AbundantSequences
            ├── Bowtie2Index
            ├── BowtieIndex
            ├── BWAIndex
            ├── Chromosomes
            └── WholeGenomeFasta

tree /path/to/genomes/metadata/hg38
/path/to/genomes/metadata/hg38
├── ENCODE
│   └── hg38-blacklist.v2.bed.gz
├── UCSC
│   ├── gap.txt.gz
│   └── hg38.gc5Base.wigVarStep.gz
```

### Download iGenome

Use the link above to find and copy the URL of your required genome archive 
and use the `download iGenomes` pipeline action to retrieve it to your server.

The following commands will download the UCSC hg38 reference genome
to the current working directory.

```sh
mkdir iGenomes
export URL=http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Homo_sapiens/UCSC/hg38/Homo_sapiens_UCSC_hg38.tar.gz
mdi download iGenomes --urls ${URL} --output-dir ${PWD}/iGenomes --data-name hg38
```

### Download metadata

The following commands will download the associated hg38 metadata files
to the current working directory.

```sh
mkdir -p download/metadata
mdi download metadata --output-dir ${PWD}/download/metadata --data-name hg38
```

### Troubleshooting

If you do not have Singularity on your server, the commands above will fail
until you build the required conda runtime environment:

```sh
mdi download conda --create
```

For more information, see:

```sh
mdi download --help
mdi download iGenomes --help
mdi download conda --help
```

If you already have an appropriately formatted reference genome installation you can use it,
but it must conform to the folder specifications above.
