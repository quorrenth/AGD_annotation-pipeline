# Pipeline for plant gene functional annotation

## Install the Database

### Blast DB

Download the nr, swissprot, trembl, and Arabidopsis protein to `/data/share/diamond`, then run:

```bash
./diamond makedb --in uniprot_sprot.fasta -d uniprot_sprot
./diamond makedb --in uniprot_trembl.fasta -d uniprot_trembl
./diamond makedb --in Araport11_genes.201606.pep.format.fasta -d ath_pep
./diamond makedb --in nr.fasta -d nr
```

### InterPro DB

Download InterPro DB to `/data/share/interpro`.

### Install eggnog-mapper

Download eggnog-mapper to `/data/share/eggnog-mapper`.

```bash
cd /data/share/eggnog-mapper
sudo python3 setup.py install
```

## Install the pipeline

```bash
git clone https://github.com/fjl23/CGD_annotation-pipeline.git
cd CGD_annotation-pipeline
```

## Install iTAK

```bash
git clone https://github.com/kentnf/iTAK.git
```

## Install AHRD

```bash
git clone https://github.com/groupschoof/AHRD.git
```

## Install Python (>=3.8)

```bash
sudo apt-get install python3 python3-dev
```

## Install JRE & Perl

```bash
sudo apt install default-jre openjdk-11-jre-headless bioperl libbio-perl-perl build-essential cpanminus
sudo cpanm Bio::SearchIO::blastxml
```

## Setup env path

```bash
export PATH=$(pwd):$(pwd)/iTAK:$PATH
```

## Run this pipeline

Make sure there is a `protein.fa` file in the current directory.

```bash
hpg_anno.pl
```
