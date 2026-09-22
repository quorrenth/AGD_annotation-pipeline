Pipeline for plant gene functional annotation

Install the Database

Blast DB

download the nr, swissprot, trembl, and arabidopsis protein to /data/share/diamond, then makedb:

./diamond makedb --in uniprot_sprot.fasta -d uniprot_sprot
./diamond makedb --in uniprot_trembl.fasta -d uniprot_trembl
./diamond makedb --in Araport11_genes.201606.pep.format.fasta -d ath_pep
./diamond makedb --in nr.fasta -d nr
InterPro DB

download InterPro DB to /data/share/interpro

Install eggnog-mapper

download eggnog-mapper to /data/share/eggnog-mapper

cd /data/share/eggnog-mapper 
sudo python3 setup.py install
Install the pipeline

git clone https://github.com/fjl23/CGD_annotation-pipeline.git
cd CGD_annotation-pipeline
Install iTAK

git clone https://github.com/kentnf/iTAK.git
Install AHRD

git clone https://github.com/groupschoof/AHRD.git
Install python (>=3.8)

sudo apt-get install python3 python3-dev
Install JRE & Perl

sudo apt install default-jre openjdk-11-jre-headless bioperl libbio-perl-perl build-essential cpanminus 
sudo cpanm Bio::SearchIO::blastxml
setup env path

export PATH=$(pwd):$(pwd)/iTAK:$PATH
run this pipeline

Make sure there is a protein.fa file in the current directory

hpg_anno.pl 
