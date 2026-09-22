#!/usr/bin/perl

use strict;
use warnings;
use Cwd qw(cwd);
use FindBin;
use IO::File;
use Bio::SearchIO;

my $usage = qq'
USAGE: $0 

';

# check if protein sequence exist
my $pep = 'protein.fa';
die "[ERR]protein.fa is not exist\n" unless -s $pep;

my $func_anno_fdr = "func_anno";
unless (-e $func_anno_fdr) {
	system("mkdir func_anno");
}

# prepare blast command
my $threads = 90;
my $db_path = "/data/share/diamond";
my $db_sp = "$db_path/uniprot_sprot.dmnd"; 
my $db_nr = "$db_path/nr.dmnd";
my $db_at = "$db_path/ath_pep.dmnd";
my $db_tr = "$db_path/uniprot_trembl.dmnd";
my $blastp_param = "-p $threads --max-target-seqs 500 --outfmt 5 --more-sensitive";

my $blastp_nr_cmd = "diamond blastp --db $db_nr -q $pep -p $threads -e 1e-3 --max-target-seqs 20 --outfmt 5 --more-sensitive -q $pep -o $func_anno_fdr/pep.dia.nr.xml";
my $blastp_tr_cmd = "diamond blastp --db $db_tr $blastp_param -q $pep -o $func_anno_fdr/pep.dia.tr.xml";
my $blastp_sp_cmd = "diamond blastp --db $db_sp $blastp_param -q $pep -o $func_anno_fdr/pep.dia.sp.xml";
my $blastp_at_cmd = "diamond blastp --db $db_at $blastp_param -q $pep -o $func_anno_fdr/pep.dia.at.xml";

# run blast command
if (-s "$func_anno_fdr/pep.dia.nr.xml" || -s "$func_anno_fdr/pep.dia.nr.xml.gz") {
	print "[WARN]the output file pep.dia.nr.xml exist, skip...\n";
} else {
	run_cmd($blastp_nr_cmd);
}

if (-s "$func_anno_fdr/pep.dia.tr.xml" || -s "$func_anno_fdr/pep.dia.tr.xml.gz" ) {
    print "[WARN]the output file pep.dia.tr.xml exist, skip...\n";
} else {
	run_cmd($blastp_tr_cmd);
}

if (-s "$func_anno_fdr/pep.dia.sp.xml" || -s "$func_anno_fdr/pep.dia.sp.xml.gz") {
    print "[WARN]the output file pep.dia.sp.xml exist, skip...\n";
} else {
	run_cmd($blastp_sp_cmd);
}

if (-s "$func_anno_fdr/pep.dia.at.xml" || -s "$func_anno_fdr/pep.dia.at.xml.gz") {
    print "[WARN]the output file pep.dia.at.xml exist, skip...\n";
} else {
    run_cmd($blastp_at_cmd);
}

# generate GO and KEGG result using eggNOG-mapper
my $eggnog_cmd = "python3 /data/share/eggnog-mapper/emapper.py -i $pep --output $func_anno_fdr/eggNOG -m diamond --cpu $threads";
if (-s "$func_anno_fdr/eggNOG.emapper.annotations") {
        print "[WARN]the eggNOG result exist, skip...\n";
} else {
        run_cmd($eggnog_cmd);
}

# run iTAK
my $itak_cmd = "iTAK.pl -o $func_anno_fdr/itak $pep";
if (-s "$func_anno_fdr/itak") {
	print "[WARN]iTAK result exist, skip ...\n";
} else {
	run_cmd($itak_cmd);
}

# run interproscan
my $interpro_db = "interproscan-5.55-88.0";
#my $interpro_db = "interproscan-5.45-80.0";
my $interpro_cmd = "/data/share/interpro/$interpro_db/interproscan.sh -goterms -dp -i $pep -f XML -d $func_anno_fdr 1>$func_anno_fdr/ipr.report.txt 2>&1";

if (-s "$func_anno_fdr/$pep.xml" || -s "$func_anno_fdr/$pep.xml.gz") {
	print "[WARN]the InterPro XML exist, skip...\n";
} else {
	run_cmd($interpro_cmd);
}

# generate ahrd analysis 
if (-s "$func_anno_fdr/ahrd.csv") {
	print "[WARN]the ahrd.csv exist, skip AHRD...\n";
}
else {
	if (-s "$func_anno_fdr/pep.dia.at.xml" && -s "$func_anno_fdr/pep.dia.sp.xml" && -s "$func_anno_fdr/pep.dia.tr.xml") {
		convert_xml_to_tab("$func_anno_fdr/pep.dia.at.xml", "$func_anno_fdr/pep.dia.at.tab");
		convert_xml_to_tab("$func_anno_fdr/pep.dia.sp.xml", "$func_anno_fdr/pep.dia.sp.tab");
		convert_xml_to_tab("$func_anno_fdr/pep.dia.tr.xml", "$func_anno_fdr/pep.dia.tr.tab");
		ahrd_yml($func_anno_fdr, "ahrd.yml", "ahrd.csv");
	}
       	elsif (-s "$func_anno_fdr/pep.dia.at.tab" && -s "$func_anno_fdr/pep.dia.sp.tab" && -s "$func_anno_fdr/pep.dia.tr.tab") {
		ahrd_yml($func_anno_fdr, "ahrd.yml", "ahrd.csv");
	}	
	else {
		print "[WARN]the at/sp/tr blast xml files are not exist, skip AHRD...\n"
	}
}

=head2
 ahrd_yml: perform ahrd analysis
=cut
sub ahrd_yml {
	my ($func_anno_fdr, $ahrd_yml, $ahrd_csv) = @_;
	my $working_dir = cwd;

	my $ahrd_yml_file = "$working_dir/$func_anno_fdr/$ahrd_yml";
	my $ahrd_csv_file = "$working_dir/$func_anno_fdr/$ahrd_csv";

	my $func_anno_at = "$working_dir/$func_anno_fdr/pep.dia.at.tab";
	my $func_anno_sp = "$working_dir/$func_anno_fdr/pep.dia.sp.tab";
	my $func_anno_tr = "$working_dir/$func_anno_fdr/pep.dia.tr.tab";

    my $ahrd_path = ${FindBin::RealBin}."/AHRD-3.3.3";

    my $fhy = IO::File->new(">".$ahrd_yml_file) || die $!;
    print $fhy qq'
proteins_fasta: $working_dir/$pep
token_score_bit_score_weight: 0.468
token_score_database_score_weight: 0.2098
token_score_overlap_score_weight: 0.3221
output: $ahrd_csv_file
blast_dbs:
  swissprot:
    weight: 653
    description_score_bit_score_weight: 2.717061
    file: $func_anno_sp
    database: /data/share/diamond/uniprot_sprot.fasta
    blacklist: $ahrd_path/test/resources/blacklist_descline.txt
    filter: $ahrd_path/test/resources/filter_descline_sprot.txt
    token_blacklist: $ahrd_path/test/resources/blacklist_token.txt

  trembl:
    weight: 904
    description_score_bit_score_weight: 2.590211
    file: $func_anno_tr
    database: /data/share/diamond/uniprot_trembl.fasta
    blacklist: $ahrd_path/test/resources/blacklist_descline.txt
    filter: $ahrd_path/test/resources/filter_descline_trembl.txt
    token_blacklist: $ahrd_path/test/resources/blacklist_token.txt

  tair:
    weight: 854
    description_score_bit_score_weight: 2.917405
    file: $func_anno_at
    database: /data/share/diamond/Araport11_genes.201606.pep.format.fasta
    fasta_header_regex: "^>(?<accession>\\\\S+)\\\\s+(?<description>.+?)\$"
    short_accession_regex: "^(?<shortAccession>.+)\$"
    blacklist: $ahrd_path/test/resources/blacklist_descline.txt
    filter: $ahrd_path/test/resources/filter_descline_tair.txt
    token_blacklist: $ahrd_path/test/resources/blacklist_token.txt
';

	$fhy->close;

    # run AHRD
    run_cmd("java -Xmx150g -jar $ahrd_path/dist/ahrd.jar $ahrd_yml_file") unless -s $ahrd_csv_file;
}

=head2
 run_cmd
=cut
sub run_cmd {
	my $cmd = shift;
	print "[CMD]:$cmd\n";
	system($cmd) && die "[ERR]CMD:$cmd\n";
}

=head2
 convert_xml_to_tab
=cut
sub convert_xml_to_tab {

	my ($blast_in, $tab_out) = @_;

	if (-s $tab_out) {
		print "[WARN]the convert blast $tab_out exist, skip ...\n";
		return 0;
	}

	my $fht = IO::File->new(">$tab_out") || die $!;
	my $searchio = Bio::SearchIO->new(-format => "blastxml", -file   => $blast_in );
	while (my $result = $searchio->next_result) {
		my $query_name   = $result->query_name;
		my $query_length = $result->query_length;
		my $query_def    = $result->query_description();
		while (my $hit = $result->next_hit) {
			my $hit_name = $hit->name;
			my $hit_desc = $hit->description;
			my $hit_len  = $hit->length;
			while (my $hsp = $hit->next_hsp) {
				my $identity     = sprintf "%0.2f", $hsp->frac_identical * 100;
				my $gaps         = $hsp->gaps;
				my $mismatch0    = $hsp->seq_inds('hit','nomatch');
				my $align_len    = $hsp->length;
				my $query_start  = $hsp->start('query');
				my $query_end    = $hsp->end('query');
				my $hit_start    = $hsp->start('hit');
				my $hit_end      = $hsp->end('hit');
				my $evalue       = $hsp->evalue;
				my $bit_score    = $hsp->bits;
				print $fht "$query_def\t$hit_name\t$identity\t$align_len\t$mismatch0\t$gaps\t$query_start\t$query_end\t$hit_start\t$hit_end\t$evalue\t$bit_score\n";
			}
    	}
	}
	$fht->close;
}
