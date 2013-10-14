#!/usr/bin/perl -w
#
use strict;
use warnings;
use Getopt::Std;

my %opts;
getopts('hk:', \%opts);

usage() if ($opts{h});
my $K = $opts{k} || 31;

my $db = shift or die usage();
my $rds = shift or die usage();

my %db_hash = ();


open IN, $db or die $!;
while (<IN>) {
	chomp;
	my $len = length $_;
	$db_hash{substr($_, $len/2-$K/2, $K)} = 1;
}
close IN;

if ($rds=~/.gz$/) {
	open IN, "gzip -dc $rds |" or die $!;
} else {
	open IN, $rds or die $!;
}

while (<IN>) {
	chomp;
	if ($. % 4 == 2) {
		query_db($_, \%db_hash, $K);
	}
}
close IN;


foreach my $k (%db_hash) {
	print join("\t", $k, $db_hash{$k}), "\n";
}
1;

sub usage {
	print STDERR << "EOF";
This program calculates known kmer frequencies using raw data.
Usage: $0 [-hk] <db_seq> <fq_file(.gz)>
	-h           This help message
	-k  <INT>    K-mer size [31]
EOF
	exit;
}

sub query_db {
	my ($seq, $db_ref, $k) = @_;
	my %db = %$db_ref;
	for (my $i = 0; $i <= length($seq)-$K; $i++) {
		$db{substr($seq, $i, $K)} ++ if ($db{substr($seq, $i, $K)});
	}

	my $rev_seq = dna_rev_seq($seq);

	for (my $i = 0; $i <= length($rev_seq)-$K; $i++) {
		$db{substr($rev_seq, $i, $K)} ++ if ($db{substr($rev_seq, $i, $K)});
	}
}

sub dna_rev_seq {
	my $seq = shift;
	my $ret = reverse $seq;
	$ret = tr/ACGTacgt/TGCAtgca/;
	return $ret;
}
