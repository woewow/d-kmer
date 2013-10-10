#!/usr/bin/perl -w
#
use strict;
use warnings;

my $db = shift or die "Usage: $0 <map file> <gene:aapos| file>\n";

my $query = shift or die "Usage: $0 <map file> <gene:aapos| file>\n";

my %hash = ();
open IN, $db or die $!;
my $pre = "";
my $pre_len = 0;
while (<IN>) {
	chomp;
	my @ele = split /\s+/, $_;
	if ($pre ne "") {
		if ($ele[0] eq $pre) {
			if (length $ele[-1] > $pre_len) {
				$hash{$ele[0]} = join "\t", @ele[1 .. 9];
				$pre_len = length $ele[-1];
			}
		} else {
			$hash{$ele[0]} = join "\t", @ele[1 .. 9];
			$pre = $ele[0];
			$pre_len = length $ele[-1];
		} 
	} else {
		$hash{$ele[0]} = join "\t", @ele[1 .. 9];
		$pre = $ele[0];
		$pre_len = length $ele[-1];
	}
}
close IN;

if (-f $query) {
	open IN, $query or die $!;
	while (<IN>) {
		chomp;
		&trans_pos($_, \%hash);
	}
	close IN;
} else {
	&trans_pos($query, \%hash);
}

1;

sub trans_pos {
	my $gene_aapos = shift;
	my $genomedb = shift;
	my ($gene, $aapos) = split /:/, $gene_aapos;
	my $rec = $$genomedb{$gene};
	my @trans = split /\s+/, $rec;
	my @exons = split /,/, $trans[7];
	my @exone = split /,/, $trans[8];
	my @cdspos = ();
	my ($s, $e);
	for (my $i = 0; $i < $trans[6]; $i ++) {
		next if ($exone[$i] < $trans[4]);
		last if ($trans[5] < $exons[$i]);
		if ($trans[4] >= $exons[$i]) {
			$s = $trans[4];
		} else {
			$s = $exons[$i];
		}
		if ($trans[5]<=$exone[$i]) {
			$e = $trans[5];
		} else {
			$e = $exone[$i];
		}
		push @cdspos, ($s .. $e);
	}
	my $aalen = scalar(@cdspos) / 3;
	if ($trans[1] eq '+') {
		print join("\t", $gene_aapos, "=>", @trans[0..1], $cdspos[$aapos*3 - 3]), "\n";
	} else {
		print join("\t", $gene_aapos, "=>", @trans[0..1], $cdspos[($aalen-$aapos)*3]), "\n";
	}
}
