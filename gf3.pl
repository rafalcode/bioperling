#!/usr/bin/env perl
# ref. http://bioperl.org/howtos/Features_and_Annotations_HOWTO.html
# gf3.pl: is able to deal with smaller gffs ... took out conditional on strand direction to allow ARS to have children too
use warnings;
use strict;
use Data::Dumper;
use Bio::Tools::GFF;

# sub routine ... printing ordered hash
sub printohchrom { # print ordered hash for hostfile
	# lexicographic sort on chromosome, secondary on start site
	my $extent;
    foreach my $el (sort { $_[0]{$a}->[1] cmp $_[0]{$b}->[1] or $_[0]{$a}->[3] <=> $_[0]{$b}->[3] } keys %{$_[0]}) {
		$extent=$_[0]{$el}->[4]-$_[0]{$el}->[3];
		print "$el\t$_[1]{$el}\t$_[0]{$el}->[0]\t$_[0]{$el}->[1]\t$_[0]{$el}->[2]\t$_[0]{$el}->[3]\t$extent\n";
    }
}

# Argument accounting ... say how many the scrit should have.
if(scalar @ARGV != 1) {
    print "Sorry, this script requires exactly one argument\n";
    die;
}

my ($gffio, $feature, $t, $t2, $ti, $ts, $te, $td, $ti2, $val);
# to keep old values
my $to='';
my $tio='';
my $tto='';
my $tso=0;
my $teo=0;
my $tdo=0;

my @tt;
my %hofn; # hash of name of the feature: most often, gene names.
my %honc; # hash of name of the feature: count of children
my %hopf; # hash of primary features.
my $nchil=0;

$gffio = Bio::Tools::GFF->new(-file => $ARGV[0], -gff_version => 3);
# loop over the input stream
while ($feature = $gffio->next_feature()) { # each feature is a line ... no connection between them, unless you match ID.
	$t=$feature->primary_tag;
	next if($t eq 'chromosome'); # not interested in full chromosomes.
	$ti=$feature->seq_id();
	$ts=$feature->location->start;
	$te=$feature->location->end;
	$td=$feature->location->strand;
	# we going to grab the Name subtag ... better than ID subtag which gives same name to gene and CDS
	if($feature->has_tag('Name')) {
		@tt=$feature->get_tag_values('Name');
	}
	print "curr dataline: $tt[0] / $t / $ti / $td / $ts / $te\n";
	print "prev dataline: $tto / $to / $tio / $tdo / $tso / $teo\n";
	# we're only interested in unique ranges:
	# print "olds $tio on strand $tdo: $tso to $teo\n";
	# print "news $ti on strand $td: $ts to $te\n";
	if( $tio ne $ti) { # if the chromosome is different, we definitely want to start a new dict.
		$hofn{$tt[0]}=[$t, $ti, $td, $ts, $te];
		$hopf{$t}++;
		$honc{$tto}=$nchil if $tto; # primfeat is new ... assign this value to previous.
		print "blk1: tto is $tto assigned $nchil while tt0 is $tt[0]\n" if defined $tto;
		$tto=$tt[0];
		$to=$t;
		$tio=$ti;
		$tso=$ts;
		$teo=$te;
		$tdo=$td;
		$nchil=0;
	} elsif ( $t eq 'gene') { # importatn enough to have its own branch conditional
		$hofn{$tt[0]}=[$t, $ti, $td, $ts, $te];
		$hopf{$t}++; # to count the primary features.
		$honc{$tto}=$nchil if $tto;
		print "blk2: tto is $tto assigned $nchil while tt0 is $tt[0]\n" if defined $tto;
		$tto=$tt[0];
		$to=$t;
		$tio=$ti;
		$tso=$ts;
		$teo=$te;
		$tdo=$td;
		$nchil=0;
	} elsif ( ($to ne $t) & ($tso < $ts) & ($teo < $te) ) { # chromosome is the same but feature is new and range is different: ignore identical ranges with differnet names.
		$hofn{$tt[0]}=[$t, $ti, $td, $ts, $te];
		$hopf{$t}++; # to count the primary features.
		$honc{$tto}=$nchil if $tto;
		print "blk2: tto is $tto assigned $nchil while tt0 is $tt[0]\n" if defined $tto;
		$tto=$tt[0];
		$to=$t;
		$tio=$ti;
		$tso=$ts;
		$teo=$te;
		$tdo=$td;
		$nchil=0;
	} elsif ( ($to ne $t) & ($tso <= $ts) & ($teo >= $te) ) {
		print "Now just count children: $tt[0] prev: $tto: ($tso <= $ts) ($teo >= $te) ($tdo == $td)\n";
		$nchil++;
	}
}
print "last primfeat $t\n";

# if($t ne 'chromosome') { # not interested in full chromosomes.
# 	if( $tio ne $ti) { # if the chromosome is different, we definitely want to start a new dict.
# 		$hofn{$tt[0]}=[$t, $ti, $td, $ts, $te];
# 		$hopf{$t}++;
# 		$honc{$tto}=$nchil if defined $tto; # primfeat is new ... assign this value to previous.
# 		print "blk1: tto is $tto assigned $nchil while tt0 is $tt[0]\n" if defined $tto;
# 		$tto=$tt[0];
# 		$nchil=0;
# 	} elsif ( ($to ne $t) & ($tso != $ts) & ($teo != $te) ) { # chromosome is the same but range is different: ignore identical ranges with differnet names.
# 		$hofn{$tt[0]}=[$t, $ti, $td, $ts, $te];
# 	} elsif ( ($to ne $t) & ($tso <= $ts) & ($teo >= $te) & ($tdo == $td) ) {
# 		print "Now just count children: $tt[0] prev: $tto: ($tso <= $ts) ($teo >= $te) ($tdo == $td)\n";
# 		$nchil++;
# 	}
# }
# assign final nchil value
$honc{$tto}=$nchil;


$gffio->close();
my $l;
my $nufeats=0;
print "hopf keys =\n";
foreach my $k (keys %hopf) {
	print "$k # $hopf{$k}\n";
	$nufeats += $hopf{$k};
}
print "\n>>> In total, quantity unique features = $nufeats\n";

# for my $k (keys %hofn) {
# 	print "$k\t$hofn{$k}->[0]\t$hofn{$k}->[1]\t$hofn{$k}->[2]\t$hofn{$k}->[3]\t$hofn{$k}->[4]\n";
# }
printohchrom(\%hofn, \%honc);

# line snippets region
# each array in hash element will have the following:
# print "$a->[0],$a->[1],$a->[2],$a->[3],$a->[4], $a->[5]) "
