#!/usr/bin/env perl
# ref. http://bioperl.org/howtos/Features_and_Annotations_HOWTO.html
# Selective: will only record unique ranges: which means no cdS and mRNA basically.
use warnings;
use strict;
use Data::Dumper;
use Bio::Tools::GFF;

# Argument accounting ... say how many the scrit should have.
if(scalar @ARGV != 1) {
    print "Sorry, this script requires exactly one argument\n";
    die;
}

my ($gffio, $feature, $t, $t2, $ti, $ts, $te, $td, $ti2, $val);
# to keep old values
my $tio='';
my $tso=0;
my $teo=0;
my $tdo=0;

my @tt;
my $nfeats=0;
my $nufeats=0; # unique features

$gffio = Bio::Tools::GFF->new(-file => $ARGV[0], -gff_version => 3);
# loop over the input stream
while ($feature = $gffio->next_feature()) { # each feature is a line ... no connection between them, unless you match ID.
	$nfeats++;
	$t=$feature->primary_tag;
	$ti=$feature->seq_id();
	$ts=$feature->location->start;
	$te=$feature->location->end;
	$td=$feature->location->strand;
	# we going to grab the Name subtag ... better than ID subtag which gives same name to gene and CDS
	if($feature->has_tag('Name')) {
		@tt=$feature->get_tag_values('Name');
	}
	# we're only interested in unique ranges:
	# print "olds $tio on strand $tdo: $tso to $teo\n";
	# print "news $ti on strand $td: $ts to $te\n";
	if( $tio ne $ti) {
		$nufeats++;
		print "Beginning new chromosome: primary tag: $t at $ti on strand $td: $ts to $te: nameval=$tt[0]\n";
	} elsif ( ($tso != $ts) & ($teo != $te) & ($tdo != $td) ) {
		$nufeats++;
		print "primary tag: $t at $ti on strand $td: $ts to $te: nameval=$tt[0]\n";
	}
	#bow out by assigning old values:
	$tio=$ti;
	$tso=$ts;
	$teo=$te;
	$tdo=$td;
}

$gffio->close();
print "Num unique features = $nufeats\n";
