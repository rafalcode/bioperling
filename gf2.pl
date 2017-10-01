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
my $nufeats=0; # unique features
my %holol; # hash of list of lists.
my $nchil=0;

$gffio = Bio::Tools::GFF->new(-file => $ARGV[0], -gff_version => 3);
# loop over the input stream
while ($feature = $gffio->next_feature()) { # each feature is a line ... no connection between them, unless you match ID.
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
		# print "Beginning new chromosome: primary tag: $t at $ti on strand $td: $ts to $te: nameval=$tt[0]\n";
		# push($holol{$t}->[$nufeats], [$tt[0], $ti, $td, $ts, $te]);
		# push($holol{$t}[$nufeats], ($tt[0], $ti, $td, $ts, $te]);
		# $holol{$t}[$nufeats]= [$tt[0], $ti, $td, $ts, $te];
		# $holol{$t}++;
		# $holol{$t}[$nufeats]= [$tt[0], $ti, $td, $ts, $te];
		push( @{$holol{$t}}, [$tt[0], $ti, $td, $ts, $te, $nchil]);
		$nufeats++;
		$nchil=0;
	} elsif ( ($tso != $ts) & ($teo != $te) & ($tdo != $td) ) {
		# print "primary tag: $t at $ti on strand $td: $ts to $te: nameval=$tt[0]\n";
		# push($holol{$t}->[$nufeats], [$tt[0], $ti, $td, $ts, $te]);
		# $holol{$t}[$nufeats]= [$tt[0], $ti, $td, $ts, $te];
		push( @{$holol{$t}}, [$tt[0], $ti, $td, $ts, $te, $nchil]);
		$nufeats++;
		$nchil=0;
	} elsif ( ($tso <= $ts) & ($teo >= $te) & ($tdo == $td) ) {
		$nchil++; # increment numb of children
	}
	#bow out by assigning old values:
	$tio=$ti;
	$tso=$ts;
	$teo=$te;
	$tdo=$td;
}
$gffio->close();
print "Num unique features = $nufeats\n";
my $l;
print "HOLOL keys = ";
foreach my $k (keys %holol) {
	$l=scalar @{$holol{$k}};
	print "$k # $l  ";
}
print "\n";
for my $k (keys %holol) {
	for my $a (@{$holol{$k}}) {
		# print "$a->[0] ";
		print "$a->[0],$a->[1],$a->[2],$a->[3],$a->[4], $a->[5]) "
	}
	print "\n";
}

# line snippets region
# each array in hash element will have the following:
# print "$a->[0],$a->[1],$a->[2],$a->[3],$a->[4], $a->[5]) "
