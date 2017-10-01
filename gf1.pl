#!/usr/bin/env perl
# ref. http://bioperl.org/howtos/Features_and_Annotations_HOWTO.html
# Selective printing now please.
use warnings;
use strict;
use Data::Dumper;
use Bio::Tools::GFF;

# Argument accounting ... say how many the scrit should have.
if(scalar @ARGV != 1) {
    print "Sorry, this script requires exactly one argument\n";
    die;
}

my ($gffio, $feature, $t, $t2, $ti, $ts, $te, $ti2, $val);
my $nfeats=0;

$gffio = Bio::Tools::GFF->new(-file => $ARGV[0], -gff_version => 3);
# loop over the input stream
while ($feature = $gffio->next_feature()) { # each feature is a line ... no connection between them, unless you match ID.
	$nfeats++;
	$t=$feature->primary_tag;
	$ti=$feature->seq_id();
	$ts=$feature->location->start;
	$te=$feature->location->end;
	$ti2=$feature->location->end;
	print "primary tag: $t at $ti: $ts to $te\n";
	for $t2 ($feature->get_all_tags) { # subtags, actually.
		print "  tag: ", $t2, "\n";
	    for my $val ($feature->get_tag_values($t2)) {
			print "    value: ", $val, "\n";
		}
	}
}

# more compact?
# # nah, this needs an array.
# my @cdsfs = grep { $_->primary_tag eq 'CDS' } $gffio->next_feature();
# my $cdsfsz=scalar @cdsfs;


$gffio->close();
print "Num features = $nfeats\n";
# print "Sz cds array = $cdsfsz\n";
