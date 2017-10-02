#!/usr/bin/env perl
# ref. http://bioperl.org/howtos/Features_and_Annotations_HOWTO.html
# gf4.pl: printing out  the features that are children to the main features.
use warnings;
use strict;
use Data::Dumper;
use Bio::Tools::GFF;

sub printohchrom { # print ordered hash three args ..
	# lexicographic sort on chromosome, secondary on start site
	my ($extent, $i);
	print "NAME\tTYPE\tCHROM\tSTRAND\tSTART\tEND\tEXTENT\tNUMCHIL\t(TYPES OF CHILS)\n";
    foreach my $el (sort { $_[0]{$a}->[1] cmp $_[0]{$b}->[1] or $_[0]{$a}->[3] <=> $_[0]{$b}->[3] } keys %{$_[0]}) {
		$extent=$_[0]{$el}->[4]-$_[0]{$el}->[3];
		# print "$el\t$_[1]{$el}\t$_[0]{$el}->[0]\t$_[0]{$el}->[1]\t$_[0]{$el}->[2]\t$_[0]{$el}->[3]\t$extent\n";
		print "$el\t$_[0]{$el}->[0]\t$_[0]{$el}->[1]\t$_[0]{$el}->[2]\t$_[0]{$el}->[3]\t$_[0]{$el}->[4]\t$extent\t$_[1]{$el}";
		for($i=0;$i<$_[1]{$el};$i++) {
			print "\t$_[2]{$el}->[$i]";
		}
		print "\n";
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
my %hochil;
my $withinregion=0;
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
	if($feature->has_tag('ID')) { # I was using Name  s it distinguishes between mRNA and CDS, but ID is usually more dependable
		@tt=$feature->get_tag_values('ID'); 
	}
	if( $tio ne $ti) { # if the chromosome is different, we definitely want to start a new dict.
		$hofn{$tt[0]}=[$t, $ti, $td, $ts, $te];
		$hopf{$t}++;
		$honc{$tto}=$nchil if $tto; # primfeat is new ... assign this value to previous.
		$tto=$tt[0];
		$to=$t;
		$tio=$ti;
		$tso=$ts;
		$teo=$te;
		$tdo=$td;
		$nchil=0;
	} elsif ( $t eq 'region') { # these usually encapsulate a gene, but you lose the gene name ... give them their own line or filter out.
		$hofn{$tt[0]}=[$t, $ti, $td, $ts, $te];
		$hopf{$t}++; # to count the primary features.
		$honc{$tto}=$nchil if $tto;
		$tto=$tt[0];
		$to=$t;
		$tio=$ti;
		$tso=$ts;
		$teo=$te;
		$tdo=$td;
		$nchil=0;
		$withinregion=0; # set to 1 if you want gene to absorbed (and then loseq
	} elsif ( ($t eq 'gene') & !($withinregion) ) { # importatn enough to have its own branch conditional
		$hofn{$tt[0]}=[$t, $ti, $td, $ts, $te];
		$hopf{$t}++; # to count the primary features.
		$honc{$tto}=$nchil if $tto;
		$tto=$tt[0];
		$to=$t;
		$tio=$ti;
		$tso=$ts;
		$teo=$te;
		$tdo=$td;
		$nchil=0;
		$withinregion=0;
	} elsif ( ($to ne $t) & ($tso < $ts) & ($teo < $te) ) { # chromosome is the same but feature is new and range is different: ignore identical ranges with differnet names.
		$hofn{$tt[0]}=[$t, $ti, $td, $ts, $te];
		$hopf{$t}++; # to count the primary features.
		$honc{$tto}=$nchil if $tto;
		$tto=$tt[0];
		$to=$t;
		$tio=$ti;
		$tso=$ts;
		$teo=$te;
		$tdo=$td;
		$nchil=0;
	} elsif ( ($to ne $t) & ($tso <= $ts) & ($teo >= $te) ) {
		push ( @{$hochil{$tto}}, $t);
		$nchil++;
	}
}
print "last primfeat $t\n";

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

# use subroutine to print hash
printohchrom(\%hofn, \%honc, \%hochil);
