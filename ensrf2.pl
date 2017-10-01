#!/usr/bin/env perl
# trying out this script does what?
use strict;
use warnings;
# Argument accounting ... say how many the scrit should have.
if(scalar @ARGV != 1) {
    print "Sorry, this script requires exactly one argument\n";
    die;
}

# The SeqIO individual sequence object
# $seqobj->id;       # prob only the first part of the ID ... i.e. not the whole string
# $seqobj->display_id();       # the human read-able id of the sequence
# $seqobj->seq();              # string of sequence
# $seqobj->subseq(5,10);       # part of the sequence as a string
# $seqobj->accession_number(); # when there, the accession number
# $seqobj->alphabet();         # one of 'dna','rna','protein'
# $seqobj->primary_id();       # a unique id for this sequence irregardless of its display_id or accession number
# $seqobj->desc();             # a description of the sequence
# $seqobj->trunc(5,10);  # truncation from 5 to 10 as new object
# $seqobj->revcom;       # reverse complements sequence
# $seqobj->translate;    # translation of the sequence
#
# get an EMBOSS factory
use Bio::Factory::EMBOSS;
use Bio::SeqIO;
my $in = Bio::SeqIO->new(-file => $ARGV[0], '-format' => 'Fasta');

my @SQA; # sequence hash
my $cou=0;
my $couex=0;
my ($I, $J);
my($T, @Q); # target and query. Query can be multisequence, while we expect target to be a single sequence, often much bigger than the individual query sequences.

# We're expecting ensemb cds plus exon fasta files, so:
while ( my $sq = $in->next_seq() ) {
    if($sq->desc =~ /cds/m) {
        $T=$sq;
    } elsif ($sq->desc =~ /exon/m) {
        # we're just going to take the first two exons.
        # if( ($couex==1) | ($couex==2) ) {
            push @Q, $sq;
            # }
        $couex++;
    }
    $cou++;
}


# Inform user of alignment intentions:

printf "Length of target sequence = %d\n", $T->length();
printf "%s %d %s\n", "No. of query (chosen exons) sequences = ", $#Q+1, ". Sizes:";
for $J (@Q) {
    print $J->length(). " ";
}
print "\n";

printf "%s %d %s\n", "We have", $#Q+1, "exons to check against target";

my $f = Bio::Factory::EMBOSS->new();

my $water = $f->program('water');
my $woutfile = 'out.water';

$water->run({-asequence => $T, -bsequence => \@Q, -gapopen => '10.0', -gapextend => '0.5', -outfile => $woutfile});

# OK, let's analyze that alignment
use Bio::AlignIO;
my $alignio_fmt = "emboss";
my $inaln = new Bio::AlignIO(-format => $alignio_fmt, -file => $woutfile);

# I want the start and end coords for each of the sequences.
# # by default bioperl doesn't give this, though it is in the output file.
my $TXTENC = ":encoding(UTF-8)"; # swanky
my $FH = undef; # file handle ... but set undefined: a declaration of intentions.
open($FH, "< $TXTENC", $woutfile);
my $AIN=0; # marker ina an alignment
my $ACOU=0; # number of alignments
my (@AT, @A, @B); # temporary array. per sequence accumulation array, overall start/stop array.

my @ISG;
while(<$FH>) {
# Identity:     400/1151 (34.8%)
    if( @AT = ($_ =~ /^# Identity:\s+(\d+)/) ) { 
        $ISG[3*$ACOU]=$AT[0];
    } elsif( @AT = ($_ =~ /^# Similarity:\s+(\d+)/) ) { 
        $ISG[3*$ACOU+1]=$AT[0];
    } elsif( @AT = ($_ =~ /^# Gaps:\s+(\d+)/) ) { 
        $ISG[3*$ACOU+2]=$AT[0];
    } elsif( @AT = ($_ =~ /^[^# ]\S+\s+(\d+)\D+(\d+)/)) { 
        if($AIN==0) {
            $AIN=1;
            $ACOU++;
        }
        push @A, @AT;
    } elsif( $AIN && /^#/) {
        $J=4*($ACOU-1);
        $B[$J]=$A[0]; # start index of 1st sequence
        $B[$J+1]=$A[$#A-2]; # end index of 1st sequence
        $B[$J+2]=$A[2]; # start index of 2nd sequence
        $B[$J+3]=$A[$#A]; # end index of 2nd sequence
        $AIN=0;
        undef @A;
    }
}
print "Num of alignments = $ACOU\n";
# for($J=0; $J<=$#B; $J+=2) {
#     print "$B[$J]:$B[$J+1] ";
# }
$ACOU=0;
my($PET, $PEQ);
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "LEN", "SCORE", "IDEN", "IPT", "SIM", "GAPS", "GPT", "TSC", "TEC", "PET", "QSC", "QEC", "PEQ";
# TEC target end Coordinate, IPT iden percentage
while ( my $aln = $inaln->next_aln ) {
    $I=3*$ACOU;
    $J=4*$ACOU;
    $PET=100*($B[$J+1] - $B[$J] +1)/$T->length();
    $PEQ=100*($B[$J+3] - $B[$J+2] +1)/$Q[$ACOU]->length();;
    printf "%d\t%4.1f\t%d\t%3.1f\t%d\t%d\t%3.1f\t%d\t%d\t%3.1f\t%d\t%d\t%3.1f\n", $aln->length(), $aln->score(), $ISG[$I], 100*$ISG[$I]/$aln->length(), $ISG[$I+1], $ISG[$I+2], 100*$ISG[$I+2]/$aln->length(), $B[$J], $B[$J+1], $PET, $B[$J+2], $B[$J+3], $PEQ;
    $ACOU++;
}
close($FH);
