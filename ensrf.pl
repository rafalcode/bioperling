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

# The access methods are farily strict
my @SQA; # sequence hash
my $cou=0;
my $couex=0;
my($T, @Q);
while ( my $sq = $in->next_seq() ) {
    if($sq->desc =~ /cds/m) {
        $T=$sq;
    } elsif ($sq->desc =~ /exon/m) {
        # we're just going to take the first two exons.
        if( ($couex==1) | ($couex==2) ) {
            push @Q, $sq;
        }
        $couex++;
    }
    $cou++;
}

# print $T->seq."\n";
print "exon count=$couex\n";
# I often can't get this to work properly
# print "We have ". $#Q+1 ." exons to check against target\n";
# you're better off with printf, something you know quite well.
printf "%s %d %s\n", "We have", $#Q+1, "exons to check against target";
# also note that I always fall forgetting thet $# is only index of last member, not size.
# print "Mean length ",$tot/$cou," Median ",$SQA[$cou/2],"\n";


my $of="uit.wa";
my $of2="uit.ma";
my $f = Bio::Factory::EMBOSS->new();

# my $water = $f->program('water');
# $water->run({'-sequencea' => $T, '-seqall' => \@Q, '-gapopen' => '10.0', '-gapextend' => '0.5', '-outfile' => $of});

# This is giving
#  sh: 1: Syntax error: "(" unexpected
#  whihc measn the system call - or whatever - to EMBOSS, isn't goign so well out on the system.

my $water = $f->program('water');
# my $wateri = $f->program_info('water');
# A string is returned from this.

# my $pma = $f->program('matcher');
# my $pmai = $f->program_info('matcher');
# $pma->run({-sequencea => $T, -sequenceb => $[0], -gapopen   => '10.0', -gapextend => '0.5', -outfile => $f});
# $pma->run({'-sequencea' => $T, '-sequenceb' => $Q[0], '-outfile' => $of});
# $pma->run({ -asequence => $T, -bsequence => $Q[0], -aformat => "pair", -alternatives => 1, -outfile => $of2});

# using cpan docs ...
my $woutfile = 'out.water';
# $water->run({-asequence => $T, -seqall => \@Q, -gapopen => '10.0', -gapextend => '0.5', -outfile => $wateroutfile});
# # I thnk the seqall option is dead.
# can easily get a "Can't call method "isa" on unblessed reference at /usr/share/perl5/Bio/Tools/Run/EMBOSSApplication.pm line 168." error here.
$water->run({-asequence => $T, -bsequence => \@Q, -gapopen => '10.0', -gapextend => '0.5', -outfile => $woutfile});
# and yes, you can give multiseq for bsequence.
#
# OK, let's analyze that alignment
#
use Bio::AlignIO;
my $alignio_fmt = "emboss";
my $inaln = new Bio::AlignIO(-format => $alignio_fmt, -file => $woutfile);
while ( my $aln = $inaln->next_aln ) {
    print $aln->score."\n";
}

# I want the start and end coords for each of the sequences.
my $TXTENC = ":encoding(UTF-8)"; # swanky
my $FH = undef; # file handle  .. but leave undefined: a declaration of intentions.
open($FH, "< $TXTENC", $woutfile);
my $AIN=0; # marker ina an alignment
my $ACOU=0; # number of alignments
my (@AT, @A, @B);
my $LAT; 
my($E, $S1S, $S2S, $S1E, $S2E);

while(<$FH>) {
    if( @AT = ($_ =~ /^[^# ]\S+\s+(\d+)\D+(\d+)/)) { 
        if($AIN==0) {
            $AIN=1;
            $ACOU++;
            $S1S=
            push @A, @AT;
        } elsif ($AIN==1) {
            push @B, @AT;
            $AIN==2;
        } elsif($AIN==2) {
            push @A, @AT;
            undef @B;
            $AIN=1;
        }
    } elsif( /^#/ && $AIN==2) {
        $S2E=$A[1];
        $AIN=0;
    }
}
print "Num of alignments = $ACOU\n";
# print "$S1S->$S1E/$S2S->$S2E\n"; 
close($FH);
