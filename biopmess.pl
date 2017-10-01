#!/usr/bin/env perl
# this script does what?
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
while ( my $sq = $in->next_seq() ) {
    #    print $sq->id."\n"; # Don't put the dstruc within quotation, must be outside them
    push @SQA, $sq; # an array of Sq objects; parts "->id", "->length"
    # print $sq->display_id."\n";
}


# Do you want to sort of lengths?
# @SQA = sort { $a->length <=> $b->length } @SQA;

my $tot = 0;
my $cou = 0;
my ($T, @Q); # the target and query sequence
foreach my $sq (@SQA) {
   $tot += $sq->length; 
   print "Seq num $cou is ". $sq->desc."\n" if $sq->desc =~ /cds/m;
   if($sq->desc =~ /cds/m) {
       $T=$sq->seq;
   } elsif ($sq->desc =~ /exon/m) {
       push @Q, $sq->seq;
   }
   $cou++;
}

print $T."\n";
print "We have $#Q exons to check against target\n";
# print "Mean length ",$tot/$cou," Median ",$SQA[$cou/2],"\n";

my $f = Bio::Factory::EMBOSS -> new();
# get an EMBOSS application  object from the factory
# $water = $f->program('water');
# my $water = $f->program_info('water');
# A string is returned from this.
my $progmatcher = $f->program_info('matcher');

print "what? $progmatcher\n";
