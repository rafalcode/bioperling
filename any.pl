#!/usr/bin/env perl
use warnings;
use strict;

use Bio::Seq;
use Bio::AlignIO;
use Bio::Factory::EMBOSS;

# my $factory = new Bio::Factory::EMBOSS;
my $f = Bio::Factory::EMBOSS->new();

# our sequences as strings:
my $seq1 = 'TGGAGGATGACAATGGAGGAGATGAAGAATGAAGCGGAGACAAACTCCATGGTATCCATG';
my $seq2 = 'ACCAACCTCATACTGCAGACCTTCAAAACCGTGGCCTAG';

# As the received advice goes, get seq obects from your sequence strings.
my $sqo1=Bio::Seq->new(-id => "seq1", -seq => $seq1);
my $sqo1r=$sqo1->revcom();
print $sqo1->seq ."\n". $sqo1r->seq ."\n";
my $sqo2=Bio::Seq->new(-id => "seq2", -seq => $seq2);

my $prog = $f->program('matcher');
my $out1='uit1.ma';
my $out2='uit2.ma';

# $prog->run({ -asequence => Bio::Seq->new(-id => "seq1", -seq => $seq1), -bsequence => Bio::Seq->new(-id => "seq2", -seq => $seq2), -aformat => "pair", -alternatives => 2, -outfile => $outfile});
$prog->run({ -asequence => $sqo1, -bsequence => $sqo2, -aformat => "markx10", -outfile => $out1});
$prog->run({ -asequence => $sqo1r, -bsequence => $sqo2, -aformat => "markx10", -outfile => $out2});
