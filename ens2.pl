#!/usr/bin/env perl
use warnings;
use strict;

# Argument accounting ... say how many the scrit should have.
if(scalar @ARGV != 1) {
    print "Sorry, this script requires exactly one argument\n";
    die;
}

use Bio::Seq;
use Bio::AlignIO;
use Bio::Factory::EMBOSS;
use Bio::SeqIO;

my $in = Bio::SeqIO->new(-file => $ARGV[0], '-format' => 'Fasta');
my @SQA; # sequence hash
my $cou=0;
my $couex=0;
my($T, @Q);
while ( my $sq = $in->next_seq() ) {
    if($sq->desc =~ /cds/m) {
        $T=$sq;
    } elsif ($sq->desc =~ /exon/m) {
        if( ($couex==1) | ($couex==2) ) { # just experimenting on two right now.
            push @Q, $sq;
        }
        $couex++;
    }
    $cou++;
}

# my $factory = new Bio::Factory::EMBOSS;
my $f = Bio::Factory::EMBOSS->new();

# our sequences as strings:
my $seq1 = 'MAVNPELAPFTLSRGIPSFDDQALSTIIQLQDCIQQAIQQLNYSTAEFLAELLYAECSILDKSSVYWSDAVYLYALSLFLNKSYHTAFQISKEFKEYHLGIAYIFGRCALQLSQGVNEAILTLLSIINVFSSNSSNTRINMVLNSNLVHIPDLATLNCLLGNLYMKLDHSKEGAFYHSEALAINPYLWESYEAICKMRATVDLKRVFFDIAGKKSNSHNNNAASSFPSTSLSHFEPRSQPSLYSKTNKNGNNNINNNVNTLFQSSNSPPSTSASSFSSIQHFSRSQQQQANTSIRTCQNKNTQTPKNPAINSKTSSALPNNISMNLVSPSSKQPTISSLAKVYNRNKLLTTPPSKLLNNDRNHQNNNNNNNNNNNNNNNNNNNNNNNNIINKTTFKTPRNLYSSTGRLTTSKKNPRSLIISNSILTSDYQITLPEIMYNFALILRSSSQYNSFKAIRLFESQIPSHIKDTMPWCLVQLGKLHFEIINYDMSLKYFNRLKDLQPARVKDMEIFSTLLWHLHDKVKSSNLANGLMDTMPNKPETWCCIGNLLSLQKDHDAAIKAFEKATQLDPNFAYAYTLQGHEHSSNDSSDSAKTCYRKALACDPQHYNAYYGLGTSAMKLGQYEEALLYFEKARSINPVNVVLICCCGGSLEKLGYKEKALQYYELACHLQPTSSLSKYKMGQLLYSMTRYNVALQTFEELVKLVPDDATAHYLLGQTYRIVGRKKDAIKELTVAMNLDPKGNQVIIDELQKCHMQE'; 
my $seq2 = 'CLIFXRLLLIQMIHPQARRAFTFLQQQEPYRIQSMEQLSTLLWHLADLPALSHLSQSLISISRSSPQAWIAVGNCFSLQKDHDEAMRCFRRATQVDEGCAYAWTLCGYEAVEMEEYERAMAFYRTAIRTDARHYNAWYVLFFFFFFFFVPGDIDSXPKKGMEWGXFISKRIDRGMRSIILKEPSKSIQLIPFFYVALVWXVGVSSYPLETMTNIDFPKKKKALEKSNDVVQALHFYERASKYAPTSAMVQFKRIRALVALQRYDEAISALVPLTHSAPDEANVFFLLGKCLLKKERRQEATMAFTNARELEPK';

# As the received advice goes, get seq obects from your sequence strings.
my $sqo1=Bio::Seq->new(-id => "seq1", -seq => $seq1);
my $sqo2=Bio::Seq->new(-id => "seq2", -seq => $seq2);

my $prog = $f->program('matcher');
my $outfile='uit3.wa';

# $prog->run({ -asequence => Bio::Seq->new(-id => "seq1", -seq => $seq1), -bsequence => Bio::Seq->new(-id => "seq2", -seq => $seq2), -aformat => "pair", -alternatives => 2, -outfile => $outfile});
$prog->run({ -asequence => $sqo1, -bsequence => $sqo2, -aformat => "markx10", -alternatives => 2, -outfile => $outfile});

# so the matcher output has been serialized to disk
my $alignio_fmt = "emboss";
my $align_io = new Bio::AlignIO(-format => $alignio_fmt, -file => $outfile);
