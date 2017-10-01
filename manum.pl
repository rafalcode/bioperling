#!/usr/bin/env perl
# gets output of matcher
use strict;
use warnings;

# Argument accounting
if($#ARGV != 0) {
    print "Usage error: this script manipulates the default output of matcher, so it requires a file with that output.\n";
    die;
}

#
# Preparing first arg which is a file name
my $FN = $ARGV[0];
my @FNR = ($ARGV[0] =~ /\/*([^\.]+)/m); # filename root
my $TXTENC = ":encoding(UTF-8)"; # swanky
my $FH = undef; # file handle  .. but leave undefined: a declaration of intentions.
# but don't open just yet ...

# declare global vars.
my $I;
my $NS0; # first space size
my $LS; # last space
my $NS2; # last space size
my $RN0; # real number first
my $RN2; # real number last
my $ND0; # num digits of first number
my $ND2; # num digits of last number
my $MARGIN=7; # matcher output only allows 6 characters for the fast ID plus a space: this seems to be hardcoded.
my @A; # our holding array.
my @AT; # temp array
my %H;

my $SEEN=0;
# OK let's open.
open($FH, "< $TXTENC", $FN);
while(<$FH>) {
    next if /^#/m;
    if( @AT = ($_ =~ /( +)(\d+)/g) && $SEEN==0) {
        $NS0=length $AT[0];
        $ND0=length $AT[1];
        $ND2=length $AT[$#AT];
        /[^ ]( +)$/m; # grab final, solitary space
        $NS2=length $1;
        $RN0=$AT[1] - $NS0 - $ND0 + $MARGIN +1;
        $RN2=$AT[$#AT] + $NS2;
        print "LastATZ $AT[$#AT]; NS2=$NS2\n";
        print "Line start idx: $RN0; Line end idx: $RN2\n";
        push @A, $RN0;
        push @A, $RN2;
        $SEEN=1;
    } elsif( /^\S+\s+(\S+)/m && $SEEN) {
        $SQ=$1;
        print "$SQ\n";
    } elsif( @AT = ($_ =~ /( +)(\d+)/g) && $SEEN) {
        $NS0=length $AT[0];
        $ND0=length $AT[1];
        $ND2=length $AT[$#AT];
        /[^ ]( +)$/m; # grab final, solitary space
        $NS2=length $1;
        $RN0=$AT[1] - $NS0 - $ND0 + $MARGIN +1;
        $RN2=$AT[$#AT] + $NS2;
        print "LastATZ $AT[$#AT]; NS2=$NS2\n";
        print "Line start idx: $RN0; Line end idx: $RN2\n";
        push @A, $RN0;
        push @A, $RN2;
        $SEEN=0;
    }
}
close($FH);

# OK so every is in our array now. It's a good idea to check
# s1 e1 s2 e2 s1 e1 s2 e2
# that e1 = s1+1 except for first and last. Here's how you do it:
for($I=1; $I<$#A-2; $I+=2) {
    if($A[$I+3] != $A[$I]+1) {
        print "$A[$I+3] != $A[$I]+1!\n";
        print "Error: unfortunately, the integrity of the collected array has suffered a disturbance\n.";
        print "This probably due to gaps ... they have character spaces, but the indices do not accumulate these.\n";
        print "Please check this script and the target matcher output and try again\n";
    }
}

# If everything is OK all we need is the start and end (+1) of seqA and seqB
printf("%s SeqA %d %d\n", $FNR[0], $A[0], $A[$#A-2]+1);
printf("%s SeqB %d %d\n", $FNR[0], $A[2], $A[$#A]+1);
# print "SeqB $A[1] ". $A[$#A]+1 ."\n";
