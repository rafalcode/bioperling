#!/usr/bin/env perl
use warnings;
use strict;

my %h;

my $el0='red';
my $el1=undef;
my $el2='green';

$h{$el0}=5 if defined $el0;
if(defined $el1) {
	$h{$el1}=6;
}
$h{$el2}=7 if defined $el2;


# print " $h{$el0} $h{$el1} $h{$el2} \n";
print " $h{$el0} $h{$el2} \n";
