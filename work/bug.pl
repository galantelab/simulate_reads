#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper::Simple;
use IntervalTree;

my @data = (
	{ low => 0,   up => 2,   data => 'AAC HE' },
	{ low => 4,   up => 12,  data => '- AAAAAAAAA' },
	{ low => 5,   up => 5,   data => 'T A' },
	{ low => 979, up => 979, data => 'A T' },
	{ low => 981, up => 983, data => 'CTT -' },
	{ low => 62,  up => 70,  data => 'GACCGCCCA TTT' }
);

my $tree = IntervalTree->new;
$tree->insert($_->{low}, $_->{up}, \$_->{data}) for @data;

my $code =  sub {
	my $node = shift;
	printf "[%d - %d] max = %d height = %d\n", $node->low, $node->high, $node->max, $node->height;
};

print "TREE: \n";
$tree->preorder($code);

print "SEARCH [0 - 99]\n";
my $nodes = $tree->search(0, 99);
printf "[%d - %d] max = %d\n", $_->low, $_->high, $_->max for @$nodes;
