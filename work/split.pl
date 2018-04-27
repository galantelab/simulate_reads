#!/usr/bin/env perl 

use strict;
use warnings;
use IntervalTree;

my $seq = "AAATTTCGCGGCTGGATATAGAGGCGGATTAGAGAGAGATCGGATATAGGGAGAGAGATATGCGGATTAGAGAGGCGCGCTTAGAGAG";

my $pos = 10;
my $subseq = "GCTGGATATA";
my $change = "TTT";
my $end = length($subseq) - 1;

my $tree = IntervalTree->new();
$tree->insert(0, length($seq) - 1, \$seq);

my $nodes = $tree->search($pos, $pos + $end);
my $merge = join "" => map { ${ $_->data } }  @$nodes;

my $split1 = substr $merge, 0, $pos;
my $split2 = substr $merge, $pos, $end + 1;
my $split3 = substr $merge, $end + 1, length($merge) - 1;

$tree->inorder(\&print_tree);
print "s1 = '$split1', s2 = '$split2', s3 = '$split3', change = '$change'\n";

$tree->delete($_->low, $_->high) for @$nodes;
$tree->insert(0, $pos - 1, \$split1); 
$tree->insert($pos, $pos + $end, \$split2);
$tree->insert($pos + $end + 1, $pos + $end + length($merge), \$split3);
$tree->insert($pos, $pos + length($change) - 1, \$change);
$tree->inorder(\&print_tree);


sub print_tree {
	my $node = shift;
	printf "[%d - %d] seq = %s\n" => $node->low, $node->high, ${ $node->data };
}
