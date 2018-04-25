#!/usr/bin/env perl 

use strict;
use warnings;
use IntervalTree;
use Data::Dumper::Simple;

my @list = (
	{ down => 1, up => 5,   }, 
	{ down => 8, up => 15,  }, 
	{ down => 2, up => 2,   }, 
	{ down => 10, up => 50, }, 
	{ down => 12, up => 18, }, 
	{ down => 25, up => 25, }, 
	{ down => 30, up => 50, }, 
	{ down => 3, up => 5,   }, 
	{ down => 5, up => 200, }, 
	{ down => 189, up => 1000}
);

my $tree = IntervalTree->new;

for (@list) {
	$tree->insert($_->{down}, $_->{up}, $_);
}

my $code = sub {
	my $node = shift;
	printf "[%d - %d] max = %d height = %d\n" => $node->low, $node->high, $node->max, $node->height;
};


my $node = $tree->left_near(1);
print Dumper($node->data);
#print "Preorder traversal\n";
#$tree->preorder($code);
#print "Inorder traversal\n";
#$tree->inorder($code);
#
#print "Remove [189 1000]\n";
#my $node = $tree->delete(189, 1000);
#print "Preorder traversal\n";
#$tree->preorder($code);
#print "Inorder traversal\n";
#$tree->inorder($code);
#
#print "Remove [1 5]\n";
#my $node = $tree->delete(1, 5);
#print "Preorder traversal\n";
#$tree->preorder($code);
#print "Inorder traversal\n";
#$tree->inorder($code);
#
#print "Remove [10 50]\n";
#my $node = $tree->delete(10, 50);
#print "Preorder traversal\n";
#$tree->preorder($code);
#print "Inorder traversal\n";
#$tree->inorder($code);
#
#print "Remove [5 200]\n";
#my $node = $tree->delete(5, 200);
#print "Preorder traversal\n";
#$tree->preorder($code);
#print "Inorder traversal\n";
#$tree->inorder($code);
#print "Postorder traversal\n";
#$tree->postorder($code);
#
#my $nodes = $tree->search(1, 1000);
#printf "[%d - %d]\n", $_->low, $_->high for @$nodes;
