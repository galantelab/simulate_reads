#!/usr/bin/env perl

use strict;
use warnings;

my @list = (1, 2, 3, 10, 20, 30, 100, 200, 300);
my $range = bs(\@list, 100, 300);
print @$range ? "[@$range]\n" : "not found!\n";

sub bs {
	my ($list, $down, $up) = @_;
	my $max_ind = $#$list;
	my $left_ind = _bs($list, 0, $max_ind, $down);
	my @range;
	for (my $i = $left_ind; $i < @$list; $i++) {
		push @range => $list->[$i] if $list->[$i] <= $up;
	}
	return \@range;
}

sub _bs {
	my ($list, $start, $end, $down) = @_;
	return $start if $start > $end;

	my $ind = int(($start + $end) / 2);
	my $guess = $list->[$ind];
	print "[$start-$end] [$ind] $guess\n";
	die "ponga" if not defined $guess;

	if ($guess < $down) {
		return _bs($list, $ind + 1, $end, $down);
	} elsif ($guess > $down) {
		return _bs($list, $start, $ind - 1, $down);
	} else {
		return $ind;
	}
}
#sub _bs {
#	my ($list, $start, $end, $down, $up) = @_;
#	return if $start > $end;
#
#	my $ind = int(($start + $end) / 2);
#	my $guess = $list->[$ind];
#	die "Not defined entry at range [$start-$end]" if not defined $guess;
#
#	if ($guess < $down) {
#		return _bs($list, $ind + 1, $end, $down, $up);
#	} elsif ($guess > $up) {
#		return _bs($list, $start, $ind - 1, $down, $up);
#	} else {
#		return $guess;
#	}
#}
