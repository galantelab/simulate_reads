#!/usr/bin/env perl

use strict;
use warnings;
use autodie;
use Data::Dumper::Simple;
use Scalar::Util 'looks_like_number';
use Devel::Size qw(size total_size);
use IntervalTree;

use constant {
	PLOIDY    => 2,
	READ_SIZE => 100
};

sub subseq {
	my ($seq_ref, $seq_len, $slice_len, $pos) = @_;
	my $read = substr $$seq_ref, $pos, $slice_len;
	return \$read;
}

sub subseq_rand {
	my ($seq_ref, $seq_len, $slice_len) = @_;
	my $usable_len = $seq_len - $slice_len;
	my $pos = int(rand($usable_len + 1));
	my $read = substr $$seq_ref, $pos, $slice_len;
	return (\$read, $pos);
}

sub index_snp {
	my $snp = shift;

	open my $fh, "<" => $snp;

	my %indexed_snp;
	my $line = 0;

	LINE:
	while (<$fh>) {
		$line++;
		chomp;
		next if /^\s*$/;
		my @fields = split;
		die "Not found all fields (SEQID, POSITION, REFERENCE, OBSERVED, PLOIDY) into file '$snp' at line $line\n"
			unless scalar @fields == 5;

		die "Second column, position, does not seem to be a number into file '$snp' at line $line\n"
			unless looks_like_number($fields[1]);

		die "Second column, position, has a value lesser or equal to zero into file '$snp' at line $line\n"
			if $fields[1] <= 0;

		$fields[1] = int($fields[1]);

		die "Third column, reference, does not seem to be a valid entry: '$fields[2]' into file '$snp' at line $line\n"
			unless $fields[2] =~ /^(\w+|-)$/;

		die "Fourth column, alteration, does not seem to be a valid entry: '$fields[3]' into file '$snp' at line $line\n"
			unless $fields[3] =~ /^(\w+|-)$/;

		die "Fifth column, ploidy, has an invalid entry: '$fields[4]' into file '$snp' at line $line. Valid ones are 'HE' or 'HO'\n"
			unless $fields[4] =~ /^(HE|HO)$/;

		if ($fields[2] eq $fields[3]) {
			warn "There is an alteration equal to the reference at '$snp' line $line. I will ignore it\n";
			next;
		}

		my %variation = (
			ref  => $fields[2],
			alt  => $fields[3],
			plo  => $fields[4],
			pos  => $fields[1] - 1
		);
		
		if (not defined $indexed_snp{$fields[0]}) {
			$indexed_snp{$fields[0]} = IntervalTree->new;
		}

		my $tree = $indexed_snp{$fields[0]};

		# Sequence inside perl begins at 0
		my $position = $fields[1] - 1;

		# Compare the alterations and reference to guess the max variation on sequence
		my $size_of_variation = ( sort { $b <=> $a } map { length } $fields[3], $fields[2] )[0];

		my $low = $position;
		my $high = $position + $size_of_variation - 1;

		my $nodes = $tree->search($low, $high);

		# My rules:
		# The biggest structural variation gains precedence.
		# If occurs an overlapping, I search to the biggest variation among
		# the saved alterations and compare it with the actual entry:
		#    *** Remove all overlapping variations if actual entry is bigger;
		#    *** Skip actual entry if it is lower than the biggest variation
		#    *** Insertionw can be before any alterations
		my @to_remove;

		NODE:
		for my $node (@$nodes) {
			my $data = $node->data;

			# Insertion after insertion
			if ($data->{ref} eq '-' && $fields[2] eq '-' && $fields[1] != $data->{pos}) {
				next NODE;

			# Alteration after insertion
			} elsif ($data->{ref} eq '-' && $fields[2] ne '-' && $fields[1] > $data->{pos}) {
				next NODE;

			# In this case, it gains the biggest one
			} else {
				my $size_of_variation_saved = $node->high - $node->low + 1;

				if ($size_of_variation_saved >= $size_of_variation) {
					warn sprintf "Early alteration [%s %d %s %s %s] masks [%s %d %s %s %s] at '%s' line %d\n"
						=> $fields[0], $data->{pos}+1, $data->{ref}, $data->{alt}, $data->{plo}, $fields[0], $position+1,
						$fields[2], $fields[3], $fields[4], $snp, $line;

					next LINE;
				} else {
					warn sprintf "Alteration [%s %d %s %s %s] masks early declaration [%s %d %s %s %s] at '%s' line %d\n"
						=> $fields[0], $position+1, $fields[2], $fields[3], $fields[4], $fields[0], $data->{pos}+1, $data->{ref},
						$data->{alt}, $data->{plo}, $snp, $line;

					push @to_remove => $node;
				}
			}
		}

		# Remove the overlapping node
		$tree->delete($_->low, $_->high) for @to_remove;

		my %variation = (
			ref  => $fields[2],
			alt  => $fields[3],
			plo  => $fields[4],
			pos  => $position
		);

		$tree->insert($low, $high, \%variation);
	}

	close $fh;

	return \%indexed_snp;
}

sub _index_fasta_tree {
	my ($fasta_index, $snp_index) = @_;
	my %fasta_tree;
	
	for my $seq_id (keys %$fasta_index) {
		my $seq_data = delete $fasta_index->{$seq_id};
		my $seq_ref = \$seq_data->{seq};
		my $size = $seq_data->{size};
		print "$seq_id: seq = $seq size = $size\n";

		my $ftree = IntervalTree->new;

		my %data = (
			seq  => $seq_ref,
			pos  => 0,
			size => $size,
			ref  => 1
		);

		$ftree->insert(0, $size - 1, \%data);
		my $stree = $snp_index->{$seq_id};

		if (@$stree) {
			_apply_str_var($ftree, $stree);
		}

		$fasta_tree{$seq_id} = $ftree;
	}

	return \%fasta_tree;
}

sub _apply_str_var {
	my ($ftree, $stree) = @_;

	my @stree_nodes;

	my $catcher = sub { push @stree_nodes => shift };
	$stree->inorder($catcher);

	for my $stree_node (@stree_nodes) {
		my $ftree_nodes = $ftree->search($stree_node->low, $stree_node->high);

		if (not @$ftree_nodes) {
			break;
		}

		my $snode_data = $stree_node->data;

		for my $ftree_node (@$ftree_nodes) {
			my $fnode_data = $ftree_node->data;

			my $seq_ref = $fnode_data->{seq};
			my $seq_size = $fnode_data->{size};

			my $split1 = substr $$seq, 0, $pos;
			my $split2 = substr $merge, $pos, $end + 1;
			my $split3 = substr $merge, $end + 1, length($merge) - 1;
		}
	}
}

sub _index_fasta {
	my $fasta = shift;

	open my $fh, "<" => $fasta;

	# indexed_genome = ID => (seq, len)
	my %indexed_fasta;

	# >ID|PID as in gencode transcripts
	my $id;

	while (<$fh>) {
		chomp;
		next if /^;/;
		if (/^>/) {
			my @fields = split /\|/;
			$id = $fields[0];
			$id =~ s/^>//;
			$id =~ s/^\s+|\s+$//g;

			# It is necessary to catch gene -> transcript relation
			# # TODO: Make a hash tarit for indexed fasta
			if (defined $fields[1]) {
				my $pid = $fields[1];
				$pid =~ s/^\s+|\s+$//g;
			}
		} else {
			die "Error reading fasta file '$fasta': Not defined id"
				unless defined $id;
			$indexed_fasta{$id}{seq} .= $_;
		}
	}

	for (keys %indexed_fasta) {
		my $size = length $indexed_fasta{$_}{seq};
		$indexed_fasta{$_}{size} = $size;
		$indexed_fasta{$_}{weight} = $size;
	}

	unless (%indexed_fasta) {
		die "Error parsing '$fasta'. Maybe the file is empty\n";
	}

	$fh->close
		or die "Cannot close file $fasta: $!\n";

	return \%indexed_fasta;
}
