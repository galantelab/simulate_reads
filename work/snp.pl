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

my $fasta = shift;
my $snp = shift;

my $fasta_index = _index_fasta($fasta);
my $snp_index = index_snp($snp);

my $code = sub {
	my $node = shift;
	printf "[%d - %d] max = %d\n" => $node->low, $node->high, $node->max;
};

print "Finito!\n";

for my $seq_id (sort keys %$snp_index) {
	my $tree = $snp_index->{$seq_id};
	print "chr = $seq_id\n";
	$tree->inorder($code);
}
#print Dumper($snp_index);
#my $pos = 0;
#my $read_r = subseq(\$fasta_index->{chr1}{seq}, $fasta_index->{chr1}{size}, READ_SIZE, $pos);
#print "read = $$read_r\nlength = ", length $$read_r, "\n";
#my $variations = bs($snp_index->{chr1}, $pos, $pos + READ_SIZE - 1);
#print Dumper($variations);
#
#for my $variation (@$variations) {
#	if ($variation->{plo} eq 'HE' && int(rand(2))) {
#		next;
#	}
#
#	# TODO: FIx HERE ->
#	# There is a problem with strctural variations. They shift
#	# the sequence so that it cannot find the next changes by position
#	# I think to solve, I need to make a tracker for the shift and when 
#	# testing and substr the sequence, I correct  the position.
#	# If the user request a big deletion, the read will be truncated and 
#	# to solve the algorithm catch the rest needed from the reference fasta.
#	# But if really it is a big deletion, or a deletion at the end of the read,
#	# It cannot be sovlved and the deletion will be "truncated".
#	my $position = $variation->{pos} - $pos;
#	my $reference = $variation->{ref};
#	my $alterations = $variation->{alt};
#	my $alteration = $alterations->[int(rand(scalar @$alterations))];
#
#	# I need to test, again, if the reference pattern
#	# still is there, because users may pass variations
#	# too close that they overlap. In short, it causes
#	# a race condition-like when the first one to occurs
#	# gain the precedence.
#	if ($reference eq '-') {
#		# Insertion
#		substr($$read_r, $position, 0) = $alteration;
#	} elsif ($alteration eq '-') {
#		# Deletion
#		if (substr($$read_r, $position, length($reference)) eq $reference) {
#			substr($$read_r, $position, length($reference)) = "";
#		}
#	} else {
#		# Change
#		if (substr($$read_r, $position, length($reference)) eq $reference) {
#			substr($$read_r, $position, length($reference)) = $alteration;
#		}
#	}
#}
#
#my $new_read_size = length $$read_r;
#
#if ($new_read_size < READ_SIZE) {
#	my $missing = READ_SIZE - $new_read_size;
#	my $beacon = $pos + READ_SIZE;
#	substr($$read_r, $new_read_size, 0) = substr($fasta_index->{chr1}{seq}, $beacon, $missing);
#	$new_read_size = length($$read_r);
#	if ($new_read_size < READ_SIZE) {
#		# May occur that the end of the fasta sequence suffers an deletion,
#		# so I cannot guess how to fill the read size. In that case,
#		# I fill the string with 'N's
#		$missing = READ_SIZE - $new_read_size;
#		substr($$read_r, $new_read_size, 0) = 'N' x $missing;
#	}
#} elsif ($new_read_size > READ_SIZE) {
#	# Just truncate!
#	$$read_r = substr($$read_r, 0, READ_SIZE);
#}
#
#print "read = $$read_r\nlength = ", length $$read_r, "\n";

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
#	my $ploidy = PLOIDY;
#	my $read_size = READ_SIZE;

	open my $fh, "<" => $snp;

	my %indexed_snp;
	my $line = 0;
	# chr pos ref @obs he

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
			unless $fields[3] =~ /^(\w+,*\w*|-)$/;

		die "Fifth column, ploidy, has an invalid entry: '$fields[4]' into file '$snp' at line $line. Valid ones are 'HE' or 'HO'\n"
			unless $fields[4] =~ /^(HE|HO)$/;

		my @alterations = split /,/ => $fields[3];

		for (my $i = 0; $i < @alterations; $i++) {
			if ($fields[2] eq $alterations[$i]) {
				warn "There is an alteration equal to the reference at '$snp' line $line. I will ignore it\n";
				splice @alterations, $i, 1;
			}
		}

		next unless @alterations;
#		my $num_of_alterations = scalar @alterations;

#		die "Fourth column has more alterations ($num_of_alterations) than the ploidy ($ploidy) into file '$snp' at line $line\n"
#			if $num_of_alterations > $ploidy;

#		die "More than one obervation ($fields[3]) for homozigosity into file '$snp' at line $line\n"
#			if $fields[4] eq 'HO' && $num_of_alterations > 1;

#		if ($fields[2] eq $fields[3]) {
#			warn "Reference and alteration patterns are equal. I will ignore this entry into file '$snp' at line $line\n";
#			next;
#		}

		if (not defined $indexed_snp{$fields[0]}) {
			$indexed_snp{$fields[0]} = IntervalTree->new;
		}

		my $tree = $indexed_snp{$fields[0]};

		# Sequence inside perl begins at 0
		my $position = $fields[1] - 1;

		# Compare the alterations and reference to guess the max variation on sequence
		my $size_of_variation = ( sort { $b <=> $a } map { length } @alterations, $fields[2] )[0];

		my $low = $position;
		my $high = $position + $size_of_variation - 1;

		my $nodes = $tree->search($low, $high);

		# My rules:
		# The biggest structural variation gains precedence.
		# If occurs an overlapping, I search to the biggest variation among
		# the saved alterations and compare it with the actual entry:
		#    *** Remove all overlapping variations if actual entry is bigger;
		#    *** Skip actual entry if it is lower than the biggest variation
		if (@$nodes) {
			# Catch the node with the biggest variation
			my $node = ( sort {($b->high - $b->low) <=> ($a->high - $a->low)} @$nodes )[0];
			my $size_of_variation_saved = $node->high - $node->low + 1;

			if ($size_of_variation_saved >= $size_of_variation) {
				my $data = $node->data;
				my $node_alt = join "," => @{ $data->{alt} };
				$node_alt =~ s/,$//;

				warn sprintf "Early alteration [%s %d %s %s %s] masks [%s %d %s %s %s] at '%s' line %d\n"
					=> $fields[0], $data->{pos}, $data->{ref}, $node_alt, $data->{plo}, $fields[0], $position, $fields[2], $fields[3], $fields[4], $snp, $line;

				next;
			} else {
				warn sprintf "Alteration [%s %d %s %s %s] masks early declaration(s) at '%s' line %d\n"
					=> $fields[0], $position, $fields[2], $fields[3], $fields[4], $snp, $line;

				# Remove all nodes overlapping
				$tree->delete($_->low, $_->high) for @$nodes;
			}
		}

		my %variation = (
			ref  => $fields[2],
			alt  => \@alterations,
			plo  => $fields[4],
			pos  => $position
		);

		$tree->insert($low, $high, \%variation);
	}

	close $fh;

	# Validate  indexed_snp on the indexed_fasta
	my $fasta_index = $fasta_index;

	for my $seq_id (keys %indexed_snp) {
		my $tree = delete $indexed_snp{$seq_id};

		my $seq = \$fasta_index->{$seq_id}{seq}
			or next;
		my $size = $fasta_index->{$seq_id}{size};

		my @nodes;
		my $catcher = sub { push @nodes => shift };

		$tree->inorder($catcher);
		my @to_remove;

		for my $node (@nodes) {
			my $data = $node->data;
			next if $data->{ref} eq '-' && $data->{pos} < $size;
			my $loc = index $$seq, $data->{ref}, $data->{pos};
			if ($loc == -1 || $loc != $data->{pos}) {
				warn "In validating '$snp': Not found reference '$data->{ref}' at fasta position $seq_id:",$data->{pos}+1,"\n";
				push @to_remove => $node;
			}
		}

		if (scalar @to_remove < scalar @nodes) {
			$indexed_snp{$seq_id} = $tree;
		}

		$tree->delete($_->low, $_->high) for (@to_remove);
	}
	
	return \%indexed_snp;
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
		$indexed_fasta{$_}{size} = length $indexed_fasta{$_}{seq};
	}

	unless (%indexed_fasta) {
		die "Error parsing '$fasta'. Maybe the file is empty\n";
	}

	$fh->close
		or die "Cannot close file $fasta: $!\n";

	return \%indexed_fasta;
}
