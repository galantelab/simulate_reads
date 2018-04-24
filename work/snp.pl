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
my $seed = shift || time();

print "My seed $seed\n";
srand $seed;

my $fasta_index = _index_fasta($fasta);
my $snp_index = index_snp($snp);

my $code = sub {
	my $node = shift;
	printf "[%d - %d] max = %d shift %d\n" => $node->low, $node->high, $node->max, $node->data->{sft};
};

print "Finito!\n";

for my $seq_id (sort keys %$snp_index) {
	my $tree = $snp_index->{$seq_id};
	print "chr = $seq_id\n";
	$tree->inorder($code);
}

my $pos = 6;
my $read_r = subseq(\$fasta_index->{chr1}{seq}, $fasta_index->{chr1}{size}, READ_SIZE, $pos);
print "read = $$read_r\nlength = ", length $$read_r, "\n";

my $end = $pos + READ_SIZE - 1;
my $variations = $snp_index->{chr1}->search($pos, $end);

my $is_ref = int(rand(2));
print "Is ref? $is_ref\n";

# The @$variations comes sorted by low.
# My appprouch is catch the last one and apply the changes
# on the sequence, next catch the penult and so on. This way
# I avoid to implement a shift tracker
for my $variation (reverse @$variations) {
	my $data = $variation->data;

	if ($is_ref && $data->{plo} eq 'HE') {
		next;
	}

	print Dumper($data);

	# TODO: FIx HERE ->
	# There is a problem with strctural variations. They shift
	# the sequence so that it cannot find the next changes by position
	# I think to solve, I need to make a tracker for the shift and when
	# testing and substr the sequence, I correct  the position.
	# If the user request a big deletion, the read will be truncated and 
	# to solve the algorithm catch the rest needed from the reference fasta.
	# But if really it is a big deletion, or a deletion at the end of the read,
	# It cannot be sovlved and the deletion will be "truncated".
	my $position = $data->{pos} - $pos;
	my $reference = $data->{ref};
	my $alteration = $data->{alt};

	# I need to test, again, if the reference pattern
	# still is there, because users may pass variations
	# too close that they overlap. In short, it causes
	# a race condition-like when the first one to occurs
	# gain the precedence.

	# Test before the 4 possibilities of intersections:
	#          [================] R
	#      [-------] V
	#                       [------] V
	#              [-----] V        
	#      [-----------------------] V
	
	my $length = 0;
	my $ref_pos = 0;

	if ($variation->low < $pos && $variation->high <= $end) {
		$position = 0;
		$length = $variation->high - $pos + 1;
		$ref_pos = $pos - $variation->low;
	} elsif ($variation->high > $end && $variation->low > $pos) {
		$position = $variation->low - $pos;
		$length = $end - $variation->low + 1;
	} elsif ($variation->high <= $end && $variation->low >= $pos) {
		$position = $variation->low - $pos;
		$length = $variation->high - $variation->low + 1;
	} else {
		$position = 0;
		$ref_pos = $pos - $variation->low;
		$length = $end - $pos + 1;
	}

	# Insertion
	if ($reference eq '-') {
#		substr($$read_r, $position, 0) = $alteration;
		substr($$read_r, $position, 0) = substr($alteration, $ref_pos, $length);

	# Deletion
	} elsif ($alteration eq '-') {
#		if (substr($$read_r, $position, length($reference)) eq $reference) {
#			substr($$read_r, $position, length($reference)) = "";
#		}
		if (substr($$read_r, $position, $length) eq substr($reference, $ref_pos, $length)) {
			substr($$read_r, $position, $length) = "";
		}

	# Change
	} else {
#		if (substr($$read_r, $position, $length) eq substr($reference, 0, $length)) {
#			substr($$read_r, $position, $length) = $alteration;
#		}
		if (substr($$read_r, $position, $length) eq substr($reference, $ref_pos, $length)) {
			substr($$read_r, $position, $length) = substr($alteration, $ref_pos, $length);
		}
	}
}

# Missing correction to the the new psition

my $new_read_size = length $$read_r;

if ($new_read_size < READ_SIZE) {
	my $missing = READ_SIZE - $new_read_size;
	my $beacon = $pos + READ_SIZE;
	substr($$read_r, $new_read_size, 0) = substr($fasta_index->{chr1}{seq}, $beacon, $missing);
	$new_read_size = length($$read_r);
	if ($new_read_size < READ_SIZE) {
		# May occur that the end of the fasta sequence suffers an deletion,
		# so I cannot guess how to fill the read size. In that case,
		# I fill the string with 'N's
		$missing = READ_SIZE - $new_read_size;
		substr($$read_r, $new_read_size, 0) = 'N' x $missing;
	}
} elsif ($new_read_size > READ_SIZE) {
	# Just truncate!
	$$read_r = substr($$read_r, 0, READ_SIZE);
}

print "read = $$read_r\nlength = ", length $$read_r, "\n";

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
			pos  => $position,
			sft  => 0
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

		# Validate node and at same time, set the shift factor,
		# in order to correct the reference genome position
		my $acm = 0;

		for my $node (@nodes) {
			my $data = $node->data;

			if ($data->{ref} ne '-' || $data->{pos} >= $size) {
				my $loc = index $$seq, $data->{ref}, $data->{pos};

				if ($loc == -1 || $loc != $data->{pos}) {
					warn "In validating '$snp': Not found reference '$data->{ref}' at fasta position $seq_id:",$data->{pos}+1,"\n";
					push @to_remove => $node;
					next;
				}
			}

			# Update shift tracker
			$acm += $node->high == $node->low
				? 0
				: length($data->{ref}) < length($data->{alt})
					? $node->high - $node->low + 1
					: $node->low - $node->high - 1;

			$data->{sft} = $acm;
		}

		my $new_size = $size + $acm;

		if ($new_size < READ_SIZE) {
			die "So many deletions on '$seq_id' resulted in a sequence lesser than the required read-size";
		}

		# I need to se the new weight based on the new size,
		# according to he INDEL patterns found
		$fasta_index->{$seq_id}{weight} = $new_size;

		if (scalar @to_remove < scalar @nodes) {
			$indexed_snp{$seq_id} = $tree;
		}

		$tree->delete($_->low, $_->high) for @to_remove;
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
