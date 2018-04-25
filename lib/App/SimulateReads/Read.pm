package App::SimulateReads::Read;
# ABSTRACT: Base class to simulate reads

use App::SimulateReads::Base 'class';

# VERSION

has 'sequencing_error' => (
	is         => 'ro',
	isa        => 'My:NumHS',
	required   => 1
);

has 'read_size' => (
	is         => 'ro',
	isa        => 'My:IntGt0',
	required   => 1
);

has '_count_base' => (
	is         => 'rw',
	isa        => 'Int',
	default    => 0
);

has '_base' => (
	is         => 'rw',
	isa        => 'Int',
	builder    => '_build_base',
	lazy_build => 1
);

has '_not_base' => (
	is         => 'ro',
	isa        => 'HashRef',
	builder    => '_build_not_base',
	lazy_build => 1
);

sub _build_not_base {
	my %not_base = (
		A => ['T', 'C', 'G'],
		a => ['t', 'c', 'g'],
		T => ['A', 'C', 'G'],
		t => ['a', 'c', 'g'],
		C => ['A', 'T', 'G'],
		c => ['a', 't', 'g'],
		G => ['A', 'T', 'C'],
		g => ['a', 't', 'c']
	);
	return \%not_base;
}

sub _build_base {
	my $self = shift;
	# If sequencing_error equal to zero, set _base to zero
	return $self->sequencing_error && int(1 / $self->sequencing_error);
}

sub subseq {
	my ($self, $seq_ref, $seq_len, $slice_len, $pos) = @_;
	my $read = substr $$seq_ref, $pos, $slice_len;
	return \$read;
}

sub subseq_rand {
	my ($self, $seq_ref, $seq_len, $slice_len) = @_;
	my $usable_len = $seq_len - $slice_len;
	my $pos = int(rand($usable_len + 1));
	my $read = substr $$seq_ref, $pos, $slice_len;
	return (\$read, $pos);
}

sub insert_sequencing_error {
	my ($self, $seq_ref) = @_;

	my $err = int($self->_count_base * $self->sequencing_error);
	my @errors;

	for (my $i = 0; $i < $err; $i++) {
		$self->update_count_base(-$self->_base);
		my $pos = $self->read_size - $self->_count_base - 1;

		my $b = substr($$seq_ref, $pos, 1);
		my $not_b = $self->_randb($b); 

		substr($$seq_ref, $pos, 1) = $not_b;
		push @errors => { b => $b, not_b => $not_b, pos => $pos };
	}

	return \@errors;
}

sub update_count_base {
	my ($self, $val) = @_;
	$self->_count_base($self->_count_base + $val);
}

sub reverse_complement {
	my ($self, $seq_ref) = @_;
	$$seq_ref = reverse $$seq_ref;
	$$seq_ref =~ tr/atcgATCG/tagcTAGC/;
}

sub _randb {
	my ($self, $base) = @_;
	return $self->_not_base->{$base}[int(rand(3))] || $base;
}

#sub insert_structural_variation {
#	my ($self, $src_ref, $src_len, $seq_ref, $seq_len, $pos, $tree) = @_;
#
#	my $max_redo = 100;
#	my $redo = 0;
#
#	# Try again if the sequence is truncated
#	REDO:
#
#	if ($redo >= 100) {
#		die "Max redo achieved: So many tries to introduce strutural variations. ",
#		"It may occur beacuse of the stochastic nature of the program. ",
#		"Maybe there are so many deletions\n";
#	}
#
#	my $end = $pos + $seq_len - 1;
#	my $variations = $tree->search($pos, $end);
#	my @var_detail;
#
#	# Raffle if this is the reference sequence
#	my $is_ref = int(rand(2));
#
#	# The @$variations comes sorted by low.
#	# My appprouch is catch the last one and apply the changes
#	# on the sequence, next catch the penult and so on. This way
#	# I avoid to implement a shift tracker
#	for my $variation (reverse @$variations) {
#		my $data = $variation->data;
#
#		if ($is_ref && $data->{plo} eq 'HE') {
#			next;
#		}
#
#		# Test before the 4 possibilities of intersections:
#		my ($alt_pos, $ref_pos, $length);
#
#		#          [================] R
#		#      [-------] V
#		if ($variation->low < $pos && $variation->high <= $end) {
#			$ref_pos = 0;
#			$alt_pos = $pos - $variation->low;
#			$length = $variation->high - $pos + 1;
#
#		#          [================] R
#		#                       [------] V
#		} elsif ($variation->low >= $pos && $variation->high > $end) {
#			$ref_pos = $variation->low - $pos;
#			$alt_pos = 0;
#			$length = $end - $variation->low + 1;
#
#		#          [================] R
#		#              [-----] V        
#		} elsif ($variation->low >= $pos && $variation->high <= $end) {
#			$ref_pos = $variation->low - $pos;
#			$alt_pos = 0;
#			$length = $variation->high - $variation->low + 1;
#
#		#          [================] R
#		#      [-----------------------] V
#		} else {
#			$ref_pos = 0;
#			$alt_pos = $pos - $variation->low;
#			$length = $end - $pos + 1;
#		}
#
#		my $reference = $data->{ref};
#		my $alteration = $data->{alt};
#
#		# Insertion
#		if ($reference eq '-') {
#			substr($$seq_ref, $ref_pos, 0) = substr($alteration, $alt_pos, $length);
#
#		# Deletion
#		} elsif ($alteration eq '-') {
#			if (substr($$seq_ref, $ref_pos, $length) eq substr($reference, $alt_pos, $length)) {
#				substr($$seq_ref, $ref_pos, $length) = "";
#			}
#
#		# Change
#		} else {
#			if (substr($$seq_ref, $ref_pos, $length) eq substr($reference, $alt_pos, $length)) {
#				substr($$seq_ref, $ref_pos, $length) = substr($alteration, $alt_pos, $length);
#			}
#		}
#
#		my %var = (
#			'ref' => $reference,
#			'alt' => $alteration,
#			'pos' => $ref_pos
#		);
#
#		push @var_detail => \%var;
#	}
#
#	# Missing correction to the the new psition
#
#	my $new_size = length $$seq_ref;
#
#	if ($new_size < $seq_len) {
#		my $missing = $seq_len - $new_size;
#		my $beacon = $pos + $seq_len;
#		# I need a recursive routine
#		substr($$seq_ref, $new_size, 0) = substr($$src_ref, $beacon, $missing);
#		$new_size = length($$seq_ref);
#		if ($new_size < $seq_len) {
#			# May occur that the end of the fasta sequence suffers an deletion,
#			# so I cannot guess how to fill the read size. In that case,
#			# In theat case, I raffle another position
#			($seq_ref, $pos) = $self->subseq_rand($src_ref, $src_len, $seq_len);
#			$redo ++;
#			goto REDO;
#		}
#	} elsif ($new_size > $seq_len) {
#		# Just truncate!
#		$$seq_ref = substr($$seq_ref, 0, $seq_len);
#	}
#
#	# Correct position
#	my $alt_pos = $variations->[0]->data->{lsht} + $pos;
#
#	my %detail = (
#		'ref_pos' => $pos,
#		'alt_pos' => $alt_pos < 0 ? 0 : $alt_pos,
#		'var'     => \@var_detail,
#		is_ref
#	);
#
#	return \%detail;
#}

sub insert_structural_variation {
	my ($self, $src_ref, $src_len, $seq_ref, $seq_len, $pos, $tree) = @_;

	# Raffle if this is the reference sequence
	my $is_ref = int(rand(2));

	my $var_apply = $self->_insert_structural_variation($src_ref, $src_len,
		$seq_ref, $seq_len, $pos, $tree, $is_ref);

	my $new_size = length $$seq_ref;
	my $shift = 0;

	while ($new_size < $seq_len) {
		my $missing = $seq_len - $new_size;
		my $beacon = $pos + $seq_len + $shift;

		my $seq_missing = substr $$src_ref, $beacon, $missing;
		my $size_missing = length $seq_missing;

		if (($new_size + $size_missing) < $seq_len) {
			# May occur that the end of the fasta sequence suffers an deletion,
			# so I cannot guess how to fill the read size. In that case,
			# I return a undef value
			return;
		}

		my $var_apply_missing = $self->_insert_structural_variation($src_ref, $src_len,
			\$seq_missing, $size_missing, $beacon, $tree, $is_ref);
		
		$$seq_ref .= $seq_missing;
		$new_size = length $$seq_ref;
		$shift += $size_missing;
		push @$var_apply => @$var_apply_missing;
	}

	if ($new_size > $seq_len) {
		# Just truncate!
		$$seq_ref = substr($$seq_ref, 0, $seq_len);
	}

	# TODO: Calculate alt_pos for all!
	# Correct position
	my $alt_pos = @$var_apply
		? $var_apply->[0]->data->{lsft} + $pos
		: $pos;

	my @var = map { $_->data } @$var_apply;

	my %detail = (
		'alt_pos' => $alt_pos < 0 ? 0 : $alt_pos,
		'is_ref'  => $is_ref,
		'var'     => \@var
	);

	return \%detail;
}

sub _insert_structural_variation {
	my ($self, $src_ref, $src_len, $seq_ref, $seq_len, $pos, $tree, $is_ref) = @_;

	my $end = $pos + $seq_len - 1;
	my $variations = $tree->search($pos, $end);
	my @var_apply;

	# The @$variations comes sorted by low.
	# My appprouch is catch the last one and apply the changes
	# on the sequence, next catch the penult and so on. This way
	# I avoid to implement a shift tracker
	for my $variation (reverse @$variations) {
		my $data = $variation->data;

		if ($is_ref && $data->{plo} eq 'HE') {
			next;
		}

		# Test before the 4 possibilities of intersections:
		my ($alt_pos, $ref_pos, $length);

		#          [================] R
		#      [-------] V
		if ($variation->low < $pos && $variation->high <= $end) {
			$ref_pos = 0;
			$alt_pos = $pos - $variation->low;
			$length = $variation->high - $pos + 1;

		#          [================] R
		#                       [------] V
		} elsif ($variation->low >= $pos && $variation->high > $end) {
			$ref_pos = $variation->low - $pos;
			$alt_pos = 0;
			$length = $end - $variation->low + 1;

		#          [================] R
		#              [-----] V        
		} elsif ($variation->low >= $pos && $variation->high <= $end) {
			$ref_pos = $variation->low - $pos;
			$alt_pos = 0;
			$length = $variation->high - $variation->low + 1;

		#          [================] R
		#      [-----------------------] V
		} else {
			$ref_pos = 0;
			$alt_pos = $pos - $variation->low;
			$length = $end - $pos + 1;
		}

		my $reference = $data->{ref};
		my $alteration = $data->{alt};

		# Insertion
		if ($reference eq '-') {
			substr($$seq_ref, $ref_pos, 0) = substr($alteration, $alt_pos, $length);

		# Deletion
		} elsif ($alteration eq '-') {
			if (substr($$seq_ref, $ref_pos, $length) eq substr($reference, $alt_pos, $length)) {
				substr($$seq_ref, $ref_pos, $length) = "";
			}

		# Change
		} else {
			if (substr($$seq_ref, $ref_pos, $length) eq substr($reference, $alt_pos, $length)) {
				substr($$seq_ref, $ref_pos, $length) = substr($alteration, $alt_pos, $length);
			}
		}

		unshift @var_apply => $variation;
	}

	return \@var_apply;
}
