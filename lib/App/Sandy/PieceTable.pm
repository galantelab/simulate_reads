package App::Sandy::PieceTable;
# ABSTRACT: Implement a piece table data structure class

use App::Sandy::Base 'class';
use App::Sandy::BTree::Interval;

with 'App::Sandy::Role::BSearch';

# VERSION

has 'orig' => (
	is         => 'ro',
	isa        => 'ScalarRef',
	required   => 1
);

has 'len' => (
	is         => 'ro',
	isa        => 'My:IntGt0',
	lazy_build => 1,
	builder    => '_build_len'
);

has 'piece_table' => (
	traits     => ['Array'],
	is         => 'ro',
	isa        => 'ArrayRef[My:Piece]',
	lazy_build => 1,
	builder    => '_build_piece_table',
	handles    => {
		_get_piece     => 'get',
		_splice_piece  => 'splice',
		_count_pieces  => 'count',
		_all_píeces    => 'elements'
	}
);

has 'logical_offset' => (
	is         => 'rw',
	isa        => 'App::Sandy::BTree::Interval',
	default    => sub { App::Sandy::BTree::Interval->new },
	handles    => {
		_add_offset    => 'insert',
		_rm_offset     => 'delete',
		_search_offset => 'search'
	}
);

has 'logical_len' => (
	is         => 'rw',
	isa        => 'My:IntGt0',
	lazy_build => 1,
	builder    => '_build_logical_len'
);

sub _build_len {
	my $self = shift;
	my $orig = $self->orig;
	return length $$orig;
}

sub _build_logical_len {
	my $self = shift;
	return $self->len;
}

sub _build_piece_table {
	my $self = shift;
	my $piece = $self->_piece_new($self->orig, 1, 0, $self->len, 0);
	return [$piece];
}

sub _piece_new {
	my ($self, $ref, $is_orig, $start, $len, $pos, $annot) = @_;

	my $piece = {
		'ref'     => $ref,     # reference to sequence
		'is_orig' => $is_orig, # set 1 if it is the orig sequence
		'start'   => $start,   # start position at reference
		'len'     => $len,     # length
		'pos'     => $pos,     # position at original sequence
		'offset'  => 0,        # position at the changed sequence
		'annot'   => $annot    # custom annotation
	};

	return $piece;
}

sub insert {
	my ($self, $ref, $pos, $annot) = @_;

	# Test if the position is inside the original sequence boundary
	if ($pos > $self->len) {
		croak "Trying to insert outside the original sequence";
	}

	# My length
	my $len = length $$ref;

	# Create piece data
	my $new_piece = $self->_piece_new($ref, 0, 0, $len, $pos, $annot);

	# Insert at end position or split
	# piece found at position 'pos'.
	my $index = $pos == $self->len
		? $self->_count_pieces
		: $self->_split_piece($pos);

	# Then insert new_piece
	$self->_splice_piece($index, 0, $new_piece);
}

sub delete {
	my ($self, $pos, $len, $annot) = @_;

	# Test if the removed region is inside the original sequence boundary
	if (($pos + $len) > $self->len) {
		croak "Trying to delete a region outside the original sequence";
	}

	# SPECIAL CASE: Delete at the very end
	if (($pos + $len) == $self->len) {
		return $self->_delete_at_end($len);
	}

	# Split piece at $pos. It will correctly fix the original
	# piece before the split and insert a new piece afterward.
	# So I need to catch tha last and fix the start and len fields
	my $index = $self->_split_piece($pos);
	my $piece = $self->_get_piece($index);

	# Fix position and len
	my $new_start = $pos + $len;
	my $new_len = $piece->{len} - $len;

	# Update!
	$piece->{start} = $piece->{pos} = $new_start;
	$piece->{len} = $new_len;
	$piece->{annot} = $annot;

	# If the new len is zero, then remove
	# this piece
	if ($new_len == 0) {
		$self->_splice_piece($index, 1);
	}
}

sub _delete_at_end {
	my ($self, $len) = @_;

	# Just catch the last piece and remove the
	# last sequence length
	my $index = $self->_count_pieces - 1;
	my $piece = $self->_get_piece($index);

	my $new_len = $piece->{len} - $len;

	# Update!
	$piece->{len} = $new_len;

	# If the new len is zero, then remove
	# this piece
	if ($new_len == 0) {
		$self->_splice_piece($index, 1);
	}
}

sub change {
	# A delete and insert operations.
	# delete from pos until len and insert ref at pos
	my ($self, $ref, $pos, $len, $annot) = @_;

	# Test if the changing region is inside the original sequence boundary
	if (($pos + $len) > $self->len) {
		croak "Trying to change a region outside the original sequence";
	}

	# My length
	my $ref_len = length $$ref;

	# Create piece data
	my $new_piece = $self->_piece_new($ref, 0, 0, $ref_len, $pos, $annot);

	# SPECIAL CASE: Change at the very end
	if (($pos + $len) == $self->len) {
		return $self->_change_at_end($new_piece, $len);
	}

	# Split piece found at position 'pos'.
	# Update old piece, insert piece and return
	# index where I can find the piece to remove
	# from and insert the change
	my $index = $self->_split_piece($pos);

	# Catch the piece from where I will remove
	my $piece = $self->_get_piece($index);

	# Fix position and len
	my $new_start = $pos + $len;
	my $new_len = $piece->{len} - $len;

	# Update!
	$piece->{start} = $piece->{pos} = $new_start;
	$piece->{len} = $new_len;

	# If the new len is zero, then remove
	# this piece
	if ($new_len == 0) {
		$self->_splice_piece($index, 1);
	}

	# Then insert new_piece
	$self->_splice_piece($index, 0, $new_piece);
}

sub _change_at_end {
	my ($self, $new_piece, $len) = @_;

	# Just catch the last piece and remove the
	# last sequence length
	my $index = $self->_count_pieces - 1;
	my $piece = $self->_get_piece($index);

	my $new_len = $piece->{len} - $len;
	$piece->{len} = $new_len;

	# If the new len is zero, then remove
	# this piece
	if ($new_len == 0) {
		$self->_splice_piece($index, 1);
	}

	# Then insert new_piece
	$self->_splice_piece($index + 1, 0, $new_piece);
}

sub calculate_logical_offset {
	# Before lookup() it is necessary to calculate
	# the positions according to the shift caused by
	# the structural variations. It will be used to
	# feed a binary tree
	my $self = shift;

	# Remove all old entries, if any, and create a new tree
	$self->logical_offset(App::Sandy::BTree::Interval->new);

	my $offset_acm = 0;

	# Insert each piece reference into a tree
	for my $piece ($self->_all_píeces) {
		# Calculate corrected piece boundaries
		my $low = $offset_acm;
		my $high = $offset_acm + $piece->{len} - 1;

		# Update piece offset
		$piece->{offset} = $offset_acm;

		# Update offset acumulator
		$offset_acm += $piece->{len};

		# Insert piece into tree
		$self->_add_offset($low, $high, $piece);
	}

	# Update logical offset with the corrected length
	$self->logical_len($offset_acm);
}

sub lookup {
	# Run 'calculate_logical_offset' before
	my ($self, $start, $len) = @_;

	# Calculate high boundary
	my $end = $start + $len - 1;

	# Return all pieces overlapping the boundaries
	return $self->_search_offset($start, $end);
}

sub _split_piece {
	my ($self, $pos) = @_;

	# Catch orig index where pos is inside
	my $index = $self->_piece_at($pos);

	# Get piece which will be updated
	my $old_piece = $self->_get_piece($index);

	# Split at the beggining of a piece,
	# or this piece has length 1
	if ($pos == $old_piece->{start}) {
		return $index;
	}

	# Calculate piece end
	my $old_end = $old_piece->{start} + $old_piece->{len} - 1;

	# Calculate the corrected length according to the split
	my $new_len = $pos - $old_piece->{start};

	# Update piece
	$old_piece->{len} = $new_len;

	# Create the second part of the split after the break position
	my $piece = $self->_piece_new($old_piece->{ref}, $old_piece->{is_orig},
		$pos, $old_end - $pos + 1, $pos);

	# Insert second part after updated piece
	$self->_splice_piece(++$index, 0, $piece);

	# return corrected index that resolves to
	# the position between the breaked piece
	return $index;
}

sub _is_pos_inside_piece {
	my ($self, $pos, $piece) = @_;
	my $end = $piece->{pos} + $piece->{len} - 1;
	return $pos >= $piece->{pos} && $pos <= $end;
}

sub _piece_at {
	my ($self, $pos) = @_;

	# State the function to compare at bsearch
	state $func = sub {
		my ($pos, $piece) = @_;
		if ($self->_is_pos_inside_piece($pos, $piece)) {
			return 0;
		} elsif ($pos > $piece->{pos}) {
			return 1;
		} else {
			return -1;
		}
	};

	# Search the piece index where $pos is inside the boundaries
	my $index = $self->with_bsearch($pos, $self->piece_table,
		$self->_count_pieces, $func);

	# Maybe it is undef. I need to take care to not
	# search to a position that was removed before.
	# I can avoid it when parsing the snv file
	if (not defined $index) {
		croak "Not found pos = $pos into piece_table. Maybe the region was removed?";
	}

	# Catch the piece at index
	my $piece = $self->_get_piece($index);

	# If I catched a non original sequence, then it must
	# be afterward
	if (not $piece->{is_orig}) {
		$piece = $self->_get_piece(++$index);
		unless ($self->_is_pos_inside_piece($piece, $pos)) {
			croak "Position is not inside the piece after non original sequence";
		}
	}

	return $index;
}
