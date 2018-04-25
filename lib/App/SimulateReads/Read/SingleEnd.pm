package App::SimulateReads::Read::SingleEnd;
# ABSTRACT: App::SimulateReads::Read subclass for simulate single-end reads.

use App::SimulateReads::Base 'class';
use constant MAX_REDO => 100;

extends 'App::SimulateReads::Read';

# VERSION

sub gen_read {
	my ($self, $seq_ref, $seq_size, $is_leader, $tree) = @_;

	if ($seq_size < $self->read_size) {
		die sprintf "seq_size (%d) must be greater or equal to read_size (%d)\n"
			=> $seq_size, $self->read_size;
	}

	my $redo = 0;

	REDO:

	if ($redo >= MAX_REDO) {
		die "Max redo achieved: So many tries to introduce strutural variations.\n",
		"It may occur beacuse of the stochastic nature of the program.\n",
		"Maybe there are so many deletions\n";
	}

	my ($read_ref, $read_pos) = $self->subseq_rand($seq_ref, $seq_size, $self->read_size);

	my $var_apply;

	if ($tree) {
		$var_apply = $self->insert_structural_variation($seq_ref, $seq_size,
			$read_ref, $self->read_size, $read_pos, $tree);

		unless ($var_apply) {
			# Try again if the sequence is truncated
			$redo++;
			goto REDO;
		}
	}

	unless ($is_leader) {
		$self->reverse_complement($read_ref);
	}

	$self->update_count_base($self->read_size);
	my $errors = $self->insert_sequencing_error($read_ref);

	my %detail = (
		'pos'    => $read_pos,
		'var'    => $var_apply,
		'error'  => $errors
	);

	return ($read_ref, \%detail);
}
