#
#===============================================================================
#
#         FILE: SimulateGenome.pm
#
#  DESCRIPTION: Simulates single-end and pair-end reads for genome
#
#        FILES: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Thiago Miller (tmiller), tmiller@mochsl.org.br
# ORGANIZATION: Group of Bioinformatics
#      VERSION: 1.0
#      CREATED: 04/25/2017 04:50:54 PM
#     REVISION: ---
#===============================================================================

package SimulateGenome;

use Moose;
use MooseX::StrictConstructor;
use My::Types;
use Fastq::SingleEnd;
use Fastq::PairedEnd;
use Carp;
use File::Cat 'cat';
use Parallel::ForkManager;
use Try::Tiny;
use feature 'state';

use namespace::autoclean;

with qw/My::Role::WeightedRaffle My::Role::IO/;

#-------------------------------------------------------------------------------
#  Moose attributes
#-------------------------------------------------------------------------------
has 'threads'         => (is => 'ro', isa => 'My:IntGt0',      required => 1);
has 'prefix'          => (is => 'ro', isa => 'Str',            required => 1);
has 'output_gzipped'  => (is => 'ro', isa => 'Bool',           required => 1);
has 'genome_file'     => (is => 'ro', isa => 'My:Fasta',       required => 1);
has 'coverage'        => (is => 'ro', isa => 'My:NumGt0',      required => 0);
has 'number_of_reads' => (is => 'ro', isa => 'My:IntGt0',      required => 0);
has 'count_loops_by'  => (is => 'ro', isa => 'My:CountLoopBy', required => 1);
has 'strand_bias'     => (is => 'ro', isa => 'My:StrandBias',  required => 1);
has 'sequence_weight' => (is => 'ro', isa => 'My:SeqWeight',   required => 1);
has 'weight_file'     => (is => 'ro', isa => 'My:File',        required => 0);
has 'fastq'           => (
	is         => 'ro',
	isa        => 'Fastq::SingleEnd | Fastq::PairedEnd',
	required   => 1,
	handles    => { get_fastq => 'fastq' }
);
has '_strand'         => (
	is         => 'ro',
	isa        => 'CodeRef',
	builder    => '_build_strand',
	lazy_build => 1
);
has '_genome'         => (
	is         => 'ro',
	isa        => 'HashRef[HashRef]',
	builder    => '_build_genome',
	lazy_build => 1
);
has '_raffle'        => (
	is         => 'ro',
	isa        => 'CodeRef',
	builder    => '_build_raffle',
	lazy_build => 1
);

sub BUILD {
	my $self = shift;

	## Validate optional attributes and/or attributes that depends of another attribute

	# If sequence_weight is 'file', then weight_file must be defined
	if ($self->sequence_weight eq 'file' and not defined $self->weight_file) {
		croak "sequence_weight=file requires a weight_file\n";
	}

	# If count_loops_by is 'coverage', then coverage must be defined. Else if
	# it is equal to 'number_of_reads', then number_of_reads must be defined
	if ($self->count_loops_by eq 'coverage' and not defined $self->coverage) {
		croak "count_loops_by=coverage requires a coverage number\n";
	} elsif ($self->count_loops_by eq 'number_of_reads' and not defined $self->number_of_reads) {
		croak "count_loops_by=number_of_reads requires a number_of_reads number\n";
	}
	
	# Just to ensure that indexed variables are built before &new returns
	$self->weights;
	$self->_genome;
}
 
sub _build_strand {
	my $self = shift;
	if ($self->strand_bias eq 'plus') {
		return sub {1};
	} elsif ($self->strand_bias eq 'minus') {
		return sub {0};
	} else {
		return sub { int(rand(2)) };
	}
}

#===  CLASS METHOD  ============================================================
#        CLASS: SimulateGenome
#       METHOD: _build_genome (BUILDER)
#   PARAMETERS: ????
#      RETURNS: ????
#  DESCRIPTION: 
#       THROWS: no exceptions
#     COMMENTS: none
#     SEE ALSO: n/a
#===============================================================================
sub _build_genome {
	my $self = shift;
	my $indexed_fasta = $self->index_fasta($self->genome_file);
	my $err;
	for my $id (keys %$indexed_fasta) {
		my $index_size = $indexed_fasta->{$id}{size};
		my $read_size = $self->fastq->read_size;
		$err .= "seqid sequence length (>$id => $index_size) lesser than required read size ($read_size)\n"
			if $index_size < $read_size;
	}
	
	croak "Error parsing '" . $self->genome_file . "':\n$err" if defined $err;
	return $indexed_fasta;
} ## --- end sub _build_genome

#===  CLASS METHOD  ============================================================
#        CLASS: SimulateGenome
#       METHOD: _build_weights (BUILDER)
#   PARAMETERS: ????
#      RETURNS: ????
#  DESCRIPTION: 
#       THROWS: no exceptions
#     COMMENTS: none
#     SEE ALSO: n/a
#===============================================================================
sub _build_weights {
	my $self = shift;
	if ($self->sequence_weight eq 'length') {
		my %chr_size = map { $_, $self->_genome->{$_}{size} } keys %{ $self->_genome };
		return $self->calculate_weights(\%chr_size);
	} elsif ($self->sequence_weight eq 'file') {
		my $indexed_file = $self->index_weight_file($self->weight_file);
		croak "Error parsing '" . $self->weight_file . "': Maybe the file is empty\n"
			unless %$indexed_file;
		my $indexed_fasta = $self->_genome;
		my $err;
		for my $id (keys %$indexed_file) {
			if (not exists $indexed_fasta->{$id}) {
				$err .= "seqid '$id' not found in '" . $self->genome_file . "'\n";
			}
		}
		croak "Error in validating '" . $self->weight_file . "':\n" . $err
			if defined $err;
		return $self->calculate_weights($indexed_file);
	}
} ## --- end sub _build_weights

sub _build_raffle {
	my $self = shift;
	if ($self->sequence_weight eq 'same') {
		my @seqids = keys %{ $self->_genome };
		return sub { 
			state @id;
			@id = @seqids if not @id;
			return $id[int(rand(scalar @id))];
		};
	} else {
		return sub { $self->weighted_raffle };
	}
}

sub _calculate_number_of_reads {
	my $self = shift;
	my $number_of_reads = 0;

	if ($self->count_loops_by eq 'coverage') {
		# It is needed to calculate the genome size
		my $genome = $self->_genome;
		my $genome_size = 0;
		$genome_size += $genome->{$_}{size} for keys %{ $genome };
		# In case it is paired-end read, divide the number of reads by 2 because Fastq::PairedEnd class
		# returns 2 reads at time
		my $read_type_factor = ref $self->fastq eq 'Fastq::PairedEnd' ? 2 : 1;
		$number_of_reads = int(($genome_size * $self->coverage) / ($self->fastq->read_size * $read_type_factor));
		# Maybe the number_of_reads is zero. It may occur due to the low coverage and/or genome_size
		if ($number_of_reads <= 0) {
			croak "Number of reads is equal to zero: Check the variables:\n" .
				"genome size: $genome_size\n" .
				"coverage: " . $self->coverage . "\n" .
				"read size: " . $self->fastq->read_size . "\n";
		}
	} else {
		$number_of_reads = $self->number_of_reads;
	}

	return $number_of_reads;
}

#===  CLASS METHOD  ============================================================
#        CLASS: SimulateGenome
#       METHOD: run_simulation
#   PARAMETERS: ????
#      RETURNS: ????
#  DESCRIPTION: 
#       THROWS: no exceptions
#     COMMENTS: none
#     SEE ALSO: n/a
#===============================================================================
sub run_simulation {
	my $self = shift;
	my $genome = $self->_genome;

	# Calculate the number of reads to be generated
	my $number_of_reads = $self->_calculate_number_of_reads;

	# Function that returns strand by strand_bias
	my $strand = $self->_strand;

	# Function that returns seqid by sequence_weight
	my $seqid = $self->_raffle;

	# File to be generated
	my $file = $self->prefix . "_simulation_seq.fastq";

	# Forks
	my $number_of_threads = $self->threads;
	my @tmp_files;
	my $pm = Parallel::ForkManager->new($number_of_threads);

	for my $tid (1..$number_of_threads) {
		#-------------------------------------------------------------------------------
		# Inside parent
		#-------------------------------------------------------------------------------
		my $file_t = "$file.${$}_$tid";
		push @tmp_files => $file_t;
		my $pid = $pm->start and next;	

		#-------------------------------------------------------------------------------
		# Inside child 
		#-------------------------------------------------------------------------------
		# Set child seed
		my $seed = time + $$;
		srand($seed);

		my $number_of_reads_t = int($number_of_reads/$number_of_threads);
		# If it is the first thread, make it work on the leftover reads of int() truncation
		$number_of_reads_t += $number_of_reads % $number_of_threads
			if $tid == 1;

		my $fh = $self->my_open_w($file_t, 0);

		# Run simualtion in child
		for my $i (1..$number_of_reads_t) {
			my $chr = $seqid->();
			my $entry;
			try {
				$entry = $self->get_fastq("SR${i}.$tid", $chr, \$genome->{$chr}{seq}, $genome->{$chr}{size}, $strand->());
			} catch {
				croak "Not defined entry for seqid '>$chr' at job $tid: $_";
			} finally {
				$fh->say($entry) unless @_;
			};
		}

		$fh->close;
		$pm->finish;
	}

	# Back to parent
	$pm->wait_all_children;

	# Concatenate all temporary files
	my $fh = $self->my_open_w($file, $self->output_gzipped);
	for my $file_t (@tmp_files) {
		cat $file_t => $fh
			or carp "Cannot concatenate $file_t to $file: $!";
	}

	$fh->close;

	# Clean up the mess
	for my $file_t (@tmp_files) {
		unlink $file_t
			or carp "Cannot remove temporary file: $file_t: $!";
	}
} ## --- end sub run_simulation

__PACKAGE__->meta->make_immutable;

1; ## --- end class SimulateGenome
