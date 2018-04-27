package App::SimulateReads::Fastq::SingleEnd;
# ABSTRACT: App::SimulateReads::Fastq subclass for simulate single-end fastq entries.

use App::SimulateReads::Base 'class';
use App::SimulateReads::Read::SingleEnd;

extends 'App::SimulateReads::Fastq';

with 'App::SimulateReads::Role::RunTimeTemplate';

# VERSION

has 'template_id' => (
	is         => 'ro',
	isa        => 'Str',
	required   => 1
);

has 'sequencing_error' => (
	is         => 'ro',
	isa        => 'My:NumHS',
	required   => 1
);

has '_gen_header' => (
	is         => 'ro',
	isa        => 'CodeRef',
	builder    => '_build_gen_header',
	lazy_build => 1
);

has '_info' => (
	traits     => ['Hash'],
	is         => 'ro',
	isa        => 'HashRef[Str]',
	builder    => '_build_info',
	lazy_build => 1,
	handles    => {
		_set_info => 'set',
		_get_info => 'get'
	}
);

has '_read' => (
	is         => 'ro',
	isa        => 'App::SimulateReads::Read::SingleEnd',
	builder    => '_build_read',
	lazy_build => 1,
	handles    => [qw{ gen_read }]
);

sub BUILD {
	my $self = shift;
	## Just to ensure that the lazy attributes are built before &new returns
	$self->_read;
}

sub _build_read {
	my $self = shift;
	App::SimulateReads::Read::SingleEnd->new(
		sequencing_error => $self->sequencing_error,
		read_size        => $self->read_size
	);
}

sub _build_gen_header {
	my $self = shift;
	my %sym_table = (
		'%q' => '$info->{quality_profile}',
		'%r' => '$info->{read_size}',
		'%e' => '$info->{sequencing_error}',
		'%c' => '$info->{seq_id}',
		'%t' => '$info->{start}',
		'%n' => '$info->{end}',
		'%i' => '$info->{instrument}',
		'%I' => '$info->{id}',
		'%R' => '$info->{read}',
		'%U' => '$info->{num}',
		'%s' => '$info->{strand}',
		'%X' => '$info->{error}',
		'%a' => '$info->{start_alt}',
		'%b' => '$info->{end_alt}',
		'%v' => '$info->{var}',
		'%z' => '$info->{ref}'
	);

	return  $self->compile_template($self->template_id, 'info', \%sym_table);
}

sub _build_info {
	my $self = shift;

	my %info = (
		instrument       => 'SR',
		quality_profile  => $self->quality_profile,
		read_size        => $self->read_size,
		sequencing_error => $self->sequencing_error
	);

	return \%info;
}

sub sprint_fastq {
	my ($self, $id, $num, $seq_id, $seq_ref, $seq_size, $is_leader, $tree) = @_;

	my ($read_ref, $detail_h) = $self->gen_read($seq_ref, $seq_size, $is_leader, $tree);

	# Set start - end regarding reference fasta
	my ($start, $end) = ($detail_h->{pos} + 1, $detail_h->{pos} + $self->read_size);

	# If there are no structural variations, set default values
	my ($start_alt, $end_alt) = ($start, $end);
	my ($ref, $vars) = ('ref', 'none');

	my $var_h = $detail_h->{var};

	# Set  structural variation details
	if ($var_h) {
		# Set start - end regarding modified fasta
		($start_alt, $end_alt) = ($var_h->{alt_pos} + 1, $var_h->{alt_pos} + $self->read_size);

		# Set if this is the alternative seq_id
		$ref = 'alt' unless $var_h->{is_ref};

		my $vars_a = $var_h->{vars};

		# Set variations string if any
		if (@$vars_a) {
			$vars = join ","
				=> map { sprintf "%d:%s/%s:%s" => $_->{rpos} + 1, $_->{ref}, $_->{alt}, $_->{plo} }
				@$vars_a;
		}
	}

	# Set defaut sequencing errors
	my $errors = 'none';
	my $error_a = $detail_h->{error};

	# Set errors if there are sequencing errors
	if (@$error_a) {
		$errors = join ","
			=> map { sprintf "%d:%s/%s" => $_->{pos} + 1, $_->{b}, $_->{not_b} }
			@$error_a;
	}

	unless ($is_leader) {
		($start, $end) = ($end, $start);
		($start_alt, $end_alt) = ($end_alt, $start_alt);
	}

	$self->_set_info(
		'id'        => $id,
		'num'       => $num,
		'seq_id'    => $seq_id,
		'start'     => $start,
		'end'       => $end,
		'read'      => 1,
		'strand'    => $is_leader ? 'P' : 'M',
		'error'     => $errors,
		'ref'       => $ref,
		'start_alt' => $start_alt,
		'end_alt'   => $end_alt,
		'var'       => $vars
	);

	my $gen_header = $self->_gen_header;
	my $header = $gen_header->($self->_info);

	return $self->fastq_template(\$header, $read_ref);
}
