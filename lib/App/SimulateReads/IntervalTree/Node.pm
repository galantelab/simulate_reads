package App::SimulateReads::IntervalTree::Node;
# ABSTRACT: Node class IntervalTree

use App::SimulateReads::Base 'class';

# VERSION

has [qw/low high/] => (
	is         => 'rw',
	isa        => 'Int',
	required   => 1
);

has 'max' => (
	is         => 'rw',
	isa        => 'Int',
	builder    => '_build_max',
	lazy_build => 1
);

has 'height' => (
	is         => 'rw',
	isa        => 'Int',
	default    => 1
);

has 'data' => (
	is         => 'rw',
	isa        => 'Maybe[Ref]',
	required   => 0
);

has [qw/left right/] => (
	is         => 'rw',
	isa        => 'Maybe[IntervalTree::Node]',
	required   => 0
);

sub _build_max {
	my $self = shift;
	return $self->high;
}
