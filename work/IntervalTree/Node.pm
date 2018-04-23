package IntervalTree::Node;

use Moose;
use MooseX::StrictConstructor;
use MooseX::UndefTolerant;
use namespace::autoclean;

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

has 'hidden' => (
	is         => 'rw',
	isa        => 'Bool',
	default    => 0
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

__PACKAGE__->meta->make_immutable;

1;
