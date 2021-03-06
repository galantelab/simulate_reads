package TestsFor::App::Sandy::PieceTable;
# ABSTRACT: Tests for 'App::Sandy::PieceTable' class

use App::Sandy::Base 'test';
#use Data::Dumper;
use base 'TestsFor';

sub startup : Tests(startup) {
	my $test = shift;
	$test->SUPER::startup;

	my $class = ref $test;
	$class->mk_classdata('default_attr');
	$class->mk_classdata('default_seq');
	$class->mk_classdata('default_table');
}

sub setup : Tests(setup) {
	my $test = shift;
	$test->SUPER::setup;

	my $seq = "A large span of text";
	my %default_attr = (orig => \$seq);

	$test->default_seq(\$seq);
	$test->default_attr(\%default_attr);
	$test->default_table($test->class_to_test->new(%default_attr));
}

sub constructor : Tests(2) {
	my $test = shift;

	my $class = $test->class_to_test;
	my $table = $test->default_table;
	my %default_attr = %{ $test->default_attr };

	while (my ($attr, $value) = each %default_attr) {
		can_ok $table, $attr;
		is $table->$attr, $value,"The value for $attr shold be correct";
	}
}

sub delete : Test(4) {
	my $test = shift;

	my $table = $test->default_table;
	my $seq = $test->default_seq;

	# Try to remove "large "
	$table->delete(2, 6);
#	diag Dumper($table->piece_table);

	my @pieces = (
		{ pos => 0, len => 2  },
		{ pos => 8, len => 12 }
	);

	my $piece_table = $table->piece_table;

	for (my $i = 0; $i < @pieces; $i++) {
		is $pieces[$i]{pos}, $piece_table->[$i]{pos},
			"table[$i]: pos should be equal to $pieces[$i]{pos}";
		is $pieces[$i]{len}, $piece_table->[$i]{len},
			"table[$i]: len should be equal to $pieces[$i]{len}";
	}
}

sub delete_at_end : Test(3) {
	my $test = shift;

	my $table = $test->default_table;
	my $seq = $test->default_seq;

	# Try to delete the five last characters
	my $pos_last = length($$seq) - 5;
	$table->delete($pos_last, 5);

	my @pieces = (
		{ pos => 0, len => $pos_last },
	);

	my $piece_table = $table->piece_table;
	is @$piece_table, 1,
		"Remove at the end should return a piece_table with 1 piece";

	for (my $i = 0; $i < @pieces; $i++) {
		is $piece_table->[$i]{pos}, $pieces[$i]{pos},
			"table[$i]: pos should be equal to $pieces[$i]{pos}";
		is $piece_table->[$i]{len}, $pieces[$i]{len},
			"table[$i]: len should be equal to $pieces[$i]{len}";
	}
}

sub delete_at_start : Test(3) {
	my $test = shift;

	my $table = $test->default_table;
	my $seq = $test->default_seq;

	# Try to delete the first five characters
	my $len = 5;
	my $pos_start = 0;
	$table->delete($pos_start, $len);

	my @pieces = (
		{ pos => $len, len => length($$seq) - $len },
	);

	my $piece_table = $table->piece_table;
	is @$piece_table, 1,
		"Remove at the start should return a piece_table with 1 piece";

	for (my $i = 0; $i < @pieces; $i++) {
		is $piece_table->[$i]{pos}, $pieces[$i]{pos},
			"table[$i]: pos should be equal to $pieces[$i]{pos}";
		is $piece_table->[$i]{len}, $pieces[$i]{len},
			"table[$i]: len should be equal to $pieces[$i]{len}";
	}
}

sub insert : Test(6) {
	my $test = shift;

	my $table = $test->default_table;
	my $seq = $test->default_seq;

	# Try to insert 'English'
	my $add = "English ";
	$table->insert(\$add, 16);
#	diag Dumper($table->piece_table);

	my @pieces = (
		{ pos => 0,  len => 16 },
		{ pos => 16, len => 8  },
		{ pos => 16, len => 4  }
	);

	my $piece_table = $table->piece_table;

	for (my $i = 0; $i < @pieces; $i++) {
		is $pieces[$i]{pos}, $piece_table->[$i]{pos},
			"table[$i]: pos should be equal to $pieces[$i]{pos}";
		is $pieces[$i]{len}, $piece_table->[$i]{len},
			"table[$i]: len should be equal to $pieces[$i]{len}";
	}
}

sub insert_at_start : Test(5) {
	my $test = shift;

	my $table = $test->default_table;
	my $seq = $test->default_seq;

	# Try to change the start A to Ponga
	my $add = "Ponga ";
	$table->insert(\$add, 0);

	my @pieces = (
		{ pos => 0, len => 6 },
		{ pos => 0, len => length($$seq) }
	);

	my $piece_table = $table->piece_table;
	is @$piece_table, 2,
		"Insert at the start should return a piece_table with 2 piece";

	for (my $i = 0; $i < @pieces; $i++) {
		is $piece_table->[$i]{pos}, $pieces[$i]{pos},
			"table[$i]: pos should be equal to $pieces[$i]{pos}";
		is $piece_table->[$i]{len}, $pieces[$i]{len},
			"table[$i]: len should be equal to $pieces[$i]{len}";
	}
}

sub insert_at_end : Test(5) {
	my $test = shift;

	my $table = $test->default_table;
	my $seq = $test->default_seq;

	# Try to change the last five characters
	my $add = "Ponga";
	my $pos = length($$seq);
	$table->insert(\$add, length($$seq));

	my @pieces = (
		{ pos => 0,    len => $pos },
		{ pos => $pos, len => 5    }
	);

	my $piece_table = $table->piece_table;
	is @$piece_table, 2,
		"Insert at the end should return a piece_table with 2 piece";

	for (my $i = 0; $i < @pieces; $i++) {
		is $piece_table->[$i]{pos}, $pieces[$i]{pos},
			"table[$i]: pos should be equal to $pieces[$i]{pos}";
		is $piece_table->[$i]{len}, $pieces[$i]{len},
			"table[$i]: len should be equal to $pieces[$i]{len}";
	}
}

sub delete_and_insert : Test(8) {
	my $test = shift;

	my $table = $test->default_table;
	my $seq = $test->default_seq;

	# Try to remove "large "
	$table->delete(2, 6);

	# Try to insert 'English'
	my $add = "English ";
	$table->insert(\$add, 16);

#	diag Dumper($table->piece_table);

	my @pieces = (
		{ pos => 0,  len => 2  },
		{ pos => 8,  len => 8  },
		{ pos => 16, len => 8  },
		{ pos => 16, len => 4  }
	);

	my $piece_table = $table->piece_table;

	for (my $i = 0; $i < @pieces; $i++) {
		is $pieces[$i]{pos}, $piece_table->[$i]{pos},
			"table[$i]: pos should be equal to $pieces[$i]{pos}";
		is $pieces[$i]{len}, $piece_table->[$i]{len},
			"table[$i]: len should be equal to $pieces[$i]{len}";
	}
}

sub lookup : Test(5) {
	my $test = shift;

	my $table = $test->default_table;
	my $seq = $test->default_seq;

	# Try to remove "large "
	$table->delete(2, 6);

	# Try to insert 'English'
	my $add = "English ";
	$table->insert(\$add, 16);

	# I cannot forget to initialize the
	# logical offsets for the lookup method
	# to work
	$table->calculate_logical_offset;

	# Look for the pieces between 16 - 28:
	# "English text"
	my $pieces = $table->lookup(16, 12);

	is scalar(@$pieces), 2,
		"lookup returned the right number of pieces";

	my @pieces = (
		{ pos => 16, len => 8  },
		{ pos => 16, len => 4  }
	);

	for (my $i = 0; $i < @pieces; $i++) {
		is $pieces[$i]{pos}, $pieces->[$i]{pos},
			"table[$i]: pos should be equal to $pieces[$i]{pos}";
		is $pieces[$i]{len}, $pieces->[$i]{len},
			"table[$i]: len should be equal to $pieces[$i]{len}";
	}
}

sub change : Test(6) {
	my $test = shift;

	my $table = $test->default_table;
	my $seq = $test->default_seq;

	# Try to change "large " to "ponga "
	my $to_remove = "large ";
	my $change = "ponga ";
	$table->change(\$change, 2, length $to_remove);
#	diag Dumper($table->piece_table);

	my @pieces = (
		{ pos => 0, len => 2  },
		{ pos => 2, len => 6  },
		{ pos => 8, len => 12 }
	);

	my $piece_table = $table->piece_table;

	for (my $i = 0; $i < @pieces; $i++) {
		is $pieces[$i]{pos}, $piece_table->[$i]{pos},
			"table[$i]: pos should be equal to $pieces[$i]{pos}";
		is $pieces[$i]{len}, $piece_table->[$i]{len},
			"table[$i]: len should be equal to $pieces[$i]{len}";
	}
}

sub change_at_start : Test(5) {
	my $test = shift;

	my $table = $test->default_table;
	my $seq = $test->default_seq;

	# Try to change the start A to Ponga
	my $add = "Ponga ";
	$table->change(\$add, 0, 1);

	my @pieces = (
		{ pos => 0, len => 6 },
		{ pos => 1, len => length($$seq) - 1 }
	);

	my $piece_table = $table->piece_table;
	is @$piece_table, 2,
		"Change at the start should return a piece_table with 2 piece";

	for (my $i = 0; $i < @pieces; $i++) {
		is $piece_table->[$i]{pos}, $pieces[$i]{pos},
			"table[$i]: pos should be equal to $pieces[$i]{pos}";
		is $piece_table->[$i]{len}, $pieces[$i]{len},
			"table[$i]: len should be equal to $pieces[$i]{len}";
	}
}

sub change_at_end : Test(5) {
	my $test = shift;

	my $table = $test->default_table;
	my $seq = $test->default_seq;

	# Try to change the last five characters
	my $pos_last = length($$seq) - 5;
	my $add = "Ponga";
	$table->change(\$add, $pos_last, 5);

	my @pieces = (
		{ pos => 0,         len => $pos_last    },
		{ pos => $pos_last, len => length($add) }
	);

	my $piece_table = $table->piece_table;
	is @$piece_table, 2,
		"Change at the end should return a piece_table with 2 piece";

	for (my $i = 0; $i < @pieces; $i++) {
		is $piece_table->[$i]{pos}, $pieces[$i]{pos},
			"table[$i]: pos should be equal to $pieces[$i]{pos}";
		is $piece_table->[$i]{len}, $pieces[$i]{len},
			"table[$i]: len should be equal to $pieces[$i]{len}";
	}
}
