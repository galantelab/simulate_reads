package App::SimulateReads::Command::Simulate::Custom;
# ABSTRACT: simulate subcommand class. Simulate a custom sequencing

use App::SimulateReads::Base 'class';

extends 'App::SimulateReads::Command::Simulate';

with 'App::SimulateReads::Role::Digest';

# VERSION

sub default_opt {
	'paired-end-id'    => '%i.%U %U',
	'single-end-id'    => '%i.%U %U',
	'seed'             => time,
	'verbose'          => 0,
	'prefix'           => 'out',
	'output-dir'       => '.',
	'jobs'             => 1,
	'gzip'             => 1,
	'count-loops-by'   => 'coverage',
	'coverage'         => 8,
	'strand-bias'      => 'random',
	'seqid-weight'     => 'length',
	'sequencing-type'  => 'paired-end',
	'fragment-mean'    => 300,
	'fragment-stdd'    => 50,
	'sequencing-error' => 0.005,
	'read-size'        => 101,
	'quality-profile'  => 'hiseq'
}

sub rm_opt {
}

__END__

=head1 SYNOPSIS

 simulate_reads simulate custom [options] <fasta-file>

 Arguments:
  a fasta-file 

 Options:
  -h, --help               brief help message
  -M, --man                full documentation
  -v, --verbose            print log messages
  -p, --prefix             prefix output [default:"out"]	
  -o, --output-dir         output directory [default:"."]
  -i, --append-id          append to the defined template id [Format]
  -I, --id                 overlap the default template id [Format]
  -j, --jobs               number of jobs [default:"1"; Integer]
  -z, --gzip               compress output file
  -s, --seed               set the seed of the base generator
  -c, --coverage           fastq-file coverage [default:"8", Number]
  -n, --number-of-reads    directly set the number of reads [Integer]
  -t, --sequencing-type    single-end or paired-end reads
                           [default:"paired-end"]
  -q, --quality-profile    illumina sequencing system profiles
                           [default:"hiseq"]
  -e, --sequencing-error   sequencing error rate
                           [default:"0.005"; Number]
  -r, --read-size          the read size [default:"101"; Integer]
  -m, --fragment-mean      the mean size fragments for paired-end reads
                           [default:"300"; Integer]
  -d, --fragment-stdd      the standard deviation for fragment sizes
                           [default:"50"; Integer]
  -b, --strand-bias        which strand to be used: plus, minus and random
                           [default:"random"]
  -w, --seqid-weight       seqid raffle type: length, same, file
                           [default: "length"]
  -f, --weight-file        an expression-matrix when seqid-weight=file

=head1 OPTIONS

=over 8

=item B<--help>

Print a brief help message and exits.

=item B<--man>

Prints the manual page and exits.

=item B<--verbose>

Prints log information to standard error

=item B<--prefix>

Concatenates the prefix to the output-file name.

=item B<--output-dir>

Creates output-file inside output-dir. If output-dir
does not exist, it is created recursively

=item B<--append-id>

Append string template to the defined template id.
See B<Format>

=item B<--id>

Overlap the default defined template id:
I<single-end> %i.%U %U and I<paired-end> %i.%U %U
e.g. SR123.1 1
See B<Format>

=item B<Format>

A string B<Format> is a combination of literal and escape characters similar to the way I<printf> works.
That way, the user has the freedom to customize the fastq sequence identifier to fit her needs. Valid
escape characteres are:

Common escape characters

	Escape       Meaning
	------       ------------------------------------------
	%i   	     instrument id composed by SR + PID
	%I           job slot number
	%q           quality profile
	%e           sequencing error
	%R           read 1, or 2 if it is the paired-end mate
	%U           read number
	%r           read size
	%c           sequence id as chromossome, ref
	%s           read or fragment strand
	%t           read start position
	%n           read end position

Paired-end specific escape characters

	Escape       Meaning
	------       ------------------------------------------
	%T           mate read start position
	%N           mate read end position
	%D           distance between the paired-reads
	%m           fragment mean
	%d           fragment standard deviation
	%f           fragment size
	%S           fragment start position
	%E           fragment end position

=item B<--jobs>

Sets the number of child jobs to be created

=item B<--gzip>

Compress the output-file with gzip algorithm. It is
possible to pass --no-output-gzip if one wants
uncompressed output-file

=item B<--seed>

Sets the seed of the base generator. The ability to set the seed is
useful for those who want reproducible simulations. Pay attention to
the number of jobs (--jobs) set, because each job receives a different
seed calculated from the I<main seed>. So, for reproducibility, the
same seed set before needs the same number of jobs set before as well.

=item B<--read-size>

Sets the read size. For now the unique valid value is 101

=item B<--coverage>

Calculates the number of reads based on the sequence
coverage: number_of_reads = (sequence_size * coverage) / read_size

=item B<--number-of-reads>

Sets directly the number of reads desired. It overrides coverage,
in case the two options are given

=item B<--sequencing-type>

Sets the sequencing type to single-end or paired-end

=item B<--fragment-mean>

If the sequencing-type is set to paired-end, it sets the
fragment mean

=item B<--fragment-stdd>

If the sequencing-type is set to paired-end, it sets the
fragment standard deviation

=item B<--sequencing-error>

Sets the sequencing error rate. Valid values are between zero and one

=item B<--quality-profile>

Sets the illumina sequencing system profile for quality. For now, the unique
valid values are hiseq and poisson

=item B<--strand-bias>

Sets which strand to use to make a read. Valid options are plus, minus and
random - if you want to randomly calculte the strand for each read

=item B<--seqid-weight>

Sets the seqid (e.g. chromossome, ensembl id) raffle behavior. Valid options are
length, same and file. If it is set to 'same', all seqid receives the same weight
when raffling. If it is set to 'length', the seqid weight is calculated based on
the seqid sequence length. And finally, if it is set to 'file', the user must set
the option --weight-file. For details, see B<--weight-file>

=item B<--weight-file>

If --seqid-weight is set to file, then this option becomes mandatory. A valid
weight file is an expression-matrix file with 2 columns. The first column is
for the seqid and the second column is for the count. The counts will be treated as weights

=back

=head1 DESCRIPTION

Simulate a custom sequencing.

=cut
