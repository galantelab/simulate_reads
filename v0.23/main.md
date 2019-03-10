### Contents at a Glance ###

1. [New Features](#new-features)
2. [Installation](#installation)
3. [Usage and Option summary](#usage-and-option-summary)
	1. [General](#general-syntax)
	2. [Main Commands](#main-commands)
		1. [Command `genome`](#command-genome)
		2. [Command `transcriptome`](#command-transcriptome)
	3. [Database Commands](#database-commands)
		1. [Command `quality`](#command-quality)
		2. [Command `expression`](#command-expression)
		3. [Command `variation`](#command-variation)
	4. [Miscellaneous Commands](#miscellaneous-commands)
		1. [Command `version`](#command-version)
		2. [Command `citation`](#command-citation)
	5. [Help Commands](#help-commands)
4. [Docker Usage](#docker-usage)
5. [Case study examples](case.md#case-study-examples)



### New Features ###

This is the version 0.23, which implements these new interesting features:

This version implents these new interesting features:
* Several ways to customize the reads' identifiers in the FASTq files for on
demand outputs.
* Ready to use built-in datasets for quality profile of sequencing, based on
experimental data from several platforms (e.g. hiseq, miseq and nextseq from
Illumina).
* Ready to use datasets for 54 tissue especific matrix expression, from GTEx
[(Xena Project)](https://xena.ucsc.edu/).
* Simulation of 3rd generation sequencing reads (PacBio and Nanopore) with
specific quality profiles
* Addition of many types of genomic variations (SNPs, INDELs, fusions an others)
in an easy-to-use way.
* Improved way to record user’s experimental data to the datasets' entries.
* Run many instances of *SANDY* in a scalable way by pulling its Docker
[image](https://hub.docker.com/r/galantelab/sandy).



### Installation ###

#### Prerequisites ####

Along with **Perl**, the user must have **zlib**, **gcc** and **make** and
**perl-doc** packages installed. To install them according to your distro, use:

* Debian/Ubuntu
```bash
$ apt-get install perl zlib1g-dev gcc make perl-doc
```

* CentOS/Fedora
```bash
$ yum install perl zlib gcc make perl-doc
```

* Archlinux
```bash
$ pacman -S perl zlib gcc make perl-doc
```

*SANDY* uses the *Comprehensive Perl Archive Network*, CPAN, as its package
manager, which allows a good control over all dependencies needed. If the user
already has Perl installed, then he may have the cpan command utility too. At
the first run, cpan will interactively configure your environment and mirror.

If the user is not sure, confirm the aforementioned prerequisites and, after
this, install *cpanminus* package manager:
```bash
$ cpan -i App::cpanminus
```

`App::cpanminus` will provide the `cpanm` utility, which has the capability of
install not only *SANDY*, but also all its dependencies, recursively.


#### Installing *SANDY* properly ####

Finally install *SANDY* with:
```bash
$ cpanm App::Sandy
```

**Important:** MacOS users must add an extra option to the command above, like
this:
```bash
	$ cpanm --force App::Sandy
```

For more details, see the [INSTALL](https://github.com/galantelab/sandy/blob/master/INSTALL)
file on *SANDY's* GitHub [repository](https://github.com/galantelab/sandy).


#### Or get *SANDY* in a Docker image ####

If the user prefer to avoid any intallation process and has Docker, you can just
pull *SANDY's* [image](https://hub.docker.com/r/galantelab/sandy) from Docker Hub with:
```bash
	$ docker pull galantelab/sandy
```

And he will take the latest version of *SANDY*, ready to rock!
So, to view some instructions about how to use *SANDY* from a docker image, see
the manual or consult this web [tutorial about SANDY usage from docker](#docker-usage).

**Important:** Docker has some strict default configurations for memory and CPU
usage on MacOS. Users of this system can change these configurations on their
behalf by accessing the [Preferences menu](https://docs.docker.com/docker-for-mac/#preferences-menu)
on the Docker icon at top right corner of their desktops.



### Usage and option summary ###

#### General Syntax ####

**Usage:**
```bash
$ sandy [options]
$ sandy help <command>
$ sandy <command> [options] <FILEs>
```
where there are basically two commands for general help, two main commands
with their own inner options, tree database management commands, some
miscellaneous commands and a specific `help` command for each of the main
commands. See:

Options							| Description
-------------------------------	| ----------------------------------------------
-h, --help						| brief help message
-u, --man						| full documentation
**Help commands:**				| 
help							| show application or command-specific help
man								| show application or command-specific documentation
**Misc commands:**				| 
version							| print the current version
citation						| export citation in BibTeX format
**Database commands:**			| 
quality							| manage quality profile database
expression						| manage expression-matrix database
variation						| manage structural variation database
**Main commands:**				| 
genome							| simulate genome sequencing
transcriptome					| simulate transcriptome sequencing

**Some examples:**
1. If the user wants to see a brief help for the `expression` command, he can
type both of the commands bellow:
	```bash
	$ sandy help expression
	```
or
	```bash
	$ sandy expression -h
	```
2. To view the verion of *SANDY* in use, just type:
	```bash
	$ sandy version
	```
3. And, to take a BibTex entry to cite *SANDY* in his work, type:
	```bash
	$ sandy citation
	```


#### Main Commands ####

##### Command `genome` #####

Use it to generate simulated FASTq-files from a given FASTA-file.
The `genome` command sets these default options for a genome sequencing
simulation:
* The strand is **randomly** chosen;
* The number of reads is calculated by the coverage;
* The chromossomes are raffled following a weighted raffle with the
sequence length as the bias;

**Usage:**
```bash
$ sandy genome [options] <fasta-file>
```
whose options' exaustive list can be consulted by `sandy genome -h` or
even `sandy help genome` commands.
At least one fasta-file must be given as the `<FILEs>` term. The results
will be one or two fastaq-files, depending on the sequencing-type option,
`-t`, for single-ended or paired-ended reads, and an additional file for
the reads-counts.

Options							| Description
-------------------------------	| ----------------------------------------------
-h, --help						| brief help message
-u, --man						| full documentation
-v, --verbose					| print log messages
-p, --prefix					| prefix output [default:"out"]
-o, --output-dir				| output directory [default:"."]
-O, --output-format				| bam, sam, fastq.gz, fastq [default:"fastq.gz"]
-1, --join-paired-ends          | merge R1 and R2 outputs in one file
-x, --compression-level         | speed compression: "1" - compress faster, <br />"9" - compress better [default:"6"; Integer]
-i, --append-id					| append to the defined template id [Format]
-I, --id						| overlap the default template id [Format]
-j, --jobs						| number of jobs [default:"1"; Integer]
-s, --seed						| set the seed of the base generator <br />[default:"time()"; Integer]
-c, --coverage					| fastq-file coverage [default:"8", Number]
-t, --sequencing-type			| single-end or paired-end reads <br />[default:"paired-end"]
-q, --quality-profile			| sequencing system profiles from quality <br />database [default:"poisson"]
-e, --sequencing-error			| sequencing error rate <br />[default:"0.005"; Number]
-m, --read-mean                 | read mean size for poisson [default:"100"; Integer]
-d, --read-stdd                 | read standard deviation size for poisson [default:"0"; Integer]
-M, --fragment-mean				| the fragment mean size for paired-end reads <br />[default:"300"; Integer]
-D, --fragment-stdd				| the fragment standard deviation size for paired-end reads <br />[default:"50"; Integer]
-a, --genomic-variation			| a list of structural variation entries from variation database. This option may be passed <br />multiple times [default:"none"]
-A, --genomic-variation-regex	| a list of perl-like regex to match structural <br />variation entries in variation database. <br />This option may be passed multiple times <br />[default:"none"]

**Some examples:**
1. These two commands, with equal effects, will produce two FASTq-files
(sequencing-type default is "paired-end"), both with a coverage of 20x
(coverage default is 8), and a plain-text reads-count file in a tab separated
fashion.
	```bash
	$ sandy genome --verbose --sequencing-type=paired-end --coverage=20 hg38.fa 2> sim.log
	```
or
	```bash
	$ sandy genome -v -t paired-end -c 20 hg38.fa 2> sim.log
	```
2. For reproducibility, the user can set an integer seed for the random raffles
with the `-s` option (seed default is environment `time()` value),
for example:
	```bash
	$ sandy genome -s 1220 my_fasta.fa
	```
3. To simulate reads from a ready registered database with a specific quality
profile other than default's one, type, for example:
	```bash
	$ sandy genome --quality-profile=hiseq_101 hg19.fa
	```
See the [quality profile](#howto) section to know how can a user register a
new profile on the database.
**Note:** If the user uses the option `-v`, by default, the log messages will
be directed to the standard error so, in the example above, it was redirected
to a file. Whithout the `-v` option, only errors messages will be printed.
4. The sequence identifier is the first and third line of a FASTq entry
beggining with a **@** token, for a read identifier, and a **+**,
for a quality identifier.
*SANDY* has the capacity to customize it, with a format string passed by
the user. This format is a combination of literal and escaped characters,
in a similar fashion used in **C** programming language's `printf`
function. For example, let's simulate a paired-end sequencing and put into it's
identifier the read length, read position and mate position:
	```bash
	$ sandy genome -s 123 --id="%i.%U read=%c:%t-%n mate=%c:%T-%N length=%r" hg38.fa
	```
In this case, results would	be:
	```bash
	$ sandy genome -s 123 --id="%i.%U read=%c:%t-%n mate=%c:%T-%N length=%r" hg38.fa
	==> Into R1
	@SR.1 read=chr6:979-880 mate=chr6:736-835 length=100
	...
	==> Into R2
	@SR.1 read=chr6:736-835 mate=chr6:979-880 length=100
	...
	```
5. To change the sequencing quality profile, use the `-q` option and a
string value (quality-profile default is "hiseq"):
	```bash
	$ sandy genome -q hiseq2 my_fasta_file.fa
	```
6. User can set the size of the reads with the `-r` option and an integer
number (reads-size default is 101):
	```bash
	$ sandy genome -r 151 my_fasta_file.fa
	```
7. User also can set the mean size of a fragment in a paired-end sequencing with
the `-m` option and an integer number (default is 300):
	```bash
	$ sandy genome -m 300 my_fasta_file.fa
	```
8. And, he can also set the standard deviation of the size of a fragment in
a paired-end sequencing with the `-D` option and an integer number
(default is 50):
	```bash
	$ sandy genome -D 30 my_fasta_file.fa
	```
9. The options above are the most frequently used ones for the `genome`
command, but many more can be found in the *SANDY's* documentation, with:
	```bash
	$ sandy genome --man
	```



##### Command `transcriptome` #####

Use it to generate simulated FASTq files from a given FASTA file,
according to an expression profile matrix file.
The `transcriptome` command sets these default options for a transcriptome
sequencing simulation as well:
* Choose the **Minus** strand;
* The number of reads is directly passed;
* The genes/transcripts are raffled following the expression matrix;

**Usage:**
```bash
$ sandy transcriptome [options] <fasta-file>
```
whose options' exaustive list can be consulted by `sandy transcriptome -h` or
even `sandy help transcriptome` commands.

Options							| Description
-------------------------------	| ----------------------------------------------
-h, --help						| brief help message
-u, --man						| full documentation
-v, --verbose					| print log messages
-p, --prefix					| prefix output [default:"out"]
-o, --output-dir				| output directory [default:"."]
-O, --output-format				| bam, sam, fastq.gz, fastq [default:"fastq.gz"]
-1, --join-paired-ends          | merge R1 and R2 outputs in one file
-x, --compression-level         | speed compression: "1" - compress faster, <br />"9" - compress better [default:"6"; Integer]
-i, --append-id					| append to the defined template id [Format]
-I, --id						| overlap the default template id [Format]
-j, --jobs						| number of jobs [default:"1"; Integer]
-s, --seed						| set the seed of the base generator <br />[default:"time()"; Integer]
-n, --number-of-reads           | set the number of reads [default:"1000000", Integer]
-t, --sequencing-type			| single-end or paired-end reads <br />[default:"paired-end"]
-q, --quality-profile			| sequencing system profiles from quality <br />database [default:"poisson"]
-e, --sequencing-error			| sequencing error rate <br />[default:"0.005"; Number]
-m, --read-mean                 | read mean size for poisson [default:"100"; Integer]
-d, --read-stdd                 | read standard deviation size for poisson [default:"0"; Integer]
-M, --fragment-mean				| the fragment mean size for paired-end reads <br />[default:"300"; Integer]
-D, --fragment-stdd				| the fragment standard deviation size for paired-end reads <br />[default:"50"; Integer]
-f, --expression-matrix         | an expression-matrix entry from database

**Some examples:**
1. The command:
	```bash
	$ sandy transcriptome --verbose --number-of-reads=1000000 --expression-matrix=brain_cortex gencode_pc_v26.fa.gz
	```
or, equivalently
	```bash
	$ sandy transcriptome -v -n 1000000 -f brain_cortex gencode_pc_v26.fa.gz
	```
will both generate a FASTq file with 1000000 reads from the *gencode_pc_v26.fa.gz*
file and a plain text file with the raw counts of the reads per gene,
according to the expression matrix provided by the *brain_cortex* entry
already registered in the database.
2. To demonstrate some other features, think about the sequencing error rate
that can be set between 0 and 1. By default, *SANDY* set this value to
0.005, which means 1 error every 200 bases.	To set it to another value,
try:
	```bash
	$ sandy transcriptome -f liver --sequencing-error=0.001 genome_pc_v26.fa.gz
	```
3. For reproducibility, the user can set the `seed` option and guarantee
the reliability of all the raffles in a later simulation.
	```bash
	$ sandy transcriptome -q hiseq_101 --seed=123 hg19.fa
	```
4. To have an idea of *SANDY's* plurality, look to how overwhelming the number
of choices could be:
	```bash
	$ sandy transcriptome \
		--expression-matrix=pancreas \
		--quality-profile=hiseq_101 \
		--sequencing-type=paired-end \
		--fragment-mean=350 \
		--fragment-stdd=100 \
		--prefix=pancreas_sim \
		--output-dir=sim_dir \
		--id="%i.%U read=%c:%t-%n mate=%c:%T-%N length=%r" \
		--verbose \
		--seed=123 \
		--jobs=30 \
		gencode_pc_v26.fa.gz
	```

**A note on paralelism:** To increase the processing speed, the simulation
can run in parallel, splitting the task among jobs. For example, type:
	```bash
	$ sandy custom -f testis -q hiseq_101 -v -i "length=%r" --jobs 15 gencode_lnc.fa.gz
	```
and *SANDY* will allocate 15 jobs. This feature works for the `genome` and
the `transcriptome` commands as well.



--------------------------------------------------------------------------------



#### Database Commands ####

##### Command `quality` #####

Use it to manage your quality profile database.
The user can add or remove your own expression profiles in the builtin database
and turn his simulations more realistic based on real experimental data.
Or he can even clean it up to restore the vendor's original entries state.
By default, *SANDY* uses a Poisson distribution when compiling the
quality entries, but like many other features, this behavior can be
overrided.

**Usage:**
```bash
$ sandy quality
$ sandy quality [options] 
$ sandy quality <sub-command>
```
whose options' exaustive list can be consulted by `sandy quality -h` or
even `sandy help quality` commands.

Options							| Description
------------------------------- | ----------------------------------------------
-h, --help						| brief help message
-u, --man						| full documentation
**Sub-Commands**				| 
add								| add a new quality profile to database
dump							| dump a quality-profle from database
remove							| remove an user quality profle from database
restore							| restore the database

**Some examples:**
1. To list the quality profiles already registered in the builtin database, you
can simply type:
	```bash
	$ sandy quality
	```
and all entries will be shown:
	```bash
	.------------------------------------------------------------------------------------------------------.
	| quality profile | mean  | stdd | error | type            | source            | provider | date       |
	+-----------------+-------+------+-------+-----------------+-------------------+----------+------------+
	| hiseq_101       |   101 |    0 | 0.001 | fragment        | 1000 genome       | vendor   | 2018-08-08 |
	| hiseq_150       |   150 |    0 | 0.001 | fragment        | SRA ID=SRR5805510 | vendor   | 2018-08-08 |
	| hiseq_51        |    51 |    0 | 0.001 | fragment        | SRA ID=SRR3185389 | vendor   | 2018-08-08 |
	| hiseq_76        |    76 |    0 | 0.001 | fragment        | SRA ID=SRR3355336 | vendor   | 2018-08-08 |
	| miseq_150       |   150 |    0 | 0.001 | fragment        | SRA ID=SRR6876696 | vendor   | 2018-08-08 |
	| miseq_301       |   301 |    0 | 0.001 | fragment        | SRA ID=SRR7089434 | vendor   | 2018-08-08 |
	| nextseq_51      |    51 |    0 | 0.001 | fragment        | SRA ID=SRR6131534 | vendor   | 2018-08-08 |
	| nextseq_85      |    85 |    0 | 0.001 | fragment        | SRA ID=SRR5445416 | vendor   | 2018-08-08 |
	| ont             | 15482 | 6195 |  0.25 | single-molecule | 10.1038/nbt.4060  | vendor   | 2018-08-08 |
	| pacbio          |  8817 | 3277 |  0.15 | single-molecule | 10.1038/nbt.2835  | vendor   | 2018-08-08 |
	| poisson         | -m    | -d   | -e    | both            | default           | vendor   | 2018-08-08 |
	'-----------------+-------+------+-------+-----------------+-------------------+----------+------------'
	```
2. To register a new probabilistic quality profile called, for example,
'my_profile.txt' to be used in the simulation of your FASTA-file.
User can use the `add` sub-command, typing:
	```bash
	$ sandy quality add -q 'my_quality_id' my_profile.txt
	```
This quality profile can be either a FASTq file or a plain text file in
a tab separated fashion (quality profile defaut density function is
*Poisson*).<br />
**Note:** Before the new entry can appear in the database's list, the new
profile needs to be validated, and if it can't, an error message will
be show. *SANDY* prevent's you before overwrite an existing entry.
3. To use a recently inserted quality profile over a given FASTA file
to simulate a transcriptomic data, use the `-q` option with the id you
registered:
	```bash
	$ sandy genome -q 'my_quality_id' my_fasta.fa
	```
4. Sometimes the user will need to update or delete some quality profile entry
(`my_profile.txt` for example) in the database. In this situation, he can
remove some actual entry and register a newer one, like this:
	```bash
	$ sandy quality remove 'my_quality_id'
	```
*SANDY* will refuse to remove any vendor's original entry from the database.
5. And, there could be times when users would wants to reset all the database to
its original state. It's a very simple command:
	```bash
	$ sandy quality restore
	```
Note that this is a dangerous command and *SANDY* will warn you about it
before make the restoration in fact.

**Note:** *SANDY* already comes with one quality profile based on the Poisson
probabilistic curve, as described by the literature
([illumina, 2018](https://www.illumina.com/content/dam/illumina-marketing/documents/products/technotes/technote_understanding_quality_scores.pdf)).


#### Command `expression` ####

The `expression` command is used to verify and update the expression matrix
database. In a transcriptome sequencing simulation, the user
must provide an expression matrix indexed into this database. *SANDY*
already comes with 52 different tissues from the GTEx project, but the
user has the freedom to include his own data as well, or even clean it up
to restore the vendor's original entries state.

**Usage:**
```bash
$ sandy expression
$ sandy expression [options] 
$ sandy expression <sub-command>
```
whose options' and sub-commands' exaustive list can be consulted by
`sandy expression -h` or even `sandy help expression` commands.

Options							| Description
------------------------------- | ----------------------------------------------
-h, --help						| brief help message
-u, --man						| full documentation
**Sub-Commands**				| 
add								| add a new expression-matrix to database
dump							| dump an expression-matrix from database
remove							| remove an user expression-matrix from database
restore							| restore the database

**Some examples:**
1. To list the expression matrixes already registered in the builtin database,
the user can simply type:
	```bash
	$ sandy expression
	```
and all registered entries will be shown:
	```bash
	.----------------------------------------------------------------------------------.
	| expression-matrix                   | source             | provider | date       |
	+-------------------------------------+--------------------+----------+------------+
	| adipose_subcutaneous                | Xena GTEx Kallisto | vendor   | 2018-08-08 |
	| adipose_visceral                    | Xena GTEx Kallisto | vendor   | 2018-08-08 |
	| adrenal_gland                       | Xena GTEx Kallisto | vendor   | 2018-08-08 |
	| artery_aorta                        | Xena GTEx Kallisto | vendor   | 2018-08-08 |
	| artery_coronary                     | Xena GTEx Kallisto | vendor   | 2018-08-08 |
	| artery_tibial                       | Xena GTEx Kallisto | vendor   | 2018-08-08 |
	| bladder                             | Xena GTEx Kallisto | vendor   | 2018-08-08 |
	| brain_amygdala                      | Xena GTEx Kallisto | vendor   | 2018-08-08 |
	| brain_anterior_cingulate_cortex     | Xena GTEx Kallisto | vendor   | 2018-08-08 |
	| brain_caudate                       | Xena GTEx Kallisto | vendor   | 2018-08-08 |
	| brain_cerebellar_hemisphere         | Xena GTEx Kallisto | vendor   | 2018-08-08 |
	| brain_cerebellum                    | Xena GTEx Kallisto | vendor   | 2018-08-08 |
	| brain_cortex                        | Xena GTEx Kallisto | vendor   | 2018-08-08 |
	| brain_frontal_cortex                | Xena GTEx Kallisto | vendor   | 2018-08-08 |
	| brain_hippocampus                   | Xena GTEx Kallisto | vendor   | 2018-08-08 |
	| brain_hypothalamus                  | Xena GTEx Kallisto | vendor   | 2018-08-08 |
	| brain_nucleus_accumbens             | Xena GTEx Kallisto | vendor   | 2018-08-08 |
	| brain_putamen                       | Xena GTEx Kallisto | vendor   | 2018-08-08 |
	| brain_spinal_cord                   | Xena GTEx Kallisto | vendor   | 2018-08-08 |
	| brain_substantia_nigra              | Xena GTEx Kallisto | vendor   | 2018-08-08 |
	| breast_mammary_tissue               | Xena GTEx Kallisto | vendor   | 2018-08-08 |
	| cells_ebv_transformed_lymphocytes   | Xena GTEx Kallisto | vendor   | 2018-08-08 |
	| cells_leukemia_cell_line            | Xena GTEx Kallisto | vendor   | 2018-08-08 |
	| cells_transformed_fibroblasts       | Xena GTEx Kallisto | vendor   | 2018-08-08 |
	| cervix_ectocervix                   | Xena GTEx Kallisto | vendor   | 2018-08-08 |
	| cervix_endocervix                   | Xena GTEx Kallisto | vendor   | 2018-08-08 |
	| colon_sigmoid                       | Xena GTEx Kallisto | vendor   | 2018-08-08 |
	| colon_transverse                    | Xena GTEx Kallisto | vendor   | 2018-08-08 |
	| esophagus_gastroesophageal_junction | Xena GTEx Kallisto | vendor   | 2018-08-08 |
	| esophagus_mucosa                    | Xena GTEx Kallisto | vendor   | 2018-08-08 |
	| esophagus_muscularis                | Xena GTEx Kallisto | vendor   | 2018-08-08 |
	| fallopian_tube                      | Xena GTEx Kallisto | vendor   | 2018-08-08 |
	| heart_atrial_appendage              | Xena GTEx Kallisto | vendor   | 2018-08-08 |
	| heart_left_ventricle                | Xena GTEx Kallisto | vendor   | 2018-08-08 |
	| kidney_cortex                       | Xena GTEx Kallisto | vendor   | 2018-08-08 |
	| liver                               | Xena GTEx Kallisto | vendor   | 2018-08-08 |
	| lung                                | Xena GTEx Kallisto | vendor   | 2018-08-08 |
	| minor_salivary_gland                | Xena GTEx Kallisto | vendor   | 2018-08-08 |
	| muscle_skeletal                     | Xena GTEx Kallisto | vendor   | 2018-08-08 |
	| nerve_tibial                        | Xena GTEx Kallisto | vendor   | 2018-08-08 |
	| ovary                               | Xena GTEx Kallisto | vendor   | 2018-08-08 |
	| pancreas                            | Xena GTEx Kallisto | vendor   | 2018-08-08 |
	| pituitary                           | Xena GTEx Kallisto | vendor   | 2018-08-08 |
	| prostate                            | Xena GTEx Kallisto | vendor   | 2018-08-08 |
	| skin_not_sun_exposed                | Xena GTEx Kallisto | vendor   | 2018-08-08 |
	| skin_sun_exposed                    | Xena GTEx Kallisto | vendor   | 2018-08-08 |
	| small_intestine_terminal_ileum      | Xena GTEx Kallisto | vendor   | 2018-08-08 |
	| spleen                              | Xena GTEx Kallisto | vendor   | 2018-08-08 |
	| stomach                             | Xena GTEx Kallisto | vendor   | 2018-08-08 |
	| testis                              | Xena GTEx Kallisto | vendor   | 2018-08-08 |
	| thyroid                             | Xena GTEx Kallisto | vendor   | 2018-08-08 |
	| uterus                              | Xena GTEx Kallisto | vendor   | 2018-08-08 |
	| vagina                              | Xena GTEx Kallisto | vendor   | 2018-08-08 |
	| whole_blood                         | Xena GTEx Kallisto | vendor   | 2018-08-08 |
	'-------------------------------------+--------------------+----------+------------'
	```
2. But, supose an user wants to register a new expression matrix file called
`my_expression.txt` to simulate the FASTA-file according to its experimentally
annnotated data. In this case, the sub-command `add` would solve your
problem:
	```bash
	$ sandy expression add -f 'my_expression_id' my_expression.txt
	```
Note that, before the new entry can appear in the database's list, the new
matrix file needs to be validated, and if it can't, an error message will
be show. *SANDY* prevents you to overwrite an existing entry.
3. So, to use the recently added expression matrix in a transcriptome
simulation, use the `-f` option on the `transcriptome` command:
	```bash
	$ sandy expression -f 'my_expression_id' my_fasta.fa
	```
4. Sometimes an user will need to update or delete some expression-matrix entry
('my_expression.txt', for example) in the database. In this situation, he can
remove the actual entry and register a newer one, like this:
	```bash
	$ sandy expression remove 'my_expression_id'
	```
*SANDY* will refuse to remove any vendor's original entry from the database.
5. Finally, there could be times when you would want to reset all the database to
its original state. It's a very simple command:
	```bash
	$ sandy expression restore
	```
Note that this is a dangerous command and *SANDY* will warn you about it
before make the restoration in fact.


#### Command `variation` ####

**Usage:**
```bash
$ sandy variation
$ sandy variation [options] 
$ sandy variation <sub-command>
```
whose options' and sub-commands' exaustive list can be consulted by
`sandy variation -h` or even `sandy help variation` commands.

Options							| Description
------------------------------- | ----------------------------------------------
-h, --help						| brief help message
-u, --man						| full documentation
**Sub-Commands:**				| 
add								| add a new structural variation to database
dump							| dump structural variation from database
remove							| remove an user structural variation from database
restore							| restore the database

**Some Examples:**
1. To show all variations entries in the database, type:
	```bash
	$ sandy variations
	```
and all entries for variations will be shown:
	```bash
	.-------------------------------------------------------------------------------------.
	| structural variation      | source                          | provider | date       |
	+---------------------------+---------------------------------+----------+------------+
	| NA12878_hg38_chr1         | IGSR - 1000 Genomes Phase 3     | vendor   | 2018-08-08 |
	| NA12878_hg38_chr10        | IGSR - 1000 Genomes Phase 3     | vendor   | 2018-08-08 |
	| NA12878_hg38_chr11        | IGSR - 1000 Genomes Phase 3     | vendor   | 2018-08-08 |
	| NA12878_hg38_chr12        | IGSR - 1000 Genomes Phase 3     | vendor   | 2018-08-08 |
	| NA12878_hg38_chr13        | IGSR - 1000 Genomes Phase 3     | vendor   | 2018-08-08 |
	| NA12878_hg38_chr14        | IGSR - 1000 Genomes Phase 3     | vendor   | 2018-08-08 |
	| NA12878_hg38_chr15        | IGSR - 1000 Genomes Phase 3     | vendor   | 2018-08-08 |
	| NA12878_hg38_chr16        | IGSR - 1000 Genomes Phase 3     | vendor   | 2018-08-08 |
	| NA12878_hg38_chr17        | IGSR - 1000 Genomes Phase 3     | vendor   | 2018-08-08 |
	| NA12878_hg38_chr18        | IGSR - 1000 Genomes Phase 3     | vendor   | 2018-08-08 |
	| NA12878_hg38_chr19        | IGSR - 1000 Genomes Phase 3     | vendor   | 2018-08-08 |
	| NA12878_hg38_chr2         | IGSR - 1000 Genomes Phase 3     | vendor   | 2018-08-08 |
	| NA12878_hg38_chr20        | IGSR - 1000 Genomes Phase 3     | vendor   | 2018-08-08 |
	| NA12878_hg38_chr21        | IGSR - 1000 Genomes Phase 3     | vendor   | 2018-08-08 |
	| NA12878_hg38_chr22        | IGSR - 1000 Genomes Phase 3     | vendor   | 2018-08-08 |
	| NA12878_hg38_chr3         | IGSR - 1000 Genomes Phase 3     | vendor   | 2018-08-08 |
	| NA12878_hg38_chr4         | IGSR - 1000 Genomes Phase 3     | vendor   | 2018-08-08 |
	| NA12878_hg38_chr5         | IGSR - 1000 Genomes Phase 3     | vendor   | 2018-08-08 |
	| NA12878_hg38_chr6         | IGSR - 1000 Genomes Phase 3     | vendor   | 2018-08-08 |
	| NA12878_hg38_chr7         | IGSR - 1000 Genomes Phase 3     | vendor   | 2018-08-08 |
	| NA12878_hg38_chr8         | IGSR - 1000 Genomes Phase 3     | vendor   | 2018-08-08 |
	| NA12878_hg38_chr9         | IGSR - 1000 Genomes Phase 3     | vendor   | 2018-08-08 |
	| NA12878_hg38_chrX         | IGSR - 1000 Genomes Phase 3     | vendor   | 2018-08-08 |
	| fusion_hg38_BCR-ABL1      | COSMIC - Gene Fusions in Cancer | vendor   | 2018-08-08 |
	| fusion_hg38_CCDC6-RET     | COSMIC - Gene Fusions in Cancer | vendor   | 2018-08-08 |
	| fusion_hg38_EML4-ALK      | COSMIC - Gene Fusions in Cancer | vendor   | 2018-08-08 |
	| fusion_hg38_EWSR1-ERG     | COSMIC - Gene Fusions in Cancer | vendor   | 2018-08-08 |
	| fusion_hg38_EWSR1-FLI1    | COSMIC - Gene Fusions in Cancer | vendor   | 2018-08-08 |
	| fusion_hg38_KIAA1549-BRAF | COSMIC - Gene Fusions in Cancer | vendor   | 2018-08-08 |
	| fusion_hg38_KMT2A-AFF1    | COSMIC - Gene Fusions in Cancer | vendor   | 2018-08-08 |
	| fusion_hg38_NCOA4-RET     | COSMIC - Gene Fusions in Cancer | vendor   | 2018-08-08 |
	| fusion_hg38_NPM1-ALK      | COSMIC - Gene Fusions in Cancer | vendor   | 2018-08-08 |
	| fusion_hg38_TMPRSS2-ERG   | COSMIC - Gene Fusions in Cancer | vendor   | 2018-08-08 |
	'---------------------------+---------------------------------+----------+------------'
	```
2. To increase the database with user's own data, use the `add` sub-command,
like this:
	```bash
	sandy variation add -a 'my_variations_id' my_vatiations.txt
	```
Note that, before the new entry can appear in the database's list, the new
variation's file needs to be validated, and if it can't, an error message will
be show. *SANDY* will prevent you to overwrite any existing entry, and *SANDY*
require these variations files to be in a GTF like format, specifying
coordinates on a reference genome with one variation per line (INDELs, SNVs and
gene fusions).
3. Now, to use the recently added variations specifications in a genomic
project, the user can use the `-a` option with the id registered for that file:
	```bash
	sandy genome -a 'my_variations_id' hg38.fa
	```
3. User can remove no-vendors entries from database as well:
	```bash
	sandy variation remove 'my_vatiations_id'
	```
Note that the user can't remove any vendor's entry.
4. Also, to reset all your variation entries to the original state (only with
the vendor's data), use the `restore` sub-command.
	```bash
	sandy variation restore
	```
5. Finally, when an user wants to simulate a reference genome with a high
coverage (ex. 50x) and insert some variations in it (maybe to abtain a positive
control for some other algorithm he's using), he can try this:
	```bash
	sandy genome -c 50 -a NA12878_hg38_chrX hg38.fa
	```
In this example, he has simulated reads for the whole genome, but the variations
are only in the X chromosome of the NA12878 individue in *SANDY's* database. An
even better way to insert variations to simulations is to use a
*regular expression* to search the entire database, like this:
	```bash
	sandy genome -c 50 -A NA12878* -a fusion_hg38_BCR-ABL1 hg38.fa
	```
This way, all entries that match NA12878 variations will be taken and,
additionally, a well studied gene fusion fusion_hg38_BCR-ABL1 introduced.

See? Notwithstanding the succint way that *SANDY's* commands are constructed, it
is possible to construct near real and complex genomes with on demand variations
in it with a few words!

--------------------------------------------------------------------------------



#### Miscellaneous Commands ####

##### Command `version` #####

*SANDY* project is made in a rolling release way, so the user can easily find
the version number he's using:
```bash
sandy version
```

##### Command `citation` #####

If *SANDY* was somehow useful, please cite its authors. With the `citation`
command, you can obtain a correct BibTeX entry and/or DOI number for the
version of *SANDY* you're using:
```bash
sandy citation
```

--------------------------------------------------------------------------------



#### Help Commands ####

**Usage:**
To get a simple general help, user can type any of these commands:
```bash
$ sandy --help
```
or for short
```bash
$ sandy -h
```
or simply call it without any arguments.
```bash
$ sandy
```

But, if a more comprehensive explanation is needed, invoke *SANDY's*
manual:
```bash
$ sandy --man
```
or for short
```bash
$ sandy -M
```

For help about specific commands, its options and inputs, type:
```bash
$ sandy help <command>
```
or
```bash
$ sandy <command> -h
```

We always exhort users to get help by consulting *SANDY's* builtin documentations
with `man sandy` or `info sandy` commands in their terminals.

--------------------------------------------------------------------------------



### Docker Usage ###

The user can run many instances of *SANDY* in a scalable way by pulling its
Docker [image](https://hub.docker.com/r/galantelab/sandy) from Docker Hub in
a way aforementioned in the [Installation](#or-get-sandy-in-a-docker-image) section.

Here, we describe how to port all the commands shown above to be used in a Docker
cointainer in a very straightforward way. For example, given the command:
```bash
$ sandy help genome
```

All user has to do is substitute the word `sandy` by
`docker run --rm -ti [options] galantelab/sandy`, like:
```bash
$ docker run --rm -ti [options] galantelab/sandy help genome
```
And the *options* are about the folders which the user wants to map inside the
container.

Let's see another example, supose the user is in a directory like
`host_path/folder1/` containing the file `gencode_pc_v26.fai.gz` on which he is
trying to use the command bellow:
```bash
$ sandy transcriptome \
		--expression-matrix=pancreas \
		--quality-profile=hiseq_101 \
		--sequencing-type=paired-end \
		--fragment-mean=350 \
		--fragment-stdd=100 \
		--prefix=pancreas_sim \
		--output-dir=sim_dir \
		--id="%i.%U read=%c:%t-%n mate=%c:%T-%N length=%r" \
		--verbose \
		--seed=123 \
		--jobs=30 \
		gencode_pc_v26.fa.gz
```

So, to adapt such a command to a Docker usage looking up to the correct path of
the directories containing the data, only substitute the first line with:
```bash
$ docker run --rm -ti -v /ABSOLUTE/host_path/folder1:/ABSOLUTE/container_path/folder2 galantelab/sandy transcriptome \
		--expression-matrix=pancreas \
		--quality-profile=hiseq_101 \
		--sequencing-type=paired-end \
		--fragment-mean=350 \
		--fragment-stdd=100 \
		--prefix=pancreas_sim \
		--output-dir=sim_dir \
		--id="%i.%U read=%c:%t-%n mate=%c:%T-%N length=%r" \
		--verbose \
		--seed=123 \
		--jobs=30 \
		gencode_pc_v26.fa.gz
```

The `-v /ABSOLUTE/host_path:/ABSOLUTE/container_path/folder1` option maps the
directory `folder1` (adding its absolute path) in the host to the `folder2`, in
the container, at `/ABSOLUTE/container_path/` directory. Obviously those paths
could be something like `/home/user/dataset/`, we are just highlighting the
importance of using the absolute paths here, otherwise it won't work correctly.
Additionally, the `-v` option can be used repeatedly in the same command, as
many times as the number of directories the data needs.

See [Docker documentation](https://docs.docker.com/) for more information about
options and commands for Docker.



--------------------------------------------------------------------------------
[Back to Top](#contents-at-a-glance) | [Back to main page](../README.md)