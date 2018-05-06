## Welcome to Sandy simulator! ##

![logo.png](img/sandy_logo.png)

**Let's make a simulation today???**

If you're looking for a bioinformatics tool that provides a simple engine to generate
single-end/paired-end reads from a given FASTA file or expression matrix file,
then *Sandy* is your choice!



### Introduction ###

Many next-generation sequencing (NGS) analyses rely on hypothetical
models and principles that are not precisely satisfied in practice. Simulated
data, which provides positive controls would be a perfect way to overcome
these difficulties. Nevertheless, most of NGS simulators are extremely
complex to use, they do not cover all kinds of the desired features needed by
the users, and (some) are very slow to run in a standard computer. Here, we
present SANDY, a straightforward, easy to use, fast, complete set of tools to
generate synthetic next-generation sequencing reads. SANDY simulates
whole genome sequencing, whole exome sequencing, RNAseq reads and it
presents several features to the users manipulate the data. Sandy can be
used therefore for benchmarking results of a variety of pipelines in the
genomics or trancriptomics.

Now, project *Sandy* is in it's 0.15 version and has earned enough maturity to
simulate some realistic features, among these:
* Simulate reads from genomic FASTA-files.
* Simulate reads from transcriptomic data, based on expression matrix files.
* Ready included databases for *quality profiles* and *expression matrixes*.
* Import and record your own *expression matrixes* and *quality profiles* to
simulate future data.



### Contents at a Glance ###

1. [Introduction](#introduction)
2. [Installation](#installation)
3. [Usage and Option summary](usage.md#usage-and-option-summary)
	1. [General](usage.md#general-syntax)
    2. [Command `genome`, its options and examples](usage.md#the-genome-command)
    3. [Command `transcriptome`, its options and examples](usage.md#the-transcriptome-command)
    4. [Command `custom`, its options and examples](usage.md#the-custom-command)
    5. [Command `quality`, its options and examples](usage.md#the-quality-command)
    6. [Command `expression`, its options and examples](usage.md#the-expression-command)
    7. [Command `help`, its options and examples](usage.md#the-help-command)
4. [A case study example](case.md#a-case-study-example)
5. [Aknowledgements](#aknowledgements)
6. [Author](#author)
7. [Copyright and License](#copyright-and-license)



### Installation ###

You can install it by two different approaches.

1. If you already use `perl` and perl modules `cpanm`, the solution comes
in one line:
	```bash
		$ cpanm App::Sandy
	```

2. If you only have `perl`, as a last resort, you can manually install *Sandy*
through the command line by downloading the [tarball](https://github.com/galantelab/sandy/archive/sandy-master.tar.gz)
from GitHub, decompressing it and then building it, like this:
	```bash
		$ wget https://github.com/galantelab/sandy/archive/sandy-master.tar.gz
		$ tar xzvf sandy-master.tar.gz
		$ cd sandy
		$ perl Makefile.PL
		$ make && make test
	```
	Then install it properly with:
	```bash
		$ make install
	```

For more details, see the INSTALL file on *Sandy's* GitHub [repository](https://github.com/galantelab/sandy).



### Aknowledgements ###


### Author ###

Thiago L. A. Miller
[<tmiller@mochsl.org.br>](tmiller@mochsl.org.br)



### Copyrirht and License ###

This software is Copyright (c) 2018 by Teaching and Research Institute from Sírio-Libanês Hospital.
This is free software, licensed under:

`The GNU General Public License, Version 3, June 2007`
