      ============================================================
      *                         COMPETE                          *
      * Competitive Occupancy: Multi-factor Profile Evaluation   *
      *                of a Thermodynamic Ensemble               *
      *         COMPETE is licensed from Duke University.        *
      *    Copyright (c) 2009-2020 by Alexander J. Hartemink.    *
      *                   All rights reserved.                   *
      ============================================================
 

## Installation

1. Clone this repository:

    ```bash
    git clone git@github.com:jianlingzhong/COMPETE.git
    ```

2. Enter the `libconfig` directory, compile, and install  `libconfig`:
   
    ```bash
    cd COMPETE/COMPETE/libconfig/libconfig-1.1
    ./configure --disable-cxx --disable-shared --prefix=`pwd`/../libconfig-1.1_inst/
    make install
    ```
   
3. Return to the `COMPETE` directory and make it:

    ```bash
    cd ../.. && make compete
    ```

## Requirements

`COMPETE` has been developed for POSIX compatible operating systems, such as
Linux, Mac OS X, and various flavors of BSD. It should run in Windows within an
installation of Cygwin, but this is neither tested nor supported. `COMPETE`
includes tools that require the Ruby programming language to be installed,
version 1.8 or higher, but the binary `COMPETE` executable can be run manually
without it.

`COMPETE` can require a significant amount of memory to run, depending on the
model (i.e. how many, and which, DNA binding factors are included; nucleosomes
are by far the biggest) and the length of the DNA sequence being analyzed. This
can range from a few megabytes to multiple gigabytes, and some trial and error
will be required to determine what is feasible on a particular machine. Running
top and watching the size of memory allocated by `COMPETE` at the beginning of a
run is a good way of determining its memory usage for a given configuration.

`COMPETE` will use two CPUs (or at least, two threads) in machines with multiple
CPUs, but this is not required. 

## Run `COMPETE`

### The Easy Route: Running `COMPETE` Via Wrapper Scripts

`COMPETE` is most easily run via a front-end wrapper script named `wrapper.rb` (inside `ruby_scripts` directory).
`wrapper.rb` takes a human-readable `YAML` input file that defines the input
parameterization as described in the following (real) example input file:

```yaml
# Roman numeral indicating which chromosome to run on
chr: IV

# coordinates on above chromosome on which to run; indexed starting at 1
range: [ 740741, 761628 ]

# Boolean indication of whether to include nucleosomes in this model
nuc: true

# array of DBFs (other than nucs) to include in this model
motifs: [ fkh2, mcm1 ]

# hash of concentrations of included DBFs
conc: { nuc: 40.0, unbound: 1.0, fkh2: 0.01, mcm1: 0.01 }

# inverse temperature parameter; generally left at 1.0, meaning no effect
inverse_temp: 1.0

# intermediate files will be named based on this
output_basename: fig3b

# if given, this will specifically define the output filename; 
# otherwise it will be generated from the output_basename and the given DBFs and concentrations
output_filename: fig3b.txt
```

Given this YAML file, `COMPETE` may be run as follows:

```bash
./wrapper.rb < wrapper_input.yaml
```

Note that `wrapper.rb` reads the `YAML` input from standard input.  This allows more
adventurous users to modify elements of it on the fly if so desired.  An example
of this will be shown below when discussing how to set concentrations as
multiples of K<sub>d</sub>, as described in [Wasson and Hartemink, 2009](http://www.ncbi.nlm.nih.gov/pubmed/19720867).

It should be noted that out of the box, `COMPETE` is configured to run on the
yeast genome that all of the analyses in the original `COMPETE` publication were
run on.  Configuring `COMPETE` to run on other sequences, or for other organisms,
is described below.

### The Hardcore Route: Running `COMPETE` Manually

Note: Although `COMPETE` can certainly be run manually, and doing so is described
in this section, the `wrapper.rb` front-end has become our method of choice, and
will be the focus of examples going forward.

Running `COMPETE` manually involves a few steps: creation of the model to include
your preferred DBFs, creation of a file containing the sequence filename and
coordinates, and actual execution of the `COMPETE` binary executable:

1. Creation of the Model File

    Creating the model file is accomplished with the `construct_model_from_motifs.rb`
    Ruby script.  Running the program with no arguments will give usage information.
    Generally, the script takes a space-delimited list of DBF IDs and a `-n` to
    include nucleosomes.  The model is output to standard out, so it must be
    redirected into a file for storage.
    
    For example, to create a model including `GAL4`, `FKH2`, and nucleosomes, and save
    it into a file named model.cfg, do the following:
    
    ```bash
    ./construct_model_from_motifs.rb gal4 fkh2 -n > model.cfg
    ```
    
    To create the same model, but without nucleosomes:
    
    ```bash
    ./construct_model_from_motifs.rb gal4 fkh2 > model.cfg
    ```
    
    Valid DBF IDs can be seen by looking at the filenames in the `pbm` directory.

2. Creation of the Sequence Filenames File

    The sequence filenames file exists solely for future support of runs across
    multiple sequences.  This may be implemented if called for, or may be phased
    out.  Regardless, as of now, only one sequence can be given at a time.
    
    The file is plain-text and has a very simple format, as follows:
    
    ```txt
    chr/II.txt 267712 289412
    ```
    
    This indicates usage of the `chr/II.txt` file, positions `267712` to `289412`
    inclusively, indexed starting at `1`.  The input file needs to be produced by a
    run of the `convert_seq.rb` Ruby script, which is described in the section "Using
    Your Own Sequence Files".

3. Executing the `COMPETE` Binary Executable

    The actual `COMPETE` binary, named `compete`, takes a collection of command line
    arguments, as follows:
    
    ```txt
    usage: compete [options] model_file seq_file
      -n  nucleosome_concentration (float)
      -m  motif_concentrations (comma delimited string of floats)
      -N  motif_labels (comma delimited string of strings, for output file column headers)
      -u  unbound_concentration (float)
      -t  inverse_temperature (float)
      -s  output only probabilities of starting each DBF per postion
    ```
    
    This usage can be printed at any time by running `compete` with no arguments, or
    by passing it `-h` or `--help`.
    
    All arguments are optional. If concentrations are unspecified, nucleosomes will
    default to a concentration of `1` and all other DBFs will default to `0.01`.  Any
    example run from the command line is:
    
    ```bash
    ./compete -n 1.0 -m 0.01,0.1,0.01 -u 1.0 -t 1.0 model.cfg seq_filenames.txt > output.txt
    ```
    
    In this case, `nucleosome` concentration is set to `1.0`, the other three `DBF`
    (motif) concentrations are set to `0.01`, `0.1`, and `0.01` respectively, the `unbound`
    concentration is set to `1.0`, and the system temperature is set to `1.0`.  Like
    other programs included in the `COMPETE` package, the `compete` binary outputs to
    standard out and must be redirected into a file; in this case, `output.txt`.

## Using Your Own Sequence Files

`COMPETE` expects sequence files to be in a simple format that is easily
obtainable using the included `convert_seq.rb` Ruby script.  This scipt converts
`FASTA` formatted files into files with `ACGT` converted to `ASCII characters` 0, 1,
2, and 3.  An example is as follows:

```bash
./convert_seq.rb -H chr/I.fasta > chr/I.txt
```

In this case, the `-H` is necessary to strip the `FASTA` header line.

## Using Your Own Motifs

`COMPETE` will use any motifs in the format of the files in the `pbm` directory.
These files contain one motif each, and are of a simple format:

```txt
line 1: Comment describing the contents of the file
line 2: "A" followed by tab delimited probability of A in all N positions of motif
line 3, 4, 5: Same as line 2, but for C, G, and T respectively
```

Currently the only way to point `COMPETE` at motifs in another directory is to
edit the `construct_model_from_motifs.rb` Ruby script, but it should be easy.
Simply look for line 118 of that file, which reads:

```ruby
motifs = parse_pbm_motif_directory("pbm")
```

Change the `pbm` in quotes to any other path, and it will look in that directory
instead.

*Note:* `COMPETE` can support `TAMO` formatted motifs as well, but this is being
phased out in support of the one-motif-per-file format of the `.pwm` files,
described above.


## Plotting the Results

`COMPETE` outputs a tab delimited file with columns of data corresponding to the
DBFs included in the model file used, and rows numbering the same as the length
of input sequence.  A simple way to plot one of these for sequence positions `1`
to `500`, in `R`, would be:

```R
d = read.table("foo.txt", header = T)
plot(d[1:500, 1], type="l")
```

The data file format should be easy enough to parse with any software you prefer
for plotting; this is merely a simple example.


## Example Files

Several example files are included in the `examples` directory.  The input
files required to recreate the examples in the publication figures are included
in the `examples/wrapper_inputs` directory.

To build the examples in figures 2 and 5 in [Wasson and Hartemink, 2009](http://www.ncbi.nlm.nih.gov/pubmed/19720867), the symlink from `construct_model_from_motifs.rb` to `construct_model_from_pbm_motifs.rb` must be
pointed instead at `construct_model_from_macisaac_motifs.rb`.  This will use the
[MacIsaac et al., 2006](https://github.com/jianlingzhong/COMPETE/edit/master/COMPETE/README.md) motifs in the included `TAMO` file.


## More Information

For further information and support, please visit:

[http://www.cs.duke.edu/~amink/software/compete](http://www.cs.duke.edu/~amink/software/compete)

If you have any questions, issues, or find any bugs, please do not hesitate to
contact the developers (contact information available on the website).
