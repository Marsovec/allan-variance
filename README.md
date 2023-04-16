# Allan variance calculator

This repository tracks development on the command-line utility to calculate Allan variance from raw (byte) data, written in Python. A script to plot the calculated variances is also provided.

## Terminology

### Chunk

A chunk is a part of the input file that is calculated upon, usually in units of bytes (as chunks are unpacked into bits only afterwards). For practical purposes, each chunk can thus be considered as a separate input file. The input file is separated into chunks because it would otherwise prove impractical or otherwise unwanted to assess the entire file as a whole (like the [RAM usage problems mentioned later](#note-on-ram-usage)).

### Cluster

A cluster is a part of the chunk of data that is being assessed. Cluster sizes are denoted by the letter $T$ further below and are all in units of bits. A chunk is broken up in this way so that we can calculate per-cluster averages which are defined below (also note that it can be broken down in an [overlapping](#overlapping-allan-variance) or a [non-overlapping](#non-overlapping-allan-variance) fashion). A cluster in position $k$ is usually denoted by $\Omega_k$, though this is rarely used as we are only interested in the averages.

### Cluster average

A cluster average is the mean value of all elements in a specific cluster. This quantity is denoted by $\overline{\Omega}_k$, where $k$ is the index of the cluster. The calculation is done by way of a sum of all elements divided by the size of the cluster. The notation $\overline{\Omega}_k(T)$ is also used, where $T$ is the cluster size.

## Background

The calculator is built using formulae specified in [this article](https://www.allaboutcircuits.com/technical-articles/intro-to-allan-variance-analysis-non-overlapping-and-overlapping-allan-variance/). Essentially, a file is read as raw data and then converted to a sequence of $N$ elements with a possible value of 0 or 1 (each byte is unpacked into bits). The sequence is then split into $K$ clusters, with each cluster having a length of $T$. Here we separate between calculating overlapping and non-overlapping Allan variance.

### Overlapping Allan variance

The overlapping variant is calculated with clusters that are selected by shifting the selection interval by 1 for each selection of size $T$. This separates the sequence into $N-T+1$ clusters. The variance is then calculated 

$$\sigma^2(T) = \frac{1}{2(N-2T+1)}\sum_{k=1}^{N-2T+1}(\overline{\Omega}_{k+T}(T)-\overline{\Omega}_k(T))^2$$

The sum is over the per-cluster averages $\overline{\Omega}_k$, where $k$ is the index of the cluster.

### Non-overlapping Allan variance

The non-overlapping variant is calculated with clusters that each begin where the previous one has ended in the sequence. The sequence is thus separated into $N/T$ clusters. The variance is calculated using

$$\sigma^2(T) = \frac{1}{2(K-1)}\sum_{k=1}^{K-1}(\overline{\Omega}_{k+1}(T)-\overline{\Omega}_k(T))^2$$

### Comparison between both variants

Non-overlapping variance is calculated much faster than overlapping simply because the former has far less clusters to consider ($N/T$ vs. $N-T+1$). Also, based on experience, the variants don't differ much in practice, so this repository assumes the non-overlapping variant as default.

## Note on RAM usage

I have given special attention to reducing RAM usage but have found no stable way to keep it down with large chunk sizes (the `-b` argument). I advise caution with chunk sizes above the default of 1MB, the largest I have tried and succeeded is 10MB, but with significant memory usage.

For this reason the experimental `-s` option exists, which is supposed to sum variances from each chunk with the previous in order to have an accurate and cheap alternative to large chunks being loaded into memory, but it does not perfectly recreate the accurate calculations directly with larger chunks, so its usage should not be considered seriously.

I can imagine ways to reduce memory usage with really carefully loading data into memory and calculating the cluster sums concurrently so as to not keep them in memory, but the current implementation, with its relatively small chunk sizes that are still feasible, is already sufficient to show trends, and that should be good enough in practice.

## Requirements

Python packages requirements are listen in the [requirements.txt](./requirements.txt) file. They can be installed using `pip`

```sh
pip install -r requirements.txt
```

You can of course also do this inside `venv`

```sh
python -m venv venv_dir
source venv_dir/bin/activate
```

## Usage

Both [allan.py](./allan.py) and [plot.py](./plot.py) have basic help pages that can be accessed respectively

```sh
python allan.py --help
python plot.py --help
```

### [allan.py](./allan.py)

The script requires at least one provided `filename` argument which contains the path to a data file. Each provided file is then parsed and processed according to the following command-line arguments (which are all optional and have, I believe, sane defaults):

* `-b`, `--blocks` - length of one chunk in MB (1MB = one million bytes), any float value is allowed, defaults to 1
* `-a`, `--avars` - whether to calculate overlapping or non-overlapping Allan variance, options are `overlap`,`nonoverlap`,`both`, defaults to `nonoverlap`
* `-T`, `--tmin` - smallest cluster size exponent to calculate with (T=2^tmin), any integer value is allowed, defaults to 0 (T=1)
* `-s`, `--sum` - whether to sum each chunk with the previous in an effort to approximate larger-chunk calculations (EXPERIMENTAL), true if specified, defaults to false
* `-c`, `--count` - number of chunks to consider, any integer value is allowed, defaults to all
* `-o`, `--offset` - offset for where to start reading data in the file in MB, any float value is allowed, defaults to 0
* `-d`, `--outdir` - path to the output directory for the generated .dat files, any string value is allowed, defaults to no output

If output path directory is specified, `outdir` is populated with a `_INFO.txt` file that contains the arguments used and .dat files with variances per chunk. Each .dat file contains two columns: the variance in regards to the selected cluster size T. The cluster sizes for which calculations are done are selected from an exponential function $T=2^n$ where n ranges from 0 to $\log_2(B/2)$ (B is the chunk size as specified by `-b`).

### [plot.py](./plot.py)

This script presents the information calculated by [allan.py](./allan.py) and stored in `outdir` as .dat files. It requires at least one `filename` argument with the path to a .dat file. Output images are saved with a timestamp in the same directory as the .dat file, or in the current working directory if `-a` is selected (read below). The following optional command-line arguments are accepted:

* `-a`, `--append` - represent all data on one plot, true if specified, defaults to false (separate files)
* `-l`, `--legend` - do not draw legend, do not draw if specified, draw by default
