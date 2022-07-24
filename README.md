# Minimizer-based de Bruijn Graph De Novo Assembler

Master's thesis implementation under the supervision of prof. Mile Šikić - based on [[0]](#0).

The assembler uses a sparse representation of input reads in minimizer space - determined by
the density parameter. By using a predetermined set of minimizers, the initial input reads are 
transformed into a compact representation from which a de Bruijn graph is built and simplified.
A GFA graph is generated as output, which can be further transformed and analyzed.

Best results are achieved using highly accurate reads with high coverage such as PacBio HiFi.

The assembler supports parallel graph construction and simplification, as well as trio-binning.
Rolling hashes are used for minimizer detection and graph construction which means the
running time is independent of the specified minimizer length `l` and graph node length `k`.
This allows the use of arbitrary values for these parameters without incurring penalties.

## Building

Dependencies:
 * Intel TBB
 * Catch2 - optional, required only for tests

Other dependencies are specified as submodules in `vendor/` or automatically fetched with CMake.

```
git clone --recurse-submodules https://github.com/yatsukha/cpp-mdbg
cd cpp-mdbg
mkdir release_build && cd release_build
cmake .. -DCMAKE_BUILD_TYPE=Release && make -j 4
```

## Running

From the build directory:
```
./mdbg input_files output_prefix
```

For help just run `./mdbg`:
```
C++ minimizer-based de Bruijn graph assembler.
Usage:
  mdbg [OPTION...] input.fastx output.gfa

  -t, --threads arg       Maximum number of concurrent threads. NOTE:
                          Default of 0 means max concurrency. (default: 0)
  -k, --kmers arg         Length of window of minimizers to use for the de
                          Bruijn graph. (default: 33)
  -l, --letters arg       Length of the minimizers. (default: 14)
  -d, --density arg       Density of the universe minimizers. (default:
                          0.005)
  -a, --analysis          Exit after outputting minimizer statistics for
                          given reads. (default: 0)
      --dry-run           Dry run, do not write. (default: 0)
  -s, --sequences         Output sequences contained within minimizers in
                          output GFA. (default: 0)
      --trio-binning arg  Format: K:T:reads0.fa:reads1.fa
                          Enables trio binning using K length kmers for
                          counting; discards kmers with frequency below T.
                          K must be <= 32. (default: "")
```

### Running with a malloc proxy

In some cases, better performance can be achieved by using an allocator designed for
multithreaded usage such as the TBB malloc proxy:
```
LD_PRELOAD=/usr/lib/libtbbmalloc_proxy.so.2 ./mdbg ...
```

Depending on how you installed TBB it might be in different location such as
`/opt/intel/oneapi/`, `~/anaconda3/pkgs/`, etc.

### Running tests

Assuming Catch2 was found during the build process, from the build directory run:
```
./test
```

## Parameter tuning

Parameters can be fine tuned using `-a`, which displays statistics and exits the program.
For minimizer density `d = 0.02` and minimizer length `l = 14` the program could output
something like this:
```
minimizers per read:
  median:             398
  90th    percentile: 370
  99th    percentile: 344
  99.9th  percentile: 325
  99.99th percentile: 312
  min:                302
```

Based on this we would pick the length of the graph nodes `k` to be at most `302`.

In general there, isn't a single right way to pick parameters. A high density will
mean the reads are not compacted as much which prolongs the assembly but may
result in more accurate graphs. A lower density will assemble the graph faster
but there may be misassemblies if the gaps between minimizers are too large.
Depending on the error rate you may get away with a very low density.

For the minimizer length the default value of `l = 14` is a good starting point,
while the graph node length `k` should be picked based on statistics. `k` in the range
`[30, 200]` usually yields good results.

## README references
<a id="0">[0]</a> 
Barış Ekim, Bonnie Berger, and Rayan Chikhi.
Minimizer-space de Bruijn graphs: Whole-genome assembly of long reads in minutes on a personal computer.
Cell Systems (2021).
