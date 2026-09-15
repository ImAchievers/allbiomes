# allbiomes

Searches for the smallest ocean-bounded continent that contains all 42 overworld
biomes (MC_26_3). Built on the cubiomes fork included under `lib/`.

## Build

Make sure to set `START_SEED` and `NUM_THREADS` in the defines before compiling.

MSYS2/MinGW on Windows, or Linux:

    make

## Run

    ./allbiomes.exe                          Run the search (resumes from checkpoint)
    ./allbiomes.exe single <seed>            Trace one seed through the pipeline
    ./allbiomes.exe single <seed> <x> <z>    Measure the continent at (x, z)

Results are saved to `allbiomes_results.txt` as `[seed] | [biomes] | [size] |
[missing/ALL]`. Progress saves to `allbiomes_checkpoint.txt` every 30s and on
Ctrl-C; a rerun resumes from that point if the file exists.

## Tuning

The main settings are the `#define`s at the top of `allbiomes.c`: `START_SEED`,
`NUM_THREADS`, `MAX_AREA` (Size Threshold), and `MIN_REPORT` (Minimum biome count
for a result to be logged).

## Searching for larger or smaller continents

Just change `MAX_AREA`. The ocean-ring radius that confirms a continent is an
island (`RING_R`) is set automatically as `0.6 * sqrt(MAX_AREA)`, so it
always sits just larger than what you're searching for. 
Leave `RING_MIN` at 8. Past ~9,000,000 also increase `ISO_W` (1536 -> 2048) so a big
isolated continent is not misread as land-bridged.