# Analysis scripts for polymer simulations of chromosome organisation.

## Cite

Amith Z. Abdulla, Marco Di Stefano, Cristian Micheletti and Daniel Jost, Topological constraints...
Numerical simulation code availabe https://github.com/physical-biology-of-chromatin/LatticePoly/tree/ExtrusionPoly

This repository contains a collection of scripts for analysing 3D polymer simulation trajectories of chromosome territories. The pipeline quantifies chromosome intermingling through mixing entropy, alpha-shape analysis, computes contact maps with trans-contact ratios (TR), and generates violin-plot visualisations. The scripts are designed to run on HPC clusters (SLURM) and follow a numbered, sequential workflow.

---

## Table of Contents

1. [Overview](#overview)
2. [Requirements](#requirements)
3. [Input Data Format](#input-data-format)
4. [Pipeline Overview](#pipeline-overview)
5. [Script Reference](#script-reference)
   - [a) Mixing Entropy — Local Flory-Huggins Mixing Entropy](#mixing-entropy--local-flory-huggins-mixing-entropy)
   - [b) Alpha-Shape Volume and Particle Counting](#step-1--alpha-shape-volume-and-particle-counting)
   - [c) Alpha-Shape Visualisation](#step-2a--alpha-shape-visualisation)
   - [d) Contact Map and Trans-Contact Ratio (Per-Replica)](#step-2b--contact-map-and-trans-contact-ratio-per-replica)
   - [e) Contact Map and Trans-Contact Ratio (All Replicas)](#step-3a--contact-map-and-trans-contact-ratio-all-replicas)
   - [f) Trans-Ratio Violin Plots](#step-3b--trans-ratio-violin-plots)
   - [g) Cis/Trans Contact Count Violin Plots](#step-4--cistrans-contact-count-violin-plots)
   - [C Programs — Contact Map Computation Engine](#c-programs--contact-map-computation-engine)
   - [Utility — Knot-State Reader](#utility--knot-state-reader)
6. [Output Files](#output-files)

---

## Overview

The simulation system consists of **three polymer chains** (A, B, C), each comprising **2,000 beads**, confined inside a sphere of radius 10 lattice units. The analysis pipeline addresses three main questions:

- **Degree of mixing**: How well-mixed are the three polymer species locally? (Flory-Huggins mixing entropy, `T3mixing_Entropy.py`.)
- **Territorial organisation**: How much volume does each chain occupy, and how many beads of one chain reside inside the territory of another? (Alpha-shape analysis, Steps 1–2a.)
- **Contact patterns**: What is the balance between cis-chromosome and trans-chromosome contacts, and how does it evolve over time? (Contact-map analysis, Steps 2b–4.)

## Requirements

### Software

| Tool | Version / Notes |
|------|----------------|
| **Python 3** | With `numpy`, `scipy` (for mixing entropy and knot-state utility). A `vtkReader` module is also required for reading simulation output (see note below). |
| **R** | >= 4.0, with packages: `alphashape3d`, `rgl`, `ggplot2`, `dplyr`, `ggpubr`, `ggrepel`, `scales`, `RColorBrewer` |
| **GCC** | Any standard C compiler (for building the contact-map programs) |
| **SLURM** | Scripts include `#SBATCH` directives; adapt for other schedulers as needed |
| **Bash** | >= 4.0 |

> **Note on `vtkReader`**: `T3mixing_Entropy.py` imports a `vtkReader` module to read simulation trajectory data. This module is not included in this repository and must be provided separately. It should expose a class whose iterator yields objects with a `polyPos` attribute containing an (N, 3) NumPy array of monomer coordinates.

### Hardware

Most analysis steps are single-threaded. Memory requirements are moderate (the SLURM headers request 15 GB), driven primarily by the contact-matrix allocation inside the C programs.

## Input Data Format

### Simulation Output (for mixing entropy)

VTK-format trajectory files, read via the `vtkReader` module. Each frame must provide the 3D coordinates of all monomers for each of the three polymer chains. The chains are identified by name (e.g. `"poly0"`, `"poly1"`, `"poly2"`) and loaded separately.

### Simulation Snapshots (for alpha-shape analysis)

Plain-text files named `time_<T>_data`, where each line contains the x, y, z coordinates of a single bead, space-separated:

```
x1 y1 z1
x2 y2 z2
...
```

The file is assumed to contain **6,000 lines** in total — the first 2,000 lines correspond to chain A, lines 2,001–4,000 to chain B, and lines 4,001–6,000 to chain C. Coordinates are in lattice units.

### Simulation Snapshots (for contact-map analysis)

Plain-text files named `time_<T>_pos.res` (or `.read_data`), each containing **N_particles** lines of x, y, z coordinates. The number of particles equals `beads_per_chain * number_of_chains` (e.g. 6,000 for a 3-chain system).

### Auxiliary Files (generated automatically)

- **`_chromosomes`**: Two-column file mapping each bead index to a chromosome identifier. Created automatically by the shell wrappers.
- **`_compartments`** (optional): Two-column file mapping each bead index to a compartment label (0 = telomere, 1 = A, 2 = B). Used only by `compute_contact_map_TR_and_CS`.
- **`DATA_FILES_INPUT.txt`**: A list of absolute paths to snapshot files to be processed. Created automatically by the shell wrappers.


---

## Script Reference

### Mixing Entropy — Flory-Huggins Mixing Entropy

#### `T3mixing_Entropy.py`

**Purpose**: Computes the normalised Flory-Huggins mixing entropy for each monomer in a three-chain polymer system. For every monomer, a local neighbourhood is defined as all monomers within a given radius. The local volume fractions of the three species are used to compute a per-monomer mixing entropy, which is then normalised by the maximum possible entropy for three species (ln 3). The result quantifies how well-mixed the chains are locally: a value of 1.0 indicates perfect mixing, while values near 0 indicate segregation.

**Language**: Python 3

**Dependencies**: `numpy`, `scipy` (`scipy.spatial.cKDTree`), `vtkReader` (external module, not included)

**Command-line arguments**:

| Argument | Position | Description |
|----------|----------|-------------|
| outputDir | `$1` | Path to the directory containing the simulation output files |
| chrom1 | `$2` | Identifier for the first polymer chain (e.g. `"poly0"`) |
| chrom2 | `$3` | Identifier for the second polymer chain (e.g. `"poly1"`) |
| chrom3 | `$4` | Identifier for the third polymer chain (e.g. `"poly2"`) |
| initFrame | `$5` | Frame index to analyse (integer) |
| radius | `$6` | Radius of the local neighbourhood for entropy calculation (float, in simulation units) |

**Algorithm**:

1. Reads the 3D coordinates of all three chains for the specified frame using `vtkReader`.
2. Concatenates all monomer positions and assigns a species label (0, 1, or 2) to each.
3. Builds a KD-tree over all monomer positions.
4. For each monomer with more than 2 neighbours within the specified radius:
   - Computes the local volume fraction φ_i of each species.
   - Computes the local mixing entropy: −Σ_i (φ_i ln φ_i).
5. Normalises each per-monomer entropy by ln(3) to obtain values in the range [0, 1].

**Key parameters**:

| Parameter | Value | Description |
|-----------|-------|-------------|
| `n_polymers` | 3 | Number of distinct polymer species |
| `min_neighbors` | 2 | Minimum number of neighbours required to compute entropy for a monomer |

**Input files**: VTK trajectory files in `outputDir`, read via `vtkReader`.

**Output file**:

| File | Description |
|------|-------------|
| `r<initFrame>_localEntropy_r<radius>.res` | One value per line: the normalised mixing entropy for each monomer that had sufficient neighbours. Written to `outputDir`. |

**Idempotency**: The script checks whether the output file already exists and aborts if it does.

**Usage**:

```bash
python3 T3mixing_Entropy.py /path/to/output poly0 poly1 poly2 100 1.5
```

This computes the normalised local mixing entropy at frame 100 using a neighbourhood radius of 1.5 simulation units.

---

### Alpha-Shape Volume and Particle Counting

#### `01_compute_ashape3D_volumeNparticles.R`

**Purpose**: Computes 3D alpha-shapes for each of the three polymer chains and determines (i) the alpha-shape volume as a function of the alpha parameter, and (ii) the number of beads from the other two chains that fall inside each chain's alpha-shape.

**Language**: R

**Dependencies**: `alphashape3d`, `rgl`

**Input files** (expected in the current working directory):

| File | Description |
|------|-------------|
| `_xyzAall` | x, y, z coordinates of chain A (2,000 lines) |
| `_xyzBall` | x, y, z coordinates of chain B (2,000 lines) |
| `_xyzCall` | x, y, z coordinates of chain C (2,000 lines) |

**Key parameters (hardcoded)**:

- `factor = 0.5 * sqrt(2)` — base alpha step size
- `alphas = seq(factor, factor * (20/factor) + 1, factor)` — range of alpha values evaluated

**Output files** (written to the current working directory):

| File | Description |
|------|-------------|
| `Avolume.tab` | Alpha values and corresponding alpha-shape volumes for chain A |
| `Bvolume.tab` | Same for chain B |
| `Cvolume.tab` | Same for chain C |
| `Bparticles_in_A.tab` | Number of chain-B beads inside chain A's alpha-shape, per alpha |
| `Cparticles_in_A.tab` | Number of chain-C beads inside chain A's alpha-shape, per alpha |
| `Aparticles_in_B.tab` | Number of chain-A beads inside chain B's alpha-shape, per alpha |
| (etc.) | All six pairwise combinations |

**Output format**: Tab-separated, two columns per file — `alpha_value` and either `volume` or `particle_count`.

**Usage**: This script is not run directly. It is called by the shell wrapper `01_compute_ashape3D_volumeNparticles.sh`.

---

#### `01_compute_ashape3D_volumeNparticles.sh` / `.cmd`

**Purpose**: SLURM-compatible shell wrapper that iterates over simulation replicates and timesteps, extracts per-chain coordinate files from the raw data, and calls the R script above.

**Language**: Bash

**Parameters**:

| Argument | Position | Description |
|----------|----------|-------------|
| Condition name | `$1` | 

**Key variables (edit before running)**:

| Variable | Description | Default |
|----------|-------------|---------|
| `inputDir` | Path to the raw simulation data | 
| `workingDir` | Scratch directory for analysis output | 

**SLURM resources**: 1 core, 15 GB RAM, 4-day time limit, `Lake` partition.

**What it does**:

1. Creates a working directory for the given condition.
2. For each replicate directory, and for each `time_<T>_data` file within it:
   - Extracts lines 1–2,000 into `_xyzAall`, lines 2,001–4,000 into `_xyzBall`, and lines 4,001–6,000 into `_xyzCall`.
   - Runs `01_compute_ashape3D_volumeNparticles.R`.
3. Skips replicate/timestep combinations that have already been analysed (checks for `Avolume.tab`).

**Example submission**:

```bash
sbatch 01_compute_ashape3D_volumeNparticles.cmd Topo_Feb2025
```

---

### Alpha-Shape Visualisation

#### `02_chain_alphaShape_plots.R`

**Purpose**: Generates violin plots (with overlaid boxplots) showing the distribution of alpha-shape volumes and particle-in-hetero-shape fractions across simulation replicates.

**Language**: R

**Dependencies**: `ggplot2`, `dplyr`, `ggpubr`, `ggrepel`, `scales`, `RColorBrewer`

**Input files** (expected in the current working directory):

| File | Description |
|------|-------------|
| `chain_alphaShape_volume.tab` | Aggregated alpha-shape volumes (columns: alpha, timestep, volume) |
| `chain_particlesInHeteroShape.tab` | Aggregated particle-in-hetero-shape counts (columns: alpha, timestep, count) |

**Key parameters (hardcoded)**:

| Parameter | Value | Description |
|-----------|-------|-------------|
| `radius` | 10.0 | Radius of the confining sphere (lattice units) |
| `totBeads` | 2000 | Number of beads per chain |

Volumes are normalised by the total confining sphere volume; particle counts are normalised by `totBeads`.

**Output files**: PDF violin plots:

- `chain_alphaShapeVol_vs_alpha_time_<T>.pdf` — volume distribution vs. alpha at timestep T
- `chain_alphaShapeVol_vs_time_tMax_20_alpha_<A>.pdf` — volume distribution vs. time at alpha A
- `chain_particlesInHeteroShape_vs_alpha_time_<T>.pdf` — particle fraction vs. alpha at timestep T
- `chain_particlesInHeteroShape_vs_time_alpha_<A>.pdf` — particle fraction vs. time at alpha A
- `alphaShapeVol_vs_time.tab` — summary table of plateau values

**Usage**: Called by `02_chain_alphaShape_plots.cmd`.

---

#### `02_chain_alphaShape_plots.cmd`

**Purpose**: SLURM wrapper that aggregates per-replicate alpha-shape results into a single table and then calls the plotting R script.

**Parameters**:

| Argument | Position | Description |
|----------|----------|-------------|
| Condition name | `$1` | 

**What it does**:

1. Changes to the analysis root directory.
2. For each matching condition directory:
   - If `chain_alphaShape_volume.tab` does not yet exist, aggregates individual `?volume.tab` and `?particles_in_*.tab` files across all replicates and timesteps.
   - Runs `02_chain_alphaShape_plots.R`.

**SLURM resources**: 1 core, 15 GB RAM, 4-day time limit.

**Example submission**:

```bash
sbatch 02_chain_alphaShape_plots.cmd Topo_Feb2025_intermingled
```

---

### Contact Map and Trans-Contact Ratio (Per-Replica)

#### `02_compute_contact_map_PSMN_TR.sh`

**Purpose**: Core shell script that computes a contact map and the trans-contact ratio for a **single replicate** over a given time window, using the `compute_contact_map_TR` C program.

**Parameters**:

| Argument | Position | Description |
|----------|----------|-------------|
| Simulation directory | `$1` | 
| Replica directory | `$2` | 
| t_min | `$3` | Start timestep of the analysis window |
| t_max | `$4` | End timestep of the analysis window |

**Key hardcoded parameters**:

| Parameter | Value | Description |
|-----------|-------|-------------|
| `rc` | 1.5 | Contact distance cutoff (lattice units) |
| `nCopies` | 3 | Number of polymer chains |
| `b` | 1.0 | Monomer diameter (lattice units) |
| `modelResolution` | 1 | Map resolution in beads |

**Output files** (per replicate/time window):

| File | Description |
|------|-------------|
| `rc_<rc>nm_res_<res>bead_from_<tmin>_to_<tmax>_every_<dt>_<replica>_TR.tab` | Per-bin trans-contact ratio |
| `rc_<rc>nm_..._TR_perChrom.tab` | Per-chromosome average cis/trans fractions |
| `rc_<rc>nm_..._<replica>.tab` | Raw contact matrix |

---

#### `02_compute_contact_map_PSMN_TR.cmd`

**Purpose**: SLURM submission script that calls `02_compute_contact_map_PSMN_TR.sh` for each replicate and each time window in a given simulation condition.

**What it does**: Iterates over replicate directories and time windows (t_max from 10 to 100 in steps of 10), checking that the corresponding snapshot file exists before launching the computation.

**Example submission**:

```bash
sbatch 02_compute_contact_map_PSMN_TR.cmd
```

Note: The simulation condition is currently hardcoded to `TopoNov` in this script. Edit the `for simDir in ...` line to target different conditions.

---

#### `02_compute_contact_map_genotoul_TR_and_CS.sh` / `.cmd`

**Purpose**: Variant of the contact-map pipeline for the **Genotoul cluster**. In addition to the trans-contact ratio, this version also computes the **compartment strength** (CS), using the `compute_contact_map_TR_and_CS` program.

**Parameters for `.cmd`**:

| Argument | Position | Description |
|----------|----------|-------------|
| Number of copies | `$1` | Ploidy / number of chain copies |
| Replica index | `$2` | Which replica to process |

**Parameters for `.sh`**:

| Argument | Position | Description |
|----------|----------|-------------|
| Condition directory | `$1` | Simulation condition subdirectory |
| Volume fraction | `$2` | Phi value for selecting the run directory |
| Min replicates | `$3` | Minimum number of replicates required |
| t_min | `$4` | Start timestep |
| t_max | `$5` | End timestep |
| t_delta | `$6` | Timestep interval |
| Tag | `$7` | Quality tag (`OK` or `tmp`) |
| Replica | `$8` | Replica index |

**Additional hardcoded parameters**:

| Parameter | Value | Description |
|-----------|-------|-------------|
| `modelResolution` | 5,000 bp | Model resolution |
| `mapResolution` | 10,000 bp | Map resolution |
| `rc` | 150 nm | Contact distance cutoff |

**Additional output**: Per-compartment (A/B) trans-ratio tables and compartment-strength tables (`*_CSHomoOverAll_perChrom.tab`).

---

### Contact Map and Trans-Contact Ratio (All Replicas)

#### `03_compute_contact_map_PSMN_TR.sh` / `.cmd`

**Purpose**: A variant of Step 2b that pools **all replicates** at each individual timestep to compute a single contact map per timepoint, rather than per-replica maps over a time window. Designed for intermingled configurations.

**Parameters for `.cmd`**:

| Argument | Position | Description |
|----------|----------|-------------|
| Condition name | `$1` | Simulation condition (e.g. `Topo_Feb2025`). The suffix `_intermingled` is appended automatically. |

**Parameters for `.sh`**:

| Argument | Position | Description |
|----------|----------|-------------|
| Simulation directory | `$1` | Full name with `_intermingled` suffix |
| (unused) | `$2` | Placeholder (set to `all`) |
| t_min | `$3` | Start timestep |
| t_max | `$4` | End timestep |
| Input directory | `$5` | Path to raw simulation data |

**Output files**: Same structure as Step 2b, but named with `_at_<tmax>_` instead of `_from_<tmin>_to_<tmax>_`.

**Example submission**:

```bash
sbatch 03_compute_contact_map_PSMN_TR.cmd Topo_Feb2025
```

---

### Trans-Ratio Violin Plots

#### `03_make_TR_violinPlots.R`

**Purpose**: Generates a violin plot of the average trans-contact ratio per chromosome across different simulation conditions or timepoints.

**Language**: R

**Input file** (expected in the current working directory):

| File | Columns |
|------|---------|
| `AvgTransRatio_perChrom.tab` | `sample` (tab) `chrom` (tab) `avg_cis_ratio` (tab) `avg_trans_ratio` |

**Output**: `AvgTransRatio_perChrom.pdf`

**Y-axis range**: 0.01 to 0.75 (hardcoded).

---

#### `03_make_TR_violinPlots.cmd`

**Purpose**: SLURM wrapper that aggregates per-chromosome TR values from the contact-map output files, constructs the input table, and calls the R plotting script.

**What it does**:

1. Creates an `AvgTransRatio_perChrom.tab` header.
2. For each condition, collects `*_TR_perChrom.tab` files across timepoints.
3. Runs the R script.
4. Moves the output PDF/table to `TR_distributions/`.

**Key variables to edit**:

- The `for condition in ...` loop lists the conditions to process.
- The `for t in ...` loop defines the time range (default: 10 to 100 in steps of 10).

**Example submission**:

```bash
sbatch 03_make_TR_violinPlots.cmd
```

---

### Cis/Trans Contact Count Violin Plots

#### `04_make_cisCounts_violinPlots.R`

**Purpose**: Generates a violin plot of raw cis-chromosome contact counts over time.

**Language**: R

**Command-line argument**:

| Argument | Position | Description |
|----------|----------|-------------|
| y_max | `$1` | Upper limit of the y-axis (integer) |

**Input file**:

| File | Columns |
|------|---------|
| `cisCounts.tab` | `time` (tab) `totCounts` (tab) `cisCounts` (tab) `transCounts` |

**Output**: `cisCounts.pdf`

**Usage**:

```bash
Rscript 04_make_cisCounts_violinPlots.R 1860
```

---

#### `04_make_transCounts_violinPlots.R`

**Purpose**: Generates a violin plot of raw trans-chromosome contact counts over time. Identical structure to `04_make_cisCounts_violinPlots.R`, but plots the `transCounts` column.

**Command-line argument**:

| Argument | Position | Description |
|----------|----------|-------------|
| y_max | `$1` | Upper limit of the y-axis (integer) |

**Input file**: Same `cisCounts.tab` as above (reads the `transCounts` column).

**Output**: `transCounts.pdf`

**Usage**:

```bash
Rscript 04_make_transCounts_violinPlots.R 1800
```

---

#### `04_make_cisCounts_violinPlots.cmd`

**Purpose**: SLURM wrapper that aggregates TR output into a cis/trans counts table and calls both violin-plot R scripts (cis and trans).

**Parameters**: None (conditions and time ranges are hardcoded).

**Key variables to edit**:

- `tMax` — maximum timestep to include (default: 20)
- `for condition in ...` — list of simulation conditions to process

**What it does**:

1. For each condition, collects `*_TR.tab` files and extracts total, cis, and trans contact counts.
2. Runs `04_make_cisCounts_violinPlots.R` with y_max = 1860.
3. Runs `04_make_transCounts_violinPlots.R` with y_max = 1800.
4. Moves output PDFs to `cisCounts_distributions/` and `transCounts_distributions/`.

**Example submission**:

```bash
sbatch 04_make_cisCounts_violinPlots.cmd
```

---

### C Programs — Contact Map Computation Engine

#### `compute_contact_map_TR.c`

**Purpose**: Reads a set of 3D coordinate snapshots, computes an all-vs-all contact map based on a distance cutoff, and outputs the **trans-contact ratio** (fraction of trans-chromosome contacts) per genomic bin.

**Compilation**:

```bash
gcc -O2 -o compute_contact_map_TR compute_contact_map_TR.c -lm
```

Note: This program includes `My_Memory.c` and `My_IO.c` via `#include` directives. You must either place these files in the expected path or update the `#include` paths in the source.

**Command-line options**:

| Flag | Type | Description |
|------|------|-------------|
| `-r` | int | Resolution of the contact map in monomers (beads per bin) |
| `-p` | int | Total number of particles in the system |
| `-d` | float | Contact distance cutoff between particles (in the same units as coordinates) |
| `-b` | float | Monomer diameter (used to normalise the distance) |
| `-k` | int | Ploidy — number of chromosome copies in the system |
| `-i` | int | Inter-copy contact mode: `0` = consider all contacts, `1` = intra-copy only |

**Required input files** (in the current working directory):

| File | Description |
|------|-------------|
| `DATA_FILES_INPUT.txt` | One absolute path per line to a coordinate snapshot file |
| `_chromosomes` | Two-column mapping: `bead_index  chromosome_id` |

**Output files**:

| File | Description |
|------|-------------|
| `TR_per_bin.txt` | Per-bin trans-contact ratio. Columns: `bin_chrom_start_end_comp`, `cis_contacts`, `cis_count`, `trans_contacts`, `trans_count`, `total_contacts`, `total_count`, `cis_fraction`, `trans_fraction` |
| `contacts.tab` | Raw contact matrix (sparse format) |
| `output.log` | Computation log with timing information |

**Example**:

```bash
./compute_contact_map_TR -r 1 -p 6000 -d 1.5 -b 1.0 -k 1 -i 0
```

---

#### `compute_contact_map_TR_and_CS.c`

**Purpose**: Extended version of `compute_contact_map_TR` that additionally computes the **compartment strength** (CS) — the observed-over-expected ratio of same-compartment contacts. The CS computation is present in the source but is currently commented out; the program can be reactivated for compartment-strength analysis by uncommenting the relevant blocks.

**Compilation and usage**: Identical to `compute_contact_map_TR.c` above.

**Additional optional input file**:

| File | Description |
|------|-------------|
| `_compartments` | Two-column mapping: `bead_index  compartment_id` (0 = telomere, 1 = A, 2 = B) |

**Additional output** (when CS code is enabled):

| File | Description |
|------|-------------|
| `CS_per_bin_<chrom>.txt` | Per-bin compartment strength for each chromosome |

---

#### `My_IO.c`

**Purpose**: Utility library providing convenience wrappers for file I/O operations (`open_r` for reading, `open_w` for writing, `open_a` for appending). Used by both C programs via `#include`.

---

#### `My_Memory.c`

**Purpose**: Utility library for dynamic memory allocation and deallocation of multi-dimensional arrays (vectors, matrices, 3D and 4D tensors) for `int`, `float`, `double`, and `char` types. Used by both C programs via `#include`.

---

### Utility — Knot-State Reader

#### `read_npy.py`

**Purpose**: Reads `.npy` files containing knot-state classification data from polymer simulations and prints a summary of knot-state frequencies.

**Language**: Python 3

**Dependencies**: `numpy`

**Input**: Hardcoded `.npy` filenames:

- `ExtSing_data_matrix.npy`
- `Initial_NoneSing_data_matrix.npy`
- `NoneSing_data_matrix.npy`
- `TopoSing_data_matrix.npy`

Each `.npy` file is expected to contain a NumPy array of records, where the 4th element (index 3) of each record is the knot-state label.

**Output**: Prints the total number of records and a dictionary of knot-state counts to standard output.

**Usage**:

```bash
python3 read_npy.py
```

---

## Output Files

The pipeline produces several categories of output:

| Category | File pattern | Description |
|----------|-------------|-------------|
| Mixing entropy | `r<frame>_localEntropy_r<radius>.res` | Per-monomer normalised mixing entropy |
| Alpha-shape volumes | `?volume.tab` | Per-chain volume vs. alpha |
| Particle intermingling | `?particles_in_?.tab` | Bead counts inside other chains' territories |
| Contact matrices | `rc_*nm_res_*bead_*.tab` | Raw contact maps |
| Trans-contact ratios | `*_TR.tab`, `*_TR_perChrom.tab` | Per-bin and per-chromosome TR |
| Compartment strength | `*_CSHomoOverAll_perChrom.tab` | Per-chromosome CS (when enabled) |
| Violin plots | `*.pdf` | Visualisations |

---


## Notes

- **Path configuration**: Several scripts contain hardcoded paths (e.g. `inputDir`, `workingDir`, paths to the C binaries). Before running, update these to match your local directory structure.
- **Idempotency**: Most scripts check whether output files already exist before re-running a computation. To force recomputation, delete the corresponding output files.
- **Adapting for other schedulers**: The `#SBATCH` directives can be replaced with equivalent directives for PBS, SGE, or other schedulers.
- **Compartment strength**: The CS calculation in `compute_contact_map_TR_and_CS.c` is currently commented out. To enable it, uncomment the relevant code blocks and provide a `_compartments` file.
