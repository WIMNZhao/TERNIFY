# TERNIFY: Efficient Sampling of PROTAC-Induced Ternary Complexes

[![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](LICENSE)
[![Version](https://img.shields.io/badge/version-2025.06.03-green.svg)]()

TERNIFY is a high-performance C++ implementation for efficient sampling and prediction of PROTAC-induced ternary complex structures. It uses advanced molecular modeling techniques to predict how PROTAC molecules bring together E3 ligases and target proteins to form productive ternary complexes.

## Overview

TERNIFY employs a sophisticated computational approach that combines:

- **Molecular docking and alignment** of PROTAC components
- **Force field-based energy calculations** using MMFF94s
- **Monte Carlo sampling** with simulated annealing
- **Conformational clustering** based on RMSD analysis
- **Multi-threading** for enhanced performance

## Features

- ✅ **High-Performance Computing**: Multi-threaded implementation for fast sampling
- ✅ **Advanced Clustering**: RMSD-based conformational clustering to remove redundancy
- ✅ **Flexible Parameterization**: Customizable energy thresholds and sampling parameters
- ✅ **Multiple Output Formats**: SDF for molecules, PDB for protein complexes
- ✅ **Verbose Mode**: Detailed progress tracking and energy analysis

## Installation

### Prerequisites

- **Conda/Miniconda** or **Anaconda**
- **CMake** (version 3.30 or higher)
- **C++17** compatible compiler (GCC 7+ or Clang 5+)

### Step-by-Step Installation

#### 1. Create and Setup Conda Environment

```bash
# Create an environment for ternify
conda create -n ternify python=3.12
conda activate ternify

```

#### 2. Install Dependencies

```bash
# Install RDKit and development libraries
conda install rdkit=2025.03.2 -c conda-forge
conda install librdkit-dev=2025.03.2 -c conda-forge

# Install Boost libraries
conda install libboost-devel=1.86.0 -c conda-forge
#conda install libboost-headers=1.86.0 -c conda-forge
conda install libboost-python-devel=1.86.0 -c conda-forge  # may be optional

# Install Eigen, NLopt, CMake
conda install eigen=3.4.0 -c conda-forge
conda install nlopt=2.10.0 -c conda-forge
conda install -c conda-forge gcc_linux-64=12 gxx_linux-64=12 cmake
conda update -c conda-forge cmake   # ensures cmake >= 3.3, default is 4.0.2
```

#### 3. Configure CMakeLists.txt

Before compilation, modify the `CMakeLists.txt` file to set the correct Conda environment path if necessary:

```cmake
# Set RDKit, NLopt library paths (adjust paths to your conda environment)
set(CONDA_ENV "/opt/anaconda3/envs/ternify/")
link_directories(${CONDA_ENV}/lib)
set(RDKIT_DIR ${CONDA_ENV})
set(NLOPT_DIR ${CONDA_ENV})
```

**Note**: Replace `/opt/anaconda3/envs/ternify/` with your actual conda environment path. You can find it using:

```bash
conda info --envs
```

Set library path if necessary for compiling and running ternify afterwards:
```bash
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
```
**Note**: Replace `$CONDA_PREFIX` with your actual conda environment path.

#### 4. Compile TERNIFY

```
# modify the rdkit, rdkit lib, nlopt lib paths in `CMakeLists.txt`

# set rdkit, nlopt lib !!!!
set(CONDA_ENV "/opt/anaconda3/envs/ternify/")
link_directories(${CONDA_ENV}/lib)
set(RDKIT_DIR ${CONDA_ENV})
set(NLOPT_DIR ${CONDA_ENV})

# compile ternify
conda activate ternify


mkdir build
cd build
cmake ..
make -j 8


# Optional: Install globally
sudo cp ternify /usr/local/bin
# =======
# Run ternify
```

### Troubleshooting Installation

**CMake can't find RDKit:**

- Verify your conda environment path in CMakeLists.txt
- Ensure all RDKit packages are properly installed
- Check that the environment is activated: `conda activate ternify`

**Compilation errors:**

- Make sure you're using a C++17 compatible compiler
- Verify all dependencies are installed with correct versions
- Try cleaning build directory: `rm -rf build && mkdir build`

**Library linking issues:**

- Double-check the `link_directories()` path in CMakeLists.txt
- Ensure `librdkit-dev` is installed alongside `rdkit`

## Usage

### Basic Usage

```bash
ternify -p tcs.inp
```

### Input File Format

Create a parameter file (e.g., `tcs.inp`) with the following format:

```ini
# Input files
PROTACs: protac.sdf
Warhead_anchor: e3.sdf
Warhead_flex: poi.sdf
Protein_anchor: e3.pdb
Protein_flex: poi.pdb

# Output files
Output_protac: TC_protac.sdf
Output_protein: TC_protein.pdb

# Sampling parameters
N_ini: 10000          # Initial sampling size
N_search: 900         # Number of solutions to optimize
N_keep: 100           # Number of final structures to keep
N_processes: 8        # Number of CPU threads

# Clustering and thresholds
RMSD_cutoff: 1.0      # RMSD threshold for clustering (Å)


# Score_only mode (align to poi.sdf and e3.sdf,optimize H,then scoring )
Score_only: 0

#  Optimization only model
Local_only: 0

# Interface definition (x_min, x_max, y_min, y_max, z_min, z_max)
Interface: -15.0, 15.0, -15.0, 15.0, -15.0, 15.0

# Debug information
Verbose: 1            # Verbose level (<1 quiet, >0 moderate, >1 detailed)

```

### Parameter Description

| Parameter       | Description                          | Default | Range      |
| --------------- | ------------------------------------ | ------- | ---------- |
| `N_ini`       | Initial conformational sampling size | 10000   | 1000-50000 |
| `N_search`    | Solutions selected for optimization  | 900     | 100-5000   |
| `N_keep`      | Final structures to output           | 100     | 10-1000    |
| `N_processes` | CPU threads for parallel processing  | 1       | 1-64       |
| `RMSD_cutoff` | Clustering threshold in Ångströms  | 1.0     | 0.5-5.0    |
| `Verbose`     | Output verbosity level               | 0       | 0-2     |

### Verbose Modes

- **Verbose: 0** - Minimal output, only essential information
- **Verbose: 1** - Progress bars, energy breakdown information
- **Verbose: 2** - Structural analysis, Detailed force field parameters and debugging information

## Output Files

### PROTAC Structures (`TC_protac.sdf`)

- Contains clustered and ranked PROTAC conformations
- Sorted by binding energy (lowest first)
- Each structure represents a unique ternary complex geometry

### Protein Complex (`TC_protein.pdb`)

- Multi-model PDB file with aligned protein conformations
- Each MODEL corresponds to a PROTAC conformation
- REMARK records contain energy scores

## Algorithm Overview

1. **Initialization**: Load and prepare molecular structures
2. **Alignment**: Align PROTAC to anchor and flexible warheads
3. **Sampling**: Generate diverse conformations using MC sampling
4. **Optimization**: Refine structures with simulated annealing
5. **Clustering**: Remove redundant conformations using RMSD
6. **Ranking**: Sort by total binding energy

## Energy Components

TERNIFY calculates the following energy terms:

- **E_intra**: Intramolecular PROTAC energy (VdW + torsional)
- **E_anchor**: PROTAC-anchor protein interaction
- **E_flex**: PROTAC-flexible protein interaction
- **E_pp**: Protein-protein interaction energy

$$
E_{Total Energy} = E_{intra} + E_{anchor} + E_{flex} + E_{pp}
$$

## Performance Tips

- **CPU Cores**: Set `N_processes` to your CPU core count for optimal performance
- **Memory Usage**: Large `N_ini` values require more RAM (~1GB per 10k structures)
- **Clustering**: Smaller `RMSD_cutoff` values produce more diverse results but take longer
- **Sampling**: Balance between `N_ini` (diversity) and `N_search` (quality)

## Troubleshooting

### Common Issues

**"Failed to setup MMFF properties"**

- Ensure input molecules have proper atom types
- Check for missing hydrogens or unusual atoms

**"No match found for warhead"**

- Verify PROTAC contains both anchor and flexible warhead substructures
- Check molecular connectivity and atom mapping

**Memory errors**

- Reduce `N_ini` or `N_search` parameters
- Increase system RAM or swap space

### Debug Mode

Run with verbose mode for detailed diagnostics:

```bash
# Create debug input file with Verbose: 2
echo "Verbose: 1" >> debug.inp
ternify debug.inp
```

## Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/improvement`)
3. Commit changes (`git commit -am 'Add new feature'`)
4. Push to branch (`git push origin feature/improvement`)
5. Create a Pull Request

## License

This project is licensed under the Apache License 2.0 - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- RDKit development team for cheminformatics toolkit
- Eigen library for linear algebra operations
- Contributors and beta testers

## Contact

- **Hongtao Zhao, PhD** - Principal Developer
- **Ximing Xu, PhD** - C++ Implementation Lead xuximing@ouc.edu.cn

For questions and support, please open an issue on GitHub or contact the development team.
