# AGENTS.md

This file provides guidance to Codex (Codex.ai/code) when working with code in this repository.

## Project Overview

This is a **5G Downlink System Level Simulator** (Tbps_DLSLS) that performs wireless channel simulations for 5G networks. The simulator models base stations (BS), mobile stations (MS), wireless channels, and scheduling algorithms to evaluate network performance under various configurations.

## Build and Run

### Building
```bash
# Create build directory if it doesn't exist
mkdir -p build
cd build

# Configure and build
cmake ..
make

# Executable will be in the root directory
cd ..
```

### Running Simulations
```bash
# Run with a configuration file
./Tbps_DLSLS Tbps_DU_B.cfg

# Other available configurations:
./Tbps_DLSLS Tbps_InH_C.cfg
./Tbps_DLSLS Tbps_Rural_B.cfg
```

## Dependencies

- **Eigen3** (version 3.3+): Matrix and linear algebra operations
- **OpenMP**: Parallel processing for channel computations
- **FFTW3**: Fast Fourier Transform for channel modeling
- **C++11**: Standard required

The project uses Homebrew-installed libraries on macOS (libomp at `/opt/homebrew/Cellar/libomp/14.0.0/include`).

## Code Architecture

### Main Simulation Flow ([src/main.cpp](src/main.cpp))

The simulation follows this sequence:
1. **Initialization Phase**: Set simulation parameters, initialize system, allocate memory
2. **Drop Loop**: For each drop (snapshot of network topology):
   - Generate user/BS locations
   - Initialize channel models
   - Generate Large-Scale Parameters (LSP) and Small-Scale Parameters (SSP)
   - Configure links between BSs and MSs
   - **Time Loop**: For each time slot (`t`):
     - Update channel coefficients
     - Perform scheduling (assign resources to users)
     - Compute downlink reception and throughput
     - Collect statistics
3. **Output Results**: Write performance metrics to files

### Core Components

#### Base Station (BS)
- **Header**: [h/BS.h](h/BS.h)
- **Implementation**: [src/BS.cpp](src/BS.cpp)
- Represents a base station with location, antenna elements, and sector information
- Antenna configuration: `d_tx[sector][M][N][P][Mg][Ng]` for 3D antenna arrays

#### Mobile Station (MS)
- **Header**: [h/MS.h](h/MS.h)
- **Implementation**: [src/MS.cpp](src/MS.cpp)
- Models user equipment with methods for:
  - Channel updates and quantization
  - CQI (Channel Quality Indicator) computation
  - SINR (Signal-to-Interference-plus-Noise Ratio) calculation
  - Link-to-System (LLS) mapping (EESM, RBIR, MIB methodologies)
  - Downlink reception with MMSE receiver

#### Channel Model
- **Header**: [h/channel.h](h/channel.h)
- **Implementation**: [src/channel.cpp](src/channel.cpp), [src/channel_update.cpp](src/channel_update.cpp)
- Implements 3GPP channel models for different scenarios (InH, UMa, UMi, RMa)
- Handles:
  - Pathloss calculation
  - LOS/NLOS probability
  - Large-Scale Parameters (LSP): shadow fading, delay spread, angle spreads
  - Small-Scale Parameters (SSP): cluster delays, powers, angles of arrival/departure
  - Channel coefficient computation and updates over time

#### Link
- **Header**: [h/Link.h](h/Link.h)
- **Implementation**: [src/Link.cpp](src/Link.cpp)
- Manages BS-MS connections with:
  - Signal and interference computation
  - Antenna gain calculations (transmit and receive)
  - RSRP (Reference Signal Received Power) computation
  - Beam selection for MIMO systems

#### Sector
- **Header**: [h/sector.h](h/sector.h)
- **Implementation**: [src/sector.cpp](src/sector.cpp)
- Handles scheduling algorithms:
  - **Type 0**: Round-robin scheduling
  - **Type 1**: Proportional Fair (PF) with MU-MIMO
  - **Type 2**: mTRP (multi-TRP) with NC-JT (Non-Coherent Joint Transmission)
- Manages CQI/PMI feedback reading and MCS (Modulation and Coding Scheme) decision

### Key Modules

- **[src/generateLSP.cpp](src/generateLSP.cpp)**: Generates Large-Scale Parameters (shadow fading, delay spread, etc.)
- **[src/generateSSP.cpp](src/generateSSP.cpp)**: Generates Small-Scale Parameters (cluster delays, powers, angles)
- **[src/scheduling.cpp](src/scheduling.cpp)**: Implements multi-threaded scheduling with OpenMP
- **[src/receive_downlink.cpp](src/receive_downlink.cpp)**: Processes downlink reception and BLER computation
- **[src/measure.cpp](src/measure.cpp)**: Collects statistics and performance metrics
- **[src/allocateMemory.cpp](src/allocateMemory.cpp)**: Dynamic memory allocation for simulation arrays

### Configuration Files

Simulation parameters are read from `.cfg` files with settings for:
- **Scenario**: InH_eMBB, Dense_Urban_eMBB, Rural_eMBB (A/B/C variants)
- **Pathloss models**: InH_A/B, UMa_A/B, UMi_A/B, RMa_A/B, ETRI variants
- **Antenna configuration**: Number of elements (M, N, P), panel configuration (Mg, Ng), spacing (dH, dV)
- **Scheduling**: Type (0=RoundRobin, 1=MU-MIMO, 2=mTRP/NCJT)
- **OFDM parameters**: Numerology, bandwidth, carrier frequency, number of RBs
- **Simulation control**: Number of drops, run times, threads, calibration mode

Example: [Tbps_DU_B.cfg](Tbps_DU_B.cfg) configures a Dense Urban scenario with 275 RBs at 160 GHz.

### Global Variables

The simulation uses extensive global variables (defined in [h/common.h](h/common.h)):
- **System objects**: `bs`, `ms`, `links`, `channel`, `sector` arrays
- **Simulation state**: `drop_idx`, `t` (time index)
- **Configuration**: `num_BS`, `num_MS`, `num_SECTORS`, antenna parameters
- **Channel data**: Precalculated channel coefficients, PMI/CQI maps
- **Standards**: BLER tables, TBS (Transport Block Size) tables, SNR mappings

### Multi-threading

The simulator uses OpenMP for parallel processing:
- Channel coefficient computation is parallelized across links
- Scheduling across sectors runs in parallel (see [src/scheduling.cpp:45-52](src/scheduling.cpp))
- Number of threads controlled by `num_threads` parameter in config file

## Important Constants and Enumerations

Defined in [h/const.h](h/const.h):
- **Test environments**: `InH_eMBB_A/B/C`, `Dense_Urban_eMBB_A/B/C`, `Rural_eMBB_A/B/C`
- **Modulation types**: `QPSK`, `QAM16`, `QAM64`, `QAM256`
- **Propagation conditions**: `NLOS_propagation`, `LOS_propagation`, `OUT2IN_propagation`
- **Wraparound matrices**: `simple_wrap_mat_1tier`, `simple_wrap_mat_2tier` for hexagonal cell layouts

## Output

Results are written to directories specified in the config file (e.g., `Dense_Urban_eMBB_C_Micro_ETRI_1mTRP/`).

## Notes for Development

- **Memory management**: Allocation happens in `Allocate_memory()`, deallocation in `Delete_memory()`
- **Channel updates**: Occur every time slot via `Update_Channel_Coef()` in the main simulation loop
- **Calibration mode**: Set `Calibration_mode = 1` in config to run without full simulation (for debugging channel models)
- **BLER data**: Uses lookup tables from `BLER_5G.dat` for link-to-system mapping
- **Eigen library**: Used extensively for matrix operations in MIMO processing and precoding
