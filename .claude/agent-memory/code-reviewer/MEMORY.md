# Code Review Memory - 5G Downlink System Level Simulator

## Codebase Overview
- **Type**: High-performance 5G wireless channel simulator
- **Language**: C++11, OpenMP, Eigen3
- **LOC**: ~20,000+ lines across 20+ source files
- **Main file**: channel.cpp (5,445 lines)

## Architecture Patterns

### Memory Management
- Mixed manual allocation (`new`/`delete`) with Eigen smart types
- Lazy allocation pattern for CHIR arrays (sector-level, on-demand)
- CHIR_vec arrays commented out to save memory (~6.3 MB per channel)
- **ray_AOA/ray_AOD**: Always allocated in Allocate_memory() as [cluster][ray][2] for precomputed per-ray angles
- **ray_gain**: Conditionally allocated in Allocate_CHIR_memory() under USE_RAY_LEVEL_DOPPLER flag
- Deallocation follows reverse order with NULL guards in Delete_memory()

### Global Variable Design
- Extensive use of extern globals in common.h (~290 globals)
- System state shared across all modules: `bs`, `ms`, `links`, `channel`, `sector`
- Configuration parameters loaded from .cfg files

### Threading Model
- OpenMP parallelization at two levels: sector scheduling, MS processing
- Nested parallelism enabled (omp_set_nested)
- Thread-safe concerns: most operations read-only during parallel sections

### Channel Model Architecture
- 3GPP-compliant models: InH, UMa, UMi, RMa
- Two-phase: Large-Scale Parameters (LSP) → Small-Scale Parameters (SSP)
- Channel coefficients stored in 6D arrays: [sector][tx_port][rx_port][cluster][ray]

## Common Patterns

### Indexing Conventions
- BS-MS channels: `channel[bs_idx][ms_idx]`
- Antenna arrays: `d_tx[sector][M][N][P][Mg][Ng]` (3D physical arrays)
- Sector indexing: `adj_sector / 3` = BS index, `adj_sector % 3` = sector within BS

### Time Evolution
- Main loop: drop → time slot (t) → channel update → scheduling → reception
- Channel updates: cluster-average Doppler (fast) or ray-level precision (accurate)

## Code Quality Issues Found

### Critical
1. **channel.cpp line 3012**: Uninitialized comments with Korean encoding corruption
2. **Missing bounds checks**: No validation on array indices before access
3. **Memory leak risk**: CHIR arrays allocated but flags not always checked before deletion
4. **fmod() for angles** (channel.cpp:4336-4352): C++ fmod() has asymmetric behavior for negative inputs; must ensure `if (x < 0) x += modulus` after fmod() for correct mathematical modulo [0, modulus)

### Performance
1. **Complex nested loops**: 7-level nesting in channel coefficient computation (line ~4630)
2. **Doppler recalculation**: Random number generation inside hot loops
3. **No const correctness**: Member functions that don't modify state aren't marked const

### Maintainability
1. **Magic numbers**: Hardcoded constants (65, 70, 90) without named explanations
2. **Dead code**: Many commented-out blocks, unclear if intentional or forgotten
3. **Mixed languages**: Korean comments mixed with English

## Best Practices Observed
- Separation of initialization, generation, and computation phases
- CMake build system with precision/threading options
- Calibration mode for channel model validation
- FFT optimization to save memory
- Static helper functions for reusable algorithms (e.g., compute_circular_angle_spread)

## Refactoring Review (2026-02-12)
- See [refactoring-patterns.md](refactoring-patterns.md) for detailed v2 refactoring analysis
- **Key finding**: link_antgain field not set in v2 (backward compatibility break)
- **Bug fix**: Original line 333 typo (sector_a = sector_z) correctly fixed in v2
- **Compilation issue**: Undefined `signal` variable at line 4012 in g_comp_mode block

## Location Generation (2026-02-13)
- See [location-generation-algorithms.md](location-generation-algorithms.md) for complete BS/MS placement algorithms
- **BS layout**: InH fixed grid (2×6 or 3×8), hex grid for UMa/UMi/RMa (7-site/19-site)
- **MS placement**: InH uniform in rectangles, hex random in 60° sectors with rotation
- **Indoor probability**: Dense Urban 80% (below 20GHz), Rural 40-50%, InH 0%
- **UE heights**: Multi-floor for UMa indoor (4-8 floors), fixed 1.5m for outdoor
- **O2I loss**: Model A/B with 80/20 low/high loss split, frequency-dependent materials

## 3GPP Compliance Notes
- **Circular angular spread**: TR 25.996 Annex A requires delta optimization (minimize over angle shifts) to handle wrapping ambiguity
- **Cluster vs. ray metrics**: LSP angular spreads (ASA/ASD/ZSA/ZSD) should be computed from SSP cluster angles, NOT from raw LSP distributions
- **Subcluster timing**: Updated in 2026-02-11 review - circular angle spread now computed AFTER Set_SUBCLUSTER() at ray level (not cluster level)
- **Per-ray angles**: RX = AOA[i] + cluster_ASA * offset_angle[i][j], TX = AOD[i] + cluster_ASD * offset_angle_rand_coupling[i][j]
- **Subcluster structure**: 2 strongest clusters split into 3 subclusters each (10+6+4 rays), others remain 20 rays
- **Power splitting**: Each ray in subcluster gets power[i] / NUM_RAY_per_ClusterNUM[i]
