#ifndef _NR_L2SM_H_
#define _NR_L2SM_H_

#include "precision.h"
#include <vector>

// nr_l2sm — MATLAB 5G Toolbox L2SM BLER-table integration.
//
// Replaces the legacy (SINR x CQI) BLER table with the MATLAB AWGN code-block
// BLER curves keyed by (BGN, R*1024, Qm, Zc). TBS enters through TS 38.212
// §5.2.2 segmentation (BGN/Zc/C) and the transport-block combination
//   TB_BLER = 1 - (1 - cbBLER)^C.
// Table file "matlab_l2sm_bler.dat" is produced by export_matlab_bler.m and is
// MathWorks-proprietary data: research use only, never in the distribution.
//
// All lookups are read-only after Load_Matlab_L2SM_BLER() (called from
// StandardInitialization, before the OpenMP receive region ever runs).
// Interpolation math is double precision regardless of the Real typedef.

struct NrSegmentInfo
{
	int bgn;   // LDPC base graph number (1 or 2)
	int C;     // number of code blocks
	int Zc;    // lifting size
	int K;     // code block size (22*Zc for BG1, 10*Zc for BG2)
};

// TS 38.212 §5.2.2 segmentation (nrDLSCHInfo equivalent).
// tbs_bits = transport block size A (bits), code_rate = R as a fraction.
NrSegmentInfo nr_dlsch_segment_info(int tbs_bits, double code_rate);

// Single source of truth for MCS index -> (Qm bits, target code rate x 1024).
// TS 38.214 Table 5.1.3.1-2 subset used by this simulator (28 entries).
// qm_bits equals the QPSK/QAM16/QAM64/QAM256 constants (2/4/6/8).
void nr_mcs_to_qm_r(int mcs_idx, int& qm_bits, Real& r_x1024);

// NR TBS quantization (TS 38.214 §5.1.3.2) shared by the per-layer TBS path
// and the CQI-threshold regeneration. Mirrors the legacy math exactly.
int nr_tbs_info_bits_from_mcs(int mcs_idx, int N_RE, int layer_factor);

// Loads matlab_l2sm_bler.dat. Returns false if missing/corrupt.
bool Load_Matlab_L2SM_BLER(const char* path);
bool Matlab_L2SM_loaded(void);

// Transport-block BLER for (effective SINR, MCS, TBS). Main receive-path entry.
Real bler_lookup_matlab(Real esinr_linear, int mcs_idx, int tbs_bits);

// Core lookup on explicit (Qm, R*1024); used by threshold regeneration and the
// selftest. Returns TB BLER; optionally reports the code-block BLER and the
// segmentation actually used.
double bler_lookup_matlab_core(double esinr_dB, int qm_bits, double r_x1024, int tbs_bits,
                               double* out_cb_bler = 0, NrSegmentInfo* out_seg = 0,
                               int* out_curve_r1024 = 0);

// Overwrites the global SINR_threshold_dB[15] with SNR@ (TB BLER = 0.1) derived
// from the MATLAB tables for each CQI's mapped MCS at a reference full-band
// single-layer allocation. Prints old vs new.
void Regenerate_SINR_thresholds_from_matlab(void);

// Validates the C++ lookup against the MATLAB-generated oracle CSV
// (bler_compare/matlab_reference_tuples.csv). Returns true if all rows pass.
bool Matlab_BLER_selftest(const char* tuples_csv);

// ---- RBIR effective-SINR mapping (matches the AWGN curves' x-axis) --------
// The MATLAB curves are tabulated against an RBIR-effective SINR; feeding them
// an EESM(legacy-Beta) effective SINR biases the steep LDPC waterfalls (the
// documented R1 risk — confirmed by A/B). These implement MATLAB's
// wireless.internal.L2SM.calculateEffectiveSINR (alpha=beta=1, dB domain):
// clamp each SINR to the table range, forward-map to RBIR, average, invert.

bool Load_RBIR_ESM(const char* path);            // rbir_esm.dat from export_rbir_esm.m
bool RBIR_ESM_loaded(void);

// Per-RB linear SINRs -> RBIR-effective SINR (linear). qm_bits in {2,4,6,8}.
Real rbir_effective_sinr_linear(const std::vector<Real>& sinr_linear, int qm_bits);

// dB-domain core (used by the linear wrapper and the selftest).
double rbir_effective_sinr_dB(const double* sinr_dB, int n, int qm_bits);

// Validates against bler_compare/rbir_reference_tuples.csv.
bool Matlab_RBIR_selftest(const char* tuples_csv);

// ---- Throughput-maximizing MCS selection --------------------------------
// argmax over all 28 MCS of  TBS(m) x (1 - TB_BLER(esinr, m, TBS(m)))  using
// the MATLAB curves. Unlike the CQI->MCS map (which reaches only 15 of the 28
// MCS and inherits the CQI-grid coarseness), this uses the full MCS table and
// is TBS/segmentation aware. TBS(m) is evaluated at a REPRESENTATIVE per-UE
// allocation (n_re_ref), single layer (matches the per-layer receiver TBS).
// A dB-grid decision table is precomputed at init so the per-(UE,RB) cost in
// the scheduler is a single array read (grid verified against direct argmax
// at build time).

void Build_TputMCS_Grid(int n_re_ref);        // call once after table loads
bool TputMCS_Grid_ready(void);
int  nr_mcs_maxtput_from_dB(Real esinr_dB);   // grid lookup (clamped)

// Shared CQI(1..15) -> MCS map (mirrors determine_MCS's hardcoded chain).
// Index by CQI; element [0] is unused.
extern const int NR_CQI2MCS[16];

#endif
