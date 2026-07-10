// ENABLE_MULTITHREADING is now controlled via CMake option (-DENABLE_MULTITHREADING=ON/OFF)

#include <complex>
#include <fstream>
#include <iostream>
#include <string>
#include <new>

using namespace std;

#include "precision.h"
#include "MS.h"
#include "BS.h"
#include "rv_generate.h"
#include "Link.h"
#include "channel.h"
#include "sector.h"
#ifdef ENABLE_PROGRESSBAR
#include "indicators.hpp"
#endif
#include <random>

extern unsigned long _seed;
extern Rand randnum;
extern default_random_engine e;

extern int num_drop;
extern int run_times;
extern int test_environment;
extern int pathloss_model;
extern int num_neighbor;
extern Real grid_interval;

extern int g_year;
extern int g_mon ;
extern int g_day ;
extern int g_hour;
extern int g_min ;

//ofdm param
extern int num_rb;
extern int fft_size;
extern int numerology;
extern Real subcarrier_spacing;
extern Real bandwidth;

// 검증용 인석 추가
extern Real sector_selected_ue;
extern Real scheduled_ue_mcs;
extern Real scheduled_ue_cqi;
extern Real scheduled_ue_widebandSINR;

extern Real ue_effective_sinr;
extern Real ue_num_traffic;
extern Real ue_info_bits;
extern Real ue_bler;

extern Real sector_metric;
extern Real sector_cqi_read;
extern Real sector_cqi_avr;

extern char *cfg_name;

extern Real user_speed;
extern int CH_CAL;
extern int Calibration_mode;
extern int TDD_mode;  // TDD reciprocity (1) vs FDD codebook (0)
extern int USE_POWER_ITERATION;  // 0: SelfAdjointEigenSolver, 1: Power Iteration
extern int POWER_ITERATION_MAX_ITER;  // Number of iterations for Power Iteration
extern int USE_RAY_LEVEL_DOPPLER;  // 0: Cluster-average Doppler (fast), 1: Ray-level precision (accurate)
extern int USE_PRECODING_BASED_SINR;  // 0: Use feedback CQI (default), 1: Calculate SINR from precoding matrix
extern Real INTERCELL_INTERFERENCE_MARGIN_DB;  // Inter-cell interference margin in dB (for future use)
extern int Scheduling_Type;
extern int scenario;
extern int Configuration_Type; /// For IMT-2020
extern int Channel_Model_Type; /// For IMT-2020
extern int TYPE;
extern int num_BS;

// Single Cell Analysis Mode
// When enabled (1), only the center BS (BS 0) is scheduled and analyzed.
// Surrounding BSs (BS 1-6 in 7-site) provide interference only.
extern int single_cell_mode;

extern int num_mTRP;

extern int num_MS;
extern int num_MS_persector;
extern int num_Indoor_TRxP;
extern int num_SECTORS;

extern int num_mTRP_SECTORS;

extern int num_LINK;
extern int NUM_RX;
extern int NUM_TX;
extern int NUM_RX_Port;
extern int NUM_TX_Port;
extern Real Mechanic_tilt;
extern int BS_M, BS_N, BS_Mg, BS_Ng, BS_Mp, BS_Np, BS_P;  /// tx antenna element
extern Real BS_dH, BS_dV, BS_dgH, BS_dgV;
extern int MS_M, MS_N, MS_Mg, MS_Ng, MS_Mp, MS_Np, MS_P;  /// rx antenna element
extern Real MS_dH, MS_dV, MS_dgH, MS_dgV;

extern ComplexReal w[8][8];                              // [num_vertical steering angle][num element in a port] (maximum)
extern ComplexReal v[8][8][8];                           // [num_vertical steering angle][num_horizontal sterring angle][num element in a port]  (maximum)
extern ComplexReal virtualization_weight_wv[8][8][8][8]; // [vertical steering angle ][horizontal steering angle][vertical_element per port][horizontal_element per port]

extern ComplexReal ue_w[4][8];                              // [num_vertical steering angle][num element in a port] (maximum)
extern ComplexReal ue_v[4][8][8];                           // [num_vertical steering angle][num_horizontal sterring angle][num element in a port]  (maximum)
extern ComplexReal ue_virtualization_weight_wv[4][8][8][8]; // [vertical steering angle ][horizontal steering angle][vertical_element per port][horizontal_element per port]

extern int tilt_azimuth_angle_LCS_size;
extern int tilt_zenith_angle_LCS_size;
extern int ue_tilt_azimuth_angle_LCS_size;
extern int ue_tilt_zenith_angle_LCS_size;

extern Real bs_tilt_azimuth_angle_LCS[8];  //// pi  max8
extern Real bs_tilt_zenith_angle_LCS[8];  /// theta  max8
extern Real ue_tilt_azimuth_angle_LCS[4];  //// pi
extern Real ue_tilt_zenith_angle_LCS[2];  /// theta

extern ComplexReal * TX_LOS_gain_theta;
extern ComplexReal * TX_LOS_gain_pi;
extern ComplexReal *** TX_NLOS_gain_theta;
extern ComplexReal *** TX_NLOS_gain_pi;

extern ComplexReal * RX_LOS_gain_theta;
extern ComplexReal * RX_LOS_gain_pi;
extern ComplexReal *** RX_NLOS_gain_theta;
extern ComplexReal *** RX_NLOS_gain_pi;

extern char *    file_name;
extern char      folder_name[50];
extern Real ** Capacity_avr_per_drop;

extern Real W;
extern Real h;
extern Real carrier_freq;
extern Real Wavelength;
extern Real bs_maxpower;
extern Real micro_bs_power;

extern Real ms_maxpower;
extern Real max_antgain;
extern Real noise;
extern Real BS_noisefig;
extern Real MS_noisefig;
extern Real ue_antenna_element_gain;
extern Real Total_BS_Tx_power;

extern Real cfg_BS_Tx_Power;
extern Real cfg_UT_Noise_Figure;
extern Real cfg_UE_antenna_element_gain;
extern Real cfg_inter_site_distance;

extern int num_floor;
extern int num_propagation_condition;

extern Real bs_height;

extern Real macro_bs_height;
extern Real micro_bs_height;

extern Real ms_height;
extern Real ms_height_out;
extern Real ms_height_in;
extern Real inter_site_distance;
extern Real min_distance;
extern Real ANGLE_tilt;
extern Real MAX_v_angle;

//mu_mimo..
extern Real N_pf;
extern int cqi_history_length;
extern Real            **** ppppCQI_Map;
extern Real            **** ppppCQI_comp_Map;
extern PMI_FEEDBACK    **** ppppPMI_map;
extern MatrixXcReal               **** ppppPMI_vector_map;   // N_tx × R (R=g_type2_rank), col 0 is dominant layer
extern MatrixXcReal               **** ppppCSI_matrix_map;  // TDD: Full channel matrices for reciprocity-based precoding
//extern int SCHEDULE_TYPE;
extern int    mcs_decision;

extern int       mx_ue_mumimo;
extern int       NUM_UE_Layer;

extern Real epsilon;
extern 	Matrix4cReal WW[16];
extern Vector4cReal code_word[64];

// jhnoh 221008
extern MatrixXcReal *** codebook_W;
extern MatrixXcReal cb_W_csirs_2_layer_1[4];
extern MatrixXcReal cb_W_csirs_2_layer_2[2];

extern MatrixXcReal **Hm[2][2];
extern int ** PMI_precalc[20000][2][2];
extern ComplexReal **** CHIR_precalc[2][2];
extern Real **sinr_estimate_precalc[20000][2][2];
extern Real **sum_of_ch_coef_precalc[20000][2][2];

// CQI & Schedule Map
extern SCHEDULE_DECISION *** pppSchedule_Map;
//extern Real            *** pppCQI_Map;
//extern int               *** pppPMI_map;



// simulation run variables
extern int t;
extern int drop_idx;

// simulation configure variables
extern int SCHEDULE_DELAY;
extern int num_compute_coef;
extern int precalculate_time;
extern int anglev_mapping_value;
extern int simple_num_BS;

//standard
extern Real overhead;
extern Real BLER_vienna[613][15];
extern Real NEW5GBLER[1427][15];
extern int    TBS[51][27];
extern Real NEW_TABLE[38][7][2];
extern Real BLER_TABLE[38][7][2];
extern bool   calibration;
extern Real * SNR;
extern Real * SNR_dB;
extern Real * SNR_5G_dB;

//extern int num_path;
//extern Real * delay;
//extern ComplexReal *** CHIR;

//receive_downlink
extern int lls_mapping;
extern Real CR_curves[38];
extern Real Beta[15];
extern Real SINR_threshold_dB[15];
extern Real RBIR_value[3][95];
extern Real Threshold[29];

//logging point
extern char path[50];
#ifdef ENABLE_PROGRESSBAR
extern indicators::ProgressBar * bar;
extern int old_time_percent;
#endif

extern BS* bs;
extern MS* ms;
extern LINK* links;
extern CHANNEL** channel;
extern Sector* sector;
extern Sector* mTRP_sector;



extern int num_of_threads;

// SINR aggregation mode: 0 = average (default), 1 = minimum (conservative)
extern int g_use_min_sinr_for_mcs;

// OLLA (Outer Loop Link Adaptation) parameters
extern int  g_olla_enable;        // 0 = disable, 1 = enable
extern Real g_olla_step_down;     // SINR offset decrease when BLER too high (dB)
extern Real g_olla_step_up;       // SINR offset increase when BLER too low (dB)
extern Real g_olla_min_offset;    // Minimum OLLA offset (dB)
extern Real g_olla_max_offset;    // Maximum OLLA offset (dB)
extern Real g_olla_init_offset;   // Initial OLLA offset (dB)
extern Real g_olla_target_bler;   // Target BLER (e.g., 0.1 for 10%)
extern Real g_olla_bler_margin;   // Margin for BLER comparison (no adjustment if |measured - target| < margin)
extern int  g_olla_window_size;   // Moving window size for BLER measurement

// MU-MIMO Scheduling Algorithm Selection
extern int  g_mumimo_scheduling_algorithm;  // 0 = SUS (default), 1 = Chordal Distance
extern Real g_chordal_alpha;                // Weight factor for Chordal Distance (Frobenius norm vs orthogonality)

// Singular Value CDF Collection
extern int g_collect_singular_values;  // 0 = disable, 1 = enable

// Cluster parameter CSV dump
extern int g_dump_cluster_params;  // 0 = disable, 1 = enable

// Channel parameter version: 0=V19 (default), 1=old TR 38.901 (pre-V19)
extern int channel_param_legacy;

// Handheld UT antenna model (TR 38.901 Table 7.8-2A)
extern int handheld_mode;            // 0=off, 1=handheld
extern int handheld_num_ports;       // active ports (e.g., 4 for Config B)
extern int handheld_port_indices[8]; // 1-based antenna position indices
extern Real handheld_beta_deg;       // UT β angle (degrees)

// Type II Codebook (TS 38.214 §5.2.2.2.3) parameters
// Codebook_Type: 1 = Type I (default), 2 = Type II (Rel-15, rank 1 only), 3 = Enhanced Type II (Rel-16, rank 1-4)
extern int  g_codebook_type;
extern int  g_type2_L;                  // Number of selected beams per rotation (2/3/4)
extern int  g_type2_phase_alphabet;     // Phase alphabet size (4 = QPSK, 8 = 8-PSK, 16 = 16-PSK for eType II)
extern int  g_type2_subband_amplitude;  // 0 = WB amplitude only, 1 = SB amplitude enabled (Rel-15)

// Rel-16 Enhanced Type II parameters (TS 38.214 §5.2.2.2.5)
extern int  g_type2_rank;               // MAX rank cap R_max (1..4). Per-UE rank may be lower.
extern int  g_etype2_param_comb;        // ParameterCombination index (1..8): determines (L, p_v, β)
extern Real g_etype2_pv;                // FD basis ratio p_v (e.g. 1/4, 1/2)
extern Real g_etype2_beta;              // Sparsity ratio β (e.g. 1/4, 1/2, 3/4)

// Rank-adaptive RI selection (per-UE rank decision at MS side)
extern int  g_rank_adaptive;            // 0 = fixed rank = g_type2_rank, 1 = per-UE RI from SVD capacity
extern int  g_su_fallback;              // 0 = always MU-MIMO, 1 = compare SU metric vs MU and pick max

// Receiver-side SIC option for own UE's multiple layers
//   0 = per-stream MMSE (own UE's other layers counted as interference)
//   1 = ideal SIC (own UE's other layers cancelled; current historical behavior)
extern int  g_use_sic;

// Per-layer MCS + per-layer HARQ (each of a UE's rank-R layers gets its own MCS,
// sub-TBS, EESM BLER and independent decode/retx).
//   0 = legacy single-MCS / single-TB HARQ (all layers share one MCS)
//   1 = per-layer MCS + per-layer HARQ
extern int  g_per_layer_mcs;

// HARQ soft combining (per-layer). Requires g_per_layer_mcs=1.
//   0 = no combining (HARQ Type I): each retx is an independent trial at the same BLER
//   1 = IR/Chase: accumulate per-layer effective SINR across (re)transmissions, look up
//       BLER at the accumulated SINR → retransmissions progressively more likely to decode
extern int  g_harq_ir;

// MATLAB 5G Toolbox L2SM BLER tables (code-rate + TBS aware; see h/nr_l2sm.h).
// Both tables stay loaded so runs are A/B comparable; flags default OFF.
extern int  g_matlab_bler;            // BLER coin from MATLAB curves instead of NEW5GBLER
extern int  g_matlab_cqi_thresholds;  // regenerate SINR_threshold_dB from MATLAB curves
extern int  g_matlab_bler_selftest;   // startup oracle validation
extern int  g_matlab_rbir;            // RBIR effective SINR (curves' native axis) when matlab_bler=1
extern int  g_matlab_tput_mcs;        // throughput-maximizing MCS over all 28 MCS (requires matlab_bler=1)
extern int  g_tput_mcs_ref_rbs;       // reference per-UE RBs for the MCS grid; 0 = auto (num_rb/3)
extern int  g_matlab_esinr_fb;        // realized-ESINR feedback correction (measurement-driven link adaptation)

// Type II codebook cached DFT beams (cleared in Type2_Codebook_Gen)
//   type2_beam_v[l][m] = u_m ⊗ u_l  (length N1*N2)
extern VectorXcReal ** type2_beam_v;
extern int type2_N1, type2_N2, type2_O1, type2_O2;

// BS-side Spatial Non-Stationarity (SNS) parameters
extern int  g_sns_bs_enabled;
extern Real g_sns_mu_P_vis;
extern Real g_sns_sigma_P_vis;
extern Real g_sns_vr_A;
extern Real g_sns_vr_B;
extern Real g_sns_vr_R;
extern Real g_sns_vr_delta;
extern Real g_sns_rolloff_C;

// Row/Column analog beam allocation (per-drop population vote)
extern int  row_beam_enable;
extern int  row_beam_az_mode;
extern int  row_beam_force_uniform;
extern int  row_beam_max_cand;
extern Real row_beam_x_db;
extern Real row_beam_zenith_min_deg;
extern Real row_beam_zenith_max_deg;
extern int  row_beam_num_zenith;
extern int  row_beam_boresight_a;

extern ofstream  Calibration_Debug_info;

extern int ** simple_wrap_mat;
extern vector <int> sector_indices_schedule;


