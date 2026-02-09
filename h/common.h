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

extern ComplexReal w[4][8];                              // [num_vertical steering angle][num element in a port] (maximum)
extern ComplexReal v[4][8][8];                           // [num_vertical steering angle][num_horizontal sterring angle][num element in a port]  (maximum)
extern ComplexReal virtualization_weight_wv[4][8][8][8]; // [vertical steering angle ][horizontal steering angle][vertical_element per port][horizontal_element per port]

extern ComplexReal ue_w[4][8];                              // [num_vertical steering angle][num element in a port] (maximum)
extern ComplexReal ue_v[4][8][8];                           // [num_vertical steering angle][num_horizontal sterring angle][num element in a port]  (maximum)
extern ComplexReal ue_virtualization_weight_wv[4][8][8][8]; // [vertical steering angle ][horizontal steering angle][vertical_element per port][horizontal_element per port]

extern int tilt_azimuth_angle_LCS_size;
extern int tilt_zenith_angle_LCS_size;
extern int ue_tilt_azimuth_angle_LCS_size;
extern int ue_tilt_zenith_angle_LCS_size;

extern Real bs_tilt_azimuth_angle_LCS[6];  //// pi
extern Real bs_tilt_zenith_angle_LCS[2];  /// theta
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
extern VectorXcReal               **** ppppPMI_vector_map;
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
extern int * comp_mode;
extern int g_comp_mode;
extern int comp_ue_pct;
extern int g_static_gain_ratio_comp;

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

extern int g_mTRP_mode;

extern ofstream  Calibration_Debug_info;

extern int ** simple_wrap_mat;
extern vector <int> sector_indices_schedule;


