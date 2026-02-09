#include <iostream>
#include <fstream>

using namespace std;
#include "precision.h"
#include "rv_generate.h"
#include "BS.h"
#include "MS.h"
#include "Link.h"
#include "channel.h"
#include "sector.h"
#include <random>


unsigned long _seed = 100000;
Rand randnum((unsigned long long int)_seed);
default_random_engine e(_seed);

int num_drop;
int run_times;
int num_compute_coef;
int test_environment;
int pathloss_model;
int num_neighbor;
Real grid_interval;

int g_year;
int g_mon ;
int g_day ;
int g_hour;
int g_min ;


//ofdm param
int num_rb;
int fft_size;
int numerology;
Real subcarrier_spacing;
Real bandwidth;

// 검증용 인석 추가
Real sector_selected_ue;
Real scheduled_ue_mcs;
Real scheduled_ue_cqi;
Real scheduled_ue_widebandSINR;

Real ue_effective_sinr;
Real ue_num_traffic;
Real ue_info_bits;
Real ue_bler;

Real sector_metric;
Real sector_cqi_read;
Real sector_cqi_avr;

char *cfg_name;


int TYPE;
int num_BS;
int num_mTRP;

int num_MS;
int num_MS_persector;
int num_Indoor_TRxP;
int num_SECTORS;

int num_mTRP_SECTORS;

int num_LINK;
int BS_M, BS_N, BS_Mg, BS_Ng, BS_Mp, BS_Np, BS_P;  /// tx antenna element
Real BS_dH, BS_dV, BS_dgH, BS_dgV;
int MS_M, MS_N, MS_Mg, MS_Ng, MS_Mp, MS_Np, MS_P;  /// rx antenna element
Real MS_dH, MS_dV, MS_dgH, MS_dgV;
int NUM_RX;
int NUM_TX;
int NUM_TX_Port;
int NUM_RX_Port;
Real Mechanic_tilt;
int Scheduling_Type;
int scenario;
int Configuration_Type; /// For IMT-2020
int Channel_Model_Type; /// For IMT-2020
int CH_CAL;
int single_cell_mode = 0;  // 0: All cells, 1: Center cell only
int Calibration_mode;
int TDD_mode = 0;  // TDD reciprocity (1) vs FDD codebook (0)
int USE_POWER_ITERATION = 0;  // 0: SelfAdjointEigenSolver, 1: Power Iteration (faster but approximate)
int POWER_ITERATION_MAX_ITER = 5;  // Number of iterations for Power Iteration method

ComplexReal w[4][8];                              // [num_vertical steering angle][num element in a port] (maximum)
ComplexReal v[4][8][8];                           // [num_vertical steering angle][num_horizontal sterring angle][num element in a port]  (maximum)
ComplexReal virtualization_weight_wv[4][8][8][8]; // [vertical steering angle ][horizontal steering angle][vertical_element per port][horizontal_element per port]

ComplexReal ue_w[4][8];                              // [num_vertical steering angle][num element in a port] (maximum)
ComplexReal ue_v[4][8][8];                           // [num_vertical steering angle][num_horizontal sterring angle][num element in a port]  (maximum)
ComplexReal ue_virtualization_weight_wv[4][8][8][8]; // [vertical steering angle ][horizontal steering angle][vertical_element per port][horizontal_element per port]

int tilt_azimuth_angle_LCS_size;
int tilt_zenith_angle_LCS_size;

int ue_tilt_azimuth_angle_LCS_size;
int ue_tilt_zenith_angle_LCS_size;

Real bs_tilt_azimuth_angle_LCS[8];  //// pi  max8
Real bs_tilt_zenith_angle_LCS[4];  /// theta  max4

Real ue_tilt_azimuth_angle_LCS[4];  //// pi  max4
Real ue_tilt_zenith_angle_LCS[2];  /// theta  max2

ComplexReal * TX_LOS_gain_theta;
ComplexReal * TX_LOS_gain_pi;
ComplexReal *** TX_NLOS_gain_theta;
ComplexReal *** TX_NLOS_gain_pi;

ComplexReal * RX_LOS_gain_theta;
ComplexReal * RX_LOS_gain_pi;
ComplexReal *** RX_NLOS_gain_theta;
ComplexReal *** RX_NLOS_gain_pi;

char * file_name;
char folder_name[50];


Real ** Capacity_avr_per_drop;

Real user_speed;
Real W;
Real h;
Real carrier_freq;
Real Total_BS_Tx_power;
Real Wavelength;
Real bs_maxpower;
Real micro_bs_power;

Real ms_maxpower;
Real max_antgain;

Real noise;
Real BS_noisefig;
Real MS_noisefig;
Real ue_antenna_element_gain;
Real bs_height;

Real cfg_BS_Tx_Power;
Real cfg_UT_Noise_Figure;

Real macro_bs_height;
Real micro_bs_height;

Real ms_height;
Real ms_height_out;
Real ms_height_in;

int num_floor;
int num_propagation_condition;

Real inter_site_distance;
Real min_distance;
Real ANGLE_tilt;
Real MAX_v_angle;

void Set_Parameter();
void Initialdrop();

// MU MIMO
//int               *** pppPMI_map;
//VectorXcReal         *** pppPMI_vector_map;
//Real            *** pppCQI_Map;

// 220815 jhnoh
PMI_FEEDBACK      **** ppppPMI_map;
VectorXcReal      **** ppppPMI_vector_map;
MatrixXcReal      **** ppppCSI_matrix_map;  // TDD: Full channel matrices
Real              **** ppppCQI_Map;
Real              **** ppppCQI_comp_Map;

// simulation configure variables
int cqi_history_length;
int       mx_ue_mumimo;
int       NUM_UE_Layer;

Real N_pf;
int SCHEDULE_DELAY;

// search serving cell - only near 7 cell, not 19 cell
int simple_num_BS;


Real epsilon = 0.5;   // Semi-orthogonality threshold for SUS (typical: 0.3~0.7)


int    mcs_decision = 1;
Matrix4cReal WW[16];
Vector4cReal code_word[64];

// jhnoh 221008
MatrixXcReal *** codebook_W;
MatrixXcReal cb_W_csirs_2_layer_1[4];
MatrixXcReal cb_W_csirs_2_layer_2[2];

// simulation run variables
int t;
int drop_idx;





int precalculate_time = 2000;
int anglev_mapping_value = 3;


ComplexReal **** CHIR_precalc[2][2];
MatrixXcReal **Hm[2][2];
int ** PMI_precalc[20000][2][2];
Real **sinr_estimate_precalc[20000][2][2];
Real **sum_of_ch_coef_precalc[20000][2][2];



//standard
Real    overhead = 0.33;
bool      calibration = false;

// CQI & Schedule Map
SCHEDULE_DECISION *** pppSchedule_Map;


//int SCHEDULE_TYPE = 1;

//receive downlink
int lls_mapping = 0;

//logging point
char path[50];

BS* bs;

BS* mTRP;

MS* ms;
LINK* links;
CHANNEL** channel;
Sector* sector;
Sector* mTRP_sector;


int num_of_threads;

int * comp_mode;
int g_comp_mode;
int comp_ue_pct;
int g_static_gain_ratio_comp;

// SINR aggregation mode: 0 = average (default), 1 = minimum (conservative)
int g_use_min_sinr_for_mcs;

// OLLA (Outer Loop Link Adaptation) parameters
int  g_olla_enable;
Real g_olla_step_down;
Real g_olla_step_up;
Real g_olla_min_offset;
Real g_olla_max_offset;
Real g_olla_init_offset;
Real g_olla_target_bler;
Real g_olla_bler_margin;
int  g_olla_window_size;

// MU-MIMO Scheduling Algorithm Selection
int  g_mumimo_scheduling_algorithm;
Real g_chordal_alpha;

// Singular Value CDF Collection
int g_collect_singular_values;  // 0 = disable, 1 = enable

int g_mTRP_mode;

ofstream  Calibration_Debug_info;

int ** simple_wrap_mat;
vector <int> sector_indices_schedule;


