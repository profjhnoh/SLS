#include <cmath>
#include "precision.h"



#ifndef _CONST_
#define _CONST_
#define MIN(x,y) (((x)<=(y))?(x):(y))
#define MAX(x,y) (((x)>=(y))?(x):(y))
#define NO_UE -1
#define NO_CQI -1
#define dB2linear(x) (pow(10.0,x/10.0))
#define dBm2linear(x) (pow(10.0, (x- 30)/10.0))
#define linear2dB(x) (10*log10(x))
#define linear2dBm(x) (10*log10(x) + 30)

//const Real pi = 4.*atan(1.);
const Real pi = REAL(3.14159265358979323846);
const Real light_speed = REAL(300000000.);


enum test_environment_TYPE
{
	InH_eMBB_A, InH_eMBB_B, InH_eMBB_C,
	Dense_Urban_eMBB_A, Dense_Urban_eMBB_B, Dense_Urban_eMBB_C_Macro, Dense_Urban_eMBB_C_Micro,
	Rural_eMBB_A, Rural_eMBB_B, Rural_eMBB_C, Dense_Urban_eMBB_Single_Sector
};

enum pathloss_model_TYPE
{
	InH_A, InH_B,
	UMa_A, UMa_B,
	UMi_A, UMi_B,
	RMa_A, RMa_B,
	InH_ETRI, UMi_ETRI
};


const Real BS_feeder_loss = 2;
const Real thermal_noise = -174.; //[dBm/Hz]
//const Real Total_BS_Tx_power = 21; // [dBm] for 20 MHz
//const Real UE_BANDWIDTH = 15000.; // [Hz]
//const Real BS_noisefig = 5;
//const Real MS_noisefig = 7;
//const int num_drop = 50;

const int MAX_NUM_CLUSTERS = 24;
const int MAX_NUM_RAYS = 20;


//DFT
//const Real fft_size = 1024.;
const int num_freq_per_rbs = 12;


//standard

const int QPSK = 2;
const int QAM16 = 4;
const int QAM64 = 6;
const int QAM256 = 8;
const int num_ofdm_symbols_per_subband_per_1ms = 14;
const int inverse_TC_rate = 3;
const int num_crc = 24;

typedef struct
{
	/* data */
	int i_11;
	int i_12;
	int i_2;
} PMI_FEEDBACK;

typedef struct
{
	int ue_selected;
	int mcs_selected;
	int cqi_selected;
	Real capacity;
	int temp_cqi;
} SCHEDULE_DECISION;

typedef struct
{
	Real x;
	Real y;

}LOCATION;

typedef struct
{
	Real x;
	Real y;
	Real z;

}LOCATION3D;


typedef struct
{
	Real phi;
	Real theta;

} SPHERICAL_ANGLE;

typedef struct
{
	Real phi;
	Real theta;

} ANTENNA_FIELD;

/*
typedef struct
{
	ComplexReal phi;
	ComplexReal theta;

} ANTENNA_FIELD_COMPLEX;
*/

typedef struct
{
	int modulation_idx;
	Real effective_code_rate;
	Real snr;

} LL_CURVE_LOOK;

typedef struct
{
	Real code_rate;
	int mod_type;
	Real snr_scale_factor;
} INPUT_TO_LLS_MAPPING;


/*
const int wrap_mat[19][19] = {
{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18 },
{ 1, 8, 9, 2, 0, 6, 7, 126, 35, 34, 33, 10, 11, 3, 4, 5, 17, 18, 127},
{ 2, 9, 10, 11, 3, 0, 1, 8, 34, 33, 56, 55, 54, 12, 13, 4, 5, 6, 7 },
{ 3, 2, 11, 12, 13, 4, 0, 1, 9, 10, 55, 54, 65, 64, 75, 14, 15, 5, 6 },
{ 4, 0, 3, 13, 14, 15, 5, 6, 1, 2, 11, 12, 64, 75, 86, 85, 84, 16, 17},
{ 5, 6, 0, 4, 15, 16, 17, 18, 7, 1, 2, 3, 13, 14, 85, 84, 107, 106, 105},
{ 6, 7, 1, 0, 5, 17, 18, 127, 126, 8, 9, 2, 3, 4, 15, 16, 106, 105, 128},
{ 7, 126, 8, 1, 6, 18, 127, 117, 125, 35, 34, 9, 2, 0, 5, 17, 105, 128, 118},
{ 8, 35, 34, 9, 1, 7, 126, 125, 36, 24, 23, 33, 10, 2, 0, 6, 18, 127, 117 },
{ 9, 34, 33, 10, 2, 1, 8, 35, 24, 23, 32, 56, 55, 11, 3, 0, 6, 7, 126 },
{ 10, 33, 56, 55, 11, 2, 9, 34, 23, 32, 45, 44, 43, 54, 12, 3, 0, 1, 8 },
{ 11, 10, 55, 54, 12, 3, 2, 9, 33, 56, 44, 43, 53, 65, 64, 13, 4, 0, 1 },
{ 12, 11, 54, 65, 64, 13, 3, 2, 10, 55, 43, 53, 66, 58, 63, 75, 14, 4, 0 },
{ 13, 3, 12, 64, 75, 14, 4, 0, 2, 11, 54, 65, 58, 63, 74, 86, 85, 15, 5 },
{ 14, 4, 13, 75, 86, 85, 15, 5, 0, 3, 12, 64, 63, 74, 87, 78, 77, 84, 16 },
{ 15, 5, 4, 14, 85, 84, 16, 17, 6, 0, 3, 13, 75, 86, 78, 77, 83, 107, 106 },
{ 16, 17, 5, 15, 84, 107, 106, 105, 18, 6, 0, 4, 14, 85, 77, 83, 108, 98, 97 },
{ 17, 18, 6, 5, 16, 106, 105, 128, 127, 7, 1, 0, 4, 15, 84, 107, 98, 97, 104 },
{ 18, 127, 7, 6, 17, 105, 128, 118, 117, 126, 8, 1, 0, 5, 16, 106, 97, 104, 129 }


};*/

// ====================================================================
// 1-SITE (SINGLE CELL) CONFIGURATION - 3 sectors with 120° separation
// ====================================================================
// This configuration uses a single base station with 3 sectors arranged
// at 0°, 120°, 240° azimuths. To simulate inter-sector interference,
// the wraparound matrix places interfering sectors in opposite directions.
//
// Layout:
//     Sector 0: 0° (pointing right)
//     Sector 1: 120° (pointing upper-left)
//     Sector 2: 240° (pointing lower-left)
//
// Each sector sees interference primarily from the other two sectors
// of the same site.
// ====================================================================
const int simple_wrap_mat_1site[1][1] =
{
	{ 0 }  // Single site, no wraparound needed
};

const int simple_wrap_mat_1tier[7][7] =
{  { 0, 1, 2, 3, 4, 5, 6},
   { 1, 3, 5, 2, 0, 6, 4},
   { 2, 5, 4, 6, 3, 0, 1},
   { 3, 2, 6, 5, 1, 4, 0},
   { 4 ,0 ,3 ,1, 6, 2, 5},
   { 5, 6, 0 ,4 ,2 ,1, 3},
   { 6, 4, 1, 0, 5, 3, 2}
};

const int simple_wrap_mat_2tier[19][19] =

{
	{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18 },
	{ 1, 7, 8, 2, 0, 6, 18, 11, 15, 14, 9, 10, 3, 4, 5, 16, 17, 13, 12 },
	{ 2, 8, 9, 10, 3, 0, 1, 15, 14, 13, 17, 16, 11, 12, 4, 5, 6, 18, 7 },
	{ 3, 2, 10, 11, 12, 4, 0, 8, 9, 17, 16, 15, 7, 18, 13, 14, 5, 6, 1 },
	{ 4, 0, 3, 12, 13, 14, 5, 1, 2, 10, 11, 7, 18, 17, 9, 8, 15, 16, 6 },
	{ 5, 6, 0, 4, 14, 15, 16, 18, 1, 2, 3, 12, 13, 9, 8, 7, 11, 10, 17 },
	{ 6, 18, 1, 0, 5, 16, 17, 12, 7, 8, 2, 3, 4, 14, 15, 11, 10, 9, 13 },
	{ 7, 11, 15, 8, 1, 18, 12, 10, 16, 5, 14, 9, 2, 0, 6, 17, 13, 4, 3 },
	{ 8, 15, 14, 9, 2, 1, 7, 16, 5, 4, 13, 17, 10, 3, 0, 6, 18, 12, 11 },
	{ 9, 14, 13, 17, 10, 2, 8, 5, 4, 12, 18, 6, 16, 11, 3, 0, 1, 7, 15 },
	{ 10, 9, 17, 16, 11, 3, 2, 14, 13, 18, 6, 5, 15, 7, 12, 4, 0, 1, 8 },
	{ 11, 10, 16, 15, 7, 12, 3, 9, 17, 6, 5, 14, 8, 1, 18, 13, 4, 0, 2 },
	{ 12, 3, 11, 7, 18, 13, 4, 2, 10, 16, 15, 8, 1, 6, 17, 9, 14, 5, 0 },
	{ 13, 4, 12, 18, 17, 9, 14, 0, 3, 11, 7, 1, 6, 16, 10, 2, 8, 15, 5 },
	{ 14, 5, 4, 13, 9, 8, 15, 6, 0, 3, 12, 18, 17, 10, 2, 1, 7, 11, 16 },
	{ 15, 16, 5, 14, 8, 7, 11, 17, 6, 0, 4, 13, 9, 2, 1, 18, 12, 3, 10 },
	{ 16, 17, 6, 5, 15, 11, 10, 13, 18, 1, 0, 4, 14, 8, 7, 12, 3, 2, 9 },
	{ 17, 13, 18, 6, 16, 10, 9, 4, 12, 7, 1, 0, 5, 15, 11, 3, 2, 8, 14 },
	{ 18, 12, 7, 1, 6, 17, 13, 3, 11, 15, 8, 2, 0, 5, 16, 10, 9, 14, 4 }
};


enum CQI_MCS_TYPE
{
	QPSK0, QPSK1, QPSK2, QPSK3, QPSK4, QPSK5, QPSK6, QPSK7, QPSK8, QPSK9, QAM1610, QAM1611, QAM1612, QAM1613, QAM1614, QAM1615, QAM1616, QAM6417, QAM6418, QAM6419, QAM6420, QAM6421, QAM6422, QAM6423, QAM6424, QAM6425, QAM6426, QAM6427, QAM6428,
	NUM_CQI_MCS_TYPES
};

enum propagation_condition
{
	NLOS_propagation, LOS_propagation, OUT2IN_propagation
};

enum LOS_TLSP
{
	LOS_SF, LOS_K, LOS_DS, LOS_ASD, LOS_ASA, LOS_ZSD, LOS_ZSA
};

enum NLOS_TLSP
{
	NLOS_SF, NLOS_DS, NLOS_ASD, NLOS_ASA, NLOS_ZSD, NLOS_ZSA
};

enum OUT2IN_TLSP
{
	OUT2IN_SF, OUT2IN_DS, OUT2IN_ASD, OUT2IN_ASA, OUT2IN_ZSD, OUT2IN_ZSA
};

const Real ray_angle_mat[20] = { 0.0447, -0.0447, 0.1413, -0.1413, 0.2492, -0.2492, 0.3715, -0.3725, 0.5129, -0.5129, 0.6797, -0.6797, 0.8844, -0.8844, 1.1481, -1.1481, 1.5195, -1.5195, 2.1551, -2.1551 };

//#define WIN
#define UE_NOT_SPECIFIED -1


#endif