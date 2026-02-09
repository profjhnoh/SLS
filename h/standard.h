#include "const.h"

Real * SNR;
Real * SNR_dB;
Real * SNR_5G_dB;


char ** bler_file_name;
int NUM_OF_MCS = 38;


Real BLER_TABLE[38][7][2];
Real NEW_TABLE[38][7][2];

int TBS[51][27];

/// vienna BLER [ SINR_dB ][ CQI_idx ]
Real BLER_vienna[613][15];
Real NEW5GBLER[1427][15];


//// For EESM 
//// Reference -> "MCS Selection for Throughput Improvement in Downlink LTE Systems" Jiancun Fan, Qinye Yin, Geoffrey Ye Li, Bingguang Peng, and Xiaolong Zhu
Real Beta[15] = { 1.00, 1.40, 1.40, 1.48, 1.50, 1.62, 3.10, 4.32, 5.37, 7.71, 15.5, 19.6, 24.7, 27.6, 28.0};
//Real SINR_threshold_dB[15] = { -9.478, -6.658, -4.098, -1.798, 0.399, 2.424, 4.489, 6.367, 8.456, 10.266, 12.218, 14.122, 15.849, 17.786, 19.809 };

//// Vienna BLER curve 0.1
//Real SINR_threshold_dB[15] = { -6.9, -5.1, -3.15, -1.25, 0.8, 2.7, 4.7, 6.55, 8.6, 10.4, 12.3, 14.2, 15.9, 17.85, 19.85 };

//// 5G BLER curve 0.1
Real SINR_threshold_dB[15] = { -7.86, -6.2, -4.36, -1.9, 0.14, 2, 4.72, 6.22, 8.04, 11.04, 11.84, 13.5, 15.36, 17.44, 19.24 };

//Real re_scale_table[7] = { 1., 0.9011, 0.7001, pow(10., -1.), pow(10., -1.5), pow(10., -2.), 7.50 * pow(10., -3.) };
Real re_scale_table[8] = { 0,0,0,0,0,0,0,0 };

//                           0               1            2              3           4          5          6                    7         8             9           10             11            12          13            14              15          16           17         18      19     20
Real CR_curves[38] = { 0.083333333, 0.092447917, 0.104166667, 0.1171875, 0.130208333, 0.1484375, 0.166666667, 0.188802083, 0.213541667, 0.235677083, 0.266927083, 0.303385417, 0.341145833, 0.37890625, 0.427083333, 0.479166667, 0.541666667, 0.614583333, 0.6875, 0.75, 0.8
, 0.333333333, 0.390625, 0.427083333, 0.479166667, 0.541666667, 0.614583333, 0.666666666, 0.75, 0.8,
0.333333333, 0.458333333, 0.520833333, 0.614583333, 0.6875, 0.78125, 0.819444444, 0.916666667 };

Real RBIR_value[3][95];
Real Threshold[29];

Real log_linear_xais(Real a, Real f_a, Real b, Real f_b, Real c);

//function
Real Effective_CR_compute(int, int, int);
Real linear_cr(Real cr1, Real cr1_value, Real cr2, Real cr2_value, Real cr);
//Real linear_cr(Real, Real, Real, Real, Real);
Real linear_approx_xais(Real a, Real b, Real c, Real d, Real y);
void file_name_setting(void);
void table_setting_BLER(void);
void table_setting_TBS();
void table_setting_NEW_BLER();
void table_setting_Threshold();
void table_setting_RBIR();
void Load_BLER_vienna_file();
void Load_5G_BLER_vienna_file();