#include "standard.h"
#include "common.h"
#include "nr_l2sm.h"
#include <iostream>
#include <fstream>


using namespace std;

/*===================================================================
FUNCTION: linear_approx_xais

2009-10-27
===================================================================*/
Real linear_approx_xais(Real a, Real b, Real c, Real d, Real y)
{
	return ((a*(b - d) + (b - y)*(c - a)) / (b - d));
}

/*===================================================================
FUNCTION: linear_approx_yais

2009-10-27
===================================================================*/
Real linear_approx_yais(Real a, Real b, Real c, Real d, Real x)
{
	return ((b*(c - a) - (b - d)*(x - a)) / (c - a));
}

/*===================================================================
FUNCTION: linear_cr

2009-10-27
===================================================================*/
Real linear_cr(Real cr1, Real cr1_value, Real cr2, Real cr2_value, Real cr)
{
	return ((cr2_value*(cr - cr1) + cr1_value*(cr2 - cr)) / (cr2 - cr1));
}

/*===================================================================
FUNCTION: log_linear_xais() , CR_computer()

AUTHOR: KR

DESCRIPTION: function for NEW_BLER table & Threshold

NOTES: Updated 7/10

===================================================================*/
Real log_linear_xais(Real a, Real f_a, Real b, Real f_b, Real f_c) // return (c)
{
	Real k, t;
	k = log10(f_a) - log10(f_b);
	t = log10(f_c) - log10(f_b);

	return((k*b + t*a - t*b) / k);
}

/*===================================================================
FUNCTION: StandardInitialization()

AUTHOR: SJ, KR, DW

DESCRIPTION: Input from the Standard & Link-Level Curve

NOTES: Updated 7/10 -> NEW_TABLE, Threshold

===================================================================*/
void StandardInitialization()
{
	table_setting_BLER();
	
	table_setting_TBS(); //transport block size  // TS 36.213(Rel 8) - 7.1.7.2.1
	//table_setting_NEW_BLER();
	table_setting_Threshold();//for MCS decision
	table_setting_RBIR();

	//// vienna BLER curve
	//Load_BLER_vienna_file();
	Load_5G_BLER_vienna_file();

	// MATLAB L2SM BLER tables (code-rate + TBS aware). The legacy table above
	// stays loaded so matlab_bler=0 runs are untouched (side-by-side A/B).
	if (g_matlab_bler || g_matlab_cqi_thresholds || g_matlab_bler_selftest)
	{
		if (!Load_Matlab_L2SM_BLER("matlab_l2sm_bler.dat"))
		{
			// Abort: with the table missing every bler_lookup_matlab returns 1.0,
			// so the run would "complete" with all-failing transport blocks. The
			// old getchar() also hung headless/batch runs on a dead prompt.
			cout << "=== FATAL: matlab_l2sm_bler.dat missing/corrupt: run export_matlab_bler.m first ===" << endl;
			exit(1);
		}
	}
	if ((g_matlab_bler && g_matlab_rbir) || g_matlab_bler_selftest)
	{
		if (!Load_RBIR_ESM("rbir_esm.dat"))
		{
			cout << "=== FATAL: rbir_esm.dat missing/corrupt: run export_rbir_esm.m first ===" << endl;
			exit(1);
		}
	}
	if (g_matlab_bler_selftest)
	{
		Matlab_BLER_selftest("bler_compare/matlab_reference_tuples.csv");
		Matlab_RBIR_selftest("bler_compare/rbir_reference_tuples.csv");
	}
	if (g_matlab_cqi_thresholds)
	{
		Regenerate_SINR_thresholds_from_matlab();
	}
	if (g_matlab_bler && g_matlab_tput_mcs)
	{
		// Reference per-UE allocation for the MCS decision grid: with ~mx/RI
		// UEs sharing the band each slot, a UE typically gets ~num_rb/3 RBs.
		int ref_rbs = (g_tput_mcs_ref_rbs > 0) ? g_tput_mcs_ref_rbs : (num_rb / 3 > 0 ? num_rb / 3 : 1);
		int n_re_ref = MIN(156, num_freq_per_rbs * num_ofdm_symbols_per_subband_per_1ms) * ref_rbs;
		Build_TputMCS_Grid(n_re_ref);
	}

	SNR = new Real[95];
	for (int i = 0; i < 95; i++)
	{
		SNR[i] = 0.5 * i - 20.;
	}

	/*
	// vienna BLER x-axis(SNR_dB)// -9.5 ~ 21.05  (+0.05)
	SNR_dB = new Real[612];
	for (int i = 0; i < 612; i++)
	{
		SNR_dB[i] = 0.05 * i - 9.5;
	}
	*/

	// 5G BLER x-axis(SNR_dB)// -9 ~ 19.5  (+0.02)
	SNR_5G_dB = new Real[1426];
	for (int i = 0; i < 1426; i++)
	{
		SNR_5G_dB[i] = 0.02 * i - 9;
	}
	
}
/*===================================================================
FUNCTION: table_setting_NEW_BLER()

AUTHOR: SJ, KR, DW

DESCRIPTION:

NOTES:

===================================================================*/
void table_setting_NEW_BLER()
{
	// NEW_BLER table generator

	//target computer	
	int  scale_target[38][7] = { 0 }; 
	
	for (int mcs_index = 0; mcs_index < 38; mcs_index++)
	{
		for (int scale_idx = 0; scale_idx < 6; scale_idx++)
		{
			for (int target_idx = 0; target_idx < 6; target_idx++)
			{


				if ((BLER_TABLE[mcs_index][scale_idx][1] >= re_scale_table[target_idx + 1])
					&& (BLER_TABLE[mcs_index][scale_idx + 1][1] <  re_scale_table[target_idx + 1])
					&& BLER_TABLE[mcs_index][scale_idx + 1][1] != 0)
				{
					scale_target[mcs_index][target_idx] = scale_idx;


				}
				else if (scale_idx == 0 && target_idx == 0 && BLER_TABLE[mcs_index][scale_idx][1] < re_scale_table[target_idx + 1])
				{
					scale_target[mcs_index][target_idx] = scale_idx;


				}
				
			}
		}
	}
	
	// NEW BLER table generator
	for (int mcs_index = 0; mcs_index < 38; mcs_index++)
	{
		for (int scale_idx = 0; scale_idx < 7; scale_idx++)
		{
			NEW_TABLE[mcs_index][scale_idx][1] = re_scale_table[scale_idx];
		}
		NEW_TABLE[mcs_index][0][0] = BLER_TABLE[mcs_index][0][0];
	}
	
	for (int mcs_index = 0; mcs_index < 38; mcs_index++)
	{
		for (int scale_idx = 1; scale_idx < 7; scale_idx++)
		{
			int t = scale_target[mcs_index][scale_idx - 1];

			//cout << "t = " << t << endl;
			
			//      cout << "scale_target["<< mcs_index <<"]["<< scale_idx - 1 <<"] = " << scale_target[mcs_index][ scale_idx - 1 ] << endl;
			

			NEW_TABLE[mcs_index][scale_idx][0] = log_linear_xais(BLER_TABLE[mcs_index][t][0], BLER_TABLE[mcs_index][t][1], BLER_TABLE[mcs_index][t + 1][0], BLER_TABLE[mcs_index][t + 1][1], NEW_TABLE[mcs_index][scale_idx][1]);
			
		}
	}
	
}

/*===================================================================
FUNCTION: table_setting_Threshold()

AUTHOR: SJ, KR, DW

DESCRIPTION:

NOTES:

===================================================================*/
void table_setting_Threshold()
{
	// Threshold computer 
	Real CR_for_Threshold[29];
	Real scale_sir[29] = { 0, };
	
	for (int i = 0; i < 29; i++)  // TS 36.213(Rel 8) - 7.1.7.1�� table ���� (Modulation and TBS index table for PDSCH)
	{
		int mod_type;
		int num_PRB = 50;   // Physical Resource Block
		int index_TBS = i;

		if (i <= 9)   
		{ mod_type = QPSK; }
		else if (i <= 16)   
		{ mod_type = QAM16; }
		else if (i <= 28)   
		{ mod_type = QAM64; }

		if (i > 9)     { index_TBS -= 1; }
		if (i > 16)     { index_TBS -= 1; }

		if (Effective_CR_compute(mod_type, index_TBS, num_PRB) < 1)
		{
			CR_for_Threshold[i] = Effective_CR_compute(mod_type, index_TBS, num_PRB);
		}
		else
		{
			CR_for_Threshold[i] = 1. / 3.;
			scale_sir[i] = linear2dB(Effective_CR_compute(mod_type, index_TBS, num_PRB));
		}

		//cout << " CR_for_Threshold[" << i <<"] = " << CR_for_Threshold[i] << endl;
		//getchar();

		if (mod_type == QPSK)
		for (int j = 0; j < 21; j++)
		{
			if (CR_for_Threshold[i] >= CR_curves[j] && CR_for_Threshold[i] <= CR_curves[j + 1])
			{
				Threshold[i] = linear_cr(CR_curves[j], NEW_TABLE[j][3][0], CR_curves[j + 1], NEW_TABLE[j + 1][3][0], CR_for_Threshold[i]) - scale_sir[i];
			}
		}

		if (mod_type == QAM16)
		for (int j = 21; j < 30; j++)
		{
			if (CR_for_Threshold[i] >= CR_curves[j] && CR_for_Threshold[i] <= CR_curves[j + 1])
			{
				Threshold[i] = linear_cr(CR_curves[j], NEW_TABLE[j][3][0], CR_curves[j + 1], NEW_TABLE[j + 1][3][0], CR_for_Threshold[i]) - scale_sir[i];
			}
		}

		if (mod_type == QAM64)
		for (int j = 30; j < 38; j++)
		{
			if (CR_for_Threshold[i] >= CR_curves[j] && CR_for_Threshold[i] <= CR_curves[j + 1])
			{
				Threshold[i] = linear_cr(CR_curves[j], NEW_TABLE[j][3][0], CR_curves[j + 1], NEW_TABLE[j + 1][3][0], CR_for_Threshold[i]) - scale_sir[i];
			}
		}
	}
	
	//Threshold[25] = Threshold[24] + 0.8 ;
	//Threshold[26] = Threshold[25] + 0.8 ;
	//Threshold[27] = Threshold[26] + 0.8 ;
	//Threshold[28] = Threshold[27] + 0.8 ;
	/*
	for( int i = 0 ; i < 29 ; i++ )
	{
	cout << " CR_for_Threshold[" << i <<"] = " << CR_for_Threshold[i] <<"\t"<< " Threshold["<< i <<"] = "<< Threshold[i] << endl;
	}

	getchar();
	*/



}


/*===================================================================
FUNCTION: Effective_CR_compute()

AUTHOR: SJ

DESCRIPTION: return the code rate from the standard if it is greater than 1/3
otherwise, return the effective SNR increase factor

NOTES:
===================================================================*/
Real Effective_CR_compute(int mod_type, int index_TBS, int num_PRB)
{
	Real effective_cr;

	int num_ofdm_symbols_for_traffic;
	num_ofdm_symbols_for_traffic = (int)(num_freq_per_rbs * num_ofdm_symbols_per_subband_per_1ms * num_PRB * (1. - overhead));


	if (mod_type == QPSK) { num_ofdm_symbols_for_traffic *= 2; }
	else if (mod_type == QAM16) { num_ofdm_symbols_for_traffic *= 4; }
	else if (mod_type == QAM64) { num_ofdm_symbols_for_traffic *= 6; }

	int num_coded_symbols;
	num_coded_symbols = inverse_TC_rate * (TBS[num_PRB][index_TBS] + num_crc);

	effective_cr = (((Real)num_coded_symbols) / ((Real)inverse_TC_rate)) / ((Real)num_ofdm_symbols_for_traffic);

	if (effective_cr > 0.916666667)
	{
		effective_cr = 0.916666667;
	}

	return(effective_cr);
}







/*===================================================================
FUNCTION: table_setting_BLER()

AUTHOR: SJ, KR, DW

DESCRIPTION:

NOTES:

===================================================================*/
void table_setting_BLER()
{
	file_name_setting();

	for (int mcs = 0; mcs < NUM_OF_MCS; mcs++)
	{
		ifstream BLER_table_file;

		if (mcs < 19)
		{
			BLER_table_file.open(bler_file_name[mcs]);
			for (int i = 0; i < 7; i++)
			{
				BLER_table_file >> BLER_TABLE[mcs][i][0];
				BLER_table_file >> BLER_TABLE[mcs][i][1];
			}
		}

		if (mcs == 19)
		{
			BLER_TABLE[mcs][0][0] = 3.98;  BLER_TABLE[mcs][0][1] = 1.00*pow(10, 0.);
			BLER_TABLE[mcs][1][0] = 4.18;  BLER_TABLE[mcs][1][1] = 9.40*pow(10, -1.);
			BLER_TABLE[mcs][2][0] = 4.38;  BLER_TABLE[mcs][2][1] = 3.98*pow(10, -1.);
			BLER_TABLE[mcs][3][0] = 4.58;  BLER_TABLE[mcs][3][1] = 3.97*pow(10, -2.);
			BLER_TABLE[mcs][4][0] = 4.78;  BLER_TABLE[mcs][4][1] = 3.30*pow(10, -3.);
			BLER_TABLE[mcs][5][0] = 0;  BLER_TABLE[mcs][5][1] = 0;
			BLER_TABLE[mcs][6][0] = 0;  BLER_TABLE[mcs][6][1] = 0;
		}


		if (mcs == 20)
		{
			BLER_TABLE[mcs][0][0] = 4.66;  BLER_TABLE[mcs][0][1] = 1.00*pow(10, 0.);
			BLER_TABLE[mcs][1][0] = 4.86;  BLER_TABLE[mcs][1][1] = 9.94*pow(10, -1.);
			BLER_TABLE[mcs][2][0] = 5.06;  BLER_TABLE[mcs][2][1] = 7.28*pow(10, -1.);
			BLER_TABLE[mcs][3][0] = 5.26;  BLER_TABLE[mcs][3][1] = 1.38*pow(10, -1.);
			BLER_TABLE[mcs][4][0] = 5.46;  BLER_TABLE[mcs][4][1] = 4.97*pow(10, -3.);
			BLER_TABLE[mcs][5][0] = 0;  BLER_TABLE[mcs][5][1] = 0;
			BLER_TABLE[mcs][6][0] = 0;  BLER_TABLE[mcs][6][1] = 0;
		}

		if (mcs == 21)
		{
			BLER_TABLE[mcs][0][0] = 3.06;  BLER_TABLE[mcs][0][1] = 1.00*pow(10, 0.);
			BLER_TABLE[mcs][1][0] = 3.26;  BLER_TABLE[mcs][1][1] = 9.14*pow(10, -1.);
			BLER_TABLE[mcs][2][0] = 3.46;  BLER_TABLE[mcs][2][1] = 2.58*pow(10, -1.);
			BLER_TABLE[mcs][3][0] = 3.56;  BLER_TABLE[mcs][3][1] = 5.72*pow(10, -2.);
			BLER_TABLE[mcs][4][0] = 3.66;  BLER_TABLE[mcs][4][1] = 7.15*pow(10, -3.);
			BLER_TABLE[mcs][5][0] = 0;  BLER_TABLE[mcs][5][1] = 0;
			BLER_TABLE[mcs][6][0] = 0;  BLER_TABLE[mcs][6][1] = 0;
		}

		if (21 < mcs && mcs < 27)
		{
			BLER_table_file.open(bler_file_name[mcs - 3]);
			for (int i = 0; i < 7; i++)
			{
				BLER_table_file >> BLER_TABLE[mcs][i][0];
				BLER_table_file >> BLER_TABLE[mcs][i][1];
			}
		}

		if (mcs == 27)
		{
			BLER_TABLE[mcs][0][0] = 8.47;  BLER_TABLE[mcs][0][1] = 1.00*pow(10, 0.);
			BLER_TABLE[mcs][1][0] = 8.67;  BLER_TABLE[mcs][1][1] = 9.92*pow(10, -1.);
			BLER_TABLE[mcs][2][0] = 8.87;  BLER_TABLE[mcs][2][1] = 6.67*pow(10, -1.);
			BLER_TABLE[mcs][3][0] = 9.07;  BLER_TABLE[mcs][3][1] = 1.08*pow(10, -1.);
			BLER_TABLE[mcs][4][0] = 9.27;  BLER_TABLE[mcs][4][1] = 1.11*pow(10, -2.);
			BLER_TABLE[mcs][5][0] = 9.37;  BLER_TABLE[mcs][5][1] = 3.80*pow(10, -3.);
			BLER_TABLE[mcs][6][0] = 0;  BLER_TABLE[mcs][6][1] = 0;
		}

		if (mcs == 28)
		{
			BLER_TABLE[mcs][0][0] = 10.18;  BLER_TABLE[mcs][0][1] = 1.00*pow(10, 0.);
			BLER_TABLE[mcs][1][0] = 10.38;  BLER_TABLE[mcs][1][1] = 8.95*pow(10, -1.);
			BLER_TABLE[mcs][2][0] = 10.58;  BLER_TABLE[mcs][2][1] = 2.79*pow(10, -1.);
			BLER_TABLE[mcs][3][0] = 10.78;  BLER_TABLE[mcs][3][1] = 2.00*pow(10, -2.);
			BLER_TABLE[mcs][4][0] = 10.98;  BLER_TABLE[mcs][4][1] = 1.57*pow(10, -3.);
			BLER_TABLE[mcs][5][0] = 0;  BLER_TABLE[mcs][5][1] = 0;
			BLER_TABLE[mcs][6][0] = 0;  BLER_TABLE[mcs][6][1] = 0;
		}

		if (mcs == 29)
		{
			BLER_TABLE[mcs][0][0] = 11.07;  BLER_TABLE[mcs][0][1] = 1.00*pow(10, 0.);
			BLER_TABLE[mcs][1][0] = 11.27;  BLER_TABLE[mcs][1][1] = 9.51*pow(10, -1.);
			BLER_TABLE[mcs][2][0] = 11.47;  BLER_TABLE[mcs][2][1] = 3.60*pow(10, -1.);
			BLER_TABLE[mcs][3][0] = 11.67;  BLER_TABLE[mcs][3][1] = 2.42*pow(10, -2.);
			BLER_TABLE[mcs][4][0] = 11.77;  BLER_TABLE[mcs][4][1] = 3.30*pow(10, -3.);
			BLER_TABLE[mcs][5][0] = 0;  BLER_TABLE[mcs][5][1] = 0;
			BLER_TABLE[mcs][6][0] = 0;  BLER_TABLE[mcs][6][1] = 0;
		}

		if (mcs == 30)
		{
			BLER_TABLE[mcs][0][0] = 6.20;  BLER_TABLE[mcs][0][1] = 1.00*pow(10, 0.);
			BLER_TABLE[mcs][1][0] = 6.40;  BLER_TABLE[mcs][1][1] = 8.62*pow(10, -1.);
			BLER_TABLE[mcs][2][0] = 6.60;  BLER_TABLE[mcs][2][1] = 3.85*pow(10, -1.);
			BLER_TABLE[mcs][3][0] = 6.80;  BLER_TABLE[mcs][3][1] = 5.31*pow(10, -2.);
			BLER_TABLE[mcs][4][0] = 7.00;  BLER_TABLE[mcs][4][1] = 1.80*pow(10, -3.);
			BLER_TABLE[mcs][5][0] = 0;  BLER_TABLE[mcs][5][1] = 0;
			BLER_TABLE[mcs][6][0] = 0;  BLER_TABLE[mcs][6][1] = 0;
		}

		if (mcs >  30)
		{
			BLER_table_file.open(bler_file_name[mcs - 7]);
			for (int i = 0; i < 7; i++)
			{
				BLER_table_file >> BLER_TABLE[mcs][i][0];
				BLER_table_file >> BLER_TABLE[mcs][i][1];

				//cout << BLER_TABLE[mcs][i][0] << endl;
				//cout << BLER_TABLE[mcs][i][1] << endl;
			}
		}
	}
}

/*===================================================================
FUNCTION: Declare_Channel_Output

2009-12-04
===================================================================*/
void file_name_setting(void)
{
	//bler_file_name = new char *[31];   
	bler_file_name = new char *[38];
	for (int mcs = 0; mcs < NUM_OF_MCS; mcs++)
	{
		bler_file_name[mcs] = new char[70];
	}

	//for (int mcs = 0; mcs < 31; mcs++)
	for (int mcs = 0; mcs < 38; mcs++)
	{
		if (mcs < 10)
		{
			sprintf(bler_file_name[mcs], "./LLR_curv/Stream2-32LRU(rank-1)/MCS%01d_32RU.dat", mcs);
		}
		else
		{
			sprintf(bler_file_name[mcs], "./LLR_curv/Stream2-32LRU(rank-1)/MCS%02d_32RU.dat", mcs);
		}

		//cout << bler_file_name[mcs] << endl;
		//getchar();
	}

	//cout << bler_file_name[0] << endl;
	//cout << bler_file_name[1] << endl; 
	//cout << bler_file_name[2] << endl;

}



/*===================================================================
FUNCTION: table_setting_TBS()

AUTHOR: SJ, KR, DW

DESCRIPTION:

NOTES:

===================================================================*/
void table_setting_TBS()
{
	int n_rbs = 0;
	for (int i_tbs = 0; i_tbs < 27; i_tbs++)
	{
		TBS[n_rbs][i_tbs] = 0;
	}

	TBS[1][0] = 16; TBS[2][0] = 32; TBS[3][0] = 56; TBS[4][0] = 88; TBS[5][0] = 120;
	TBS[6][0] = 152; TBS[7][0] = 176; TBS[8][0] = 200; TBS[9][0] = 232; TBS[10][0] = 248;
	TBS[11][0] = 288; TBS[12][0] = 304; TBS[13][0] = 344; TBS[14][0] = 376; TBS[15][0] = 392;
	TBS[16][0] = 424; TBS[17][0] = 456; TBS[18][0] = 488; TBS[19][0] = 504; TBS[20][0] = 536;
	TBS[21][0] = 568; TBS[22][0] = 600; TBS[23][0] = 616; TBS[24][0] = 648; TBS[25][0] = 680;
	TBS[26][0] = 712; TBS[27][0] = 744; TBS[28][0] = 776; TBS[29][0] = 776; TBS[30][0] = 808;
	TBS[31][0] = 840; TBS[32][0] = 872; TBS[33][0] = 904; TBS[34][0] = 936; TBS[35][0] = 968;
	TBS[36][0] = 1000; TBS[37][0] = 1032; TBS[38][0] = 1032; TBS[39][0] = 1064; TBS[40][0] = 1096;
	TBS[41][0] = 1128; TBS[42][0] = 1160; TBS[43][0] = 1192; TBS[44][0] = 1224; TBS[45][0] = 1256;
	TBS[46][0] = 1256; TBS[47][0] = 1288; TBS[48][0] = 1320; TBS[49][0] = 1352; TBS[50][0] = 1384;

	TBS[1][1] = 24; TBS[2][1] = 48; TBS[3][1] = 88; TBS[4][1] = 120; TBS[5][1] = 160;
	TBS[6][1] = 200; TBS[7][1] = 232; TBS[8][1] = 272; TBS[9][1] = 304; TBS[10][1] = 344;
	TBS[11][1] = 376; TBS[12][1] = 424; TBS[13][1] = 456; TBS[14][1] = 488; TBS[15][1] = 520;
	TBS[16][1] = 568; TBS[17][1] = 600; TBS[18][1] = 632; TBS[19][1] = 680; TBS[20][1] = 712;
	TBS[21][1] = 744; TBS[22][1] = 776; TBS[23][1] = 808; TBS[24][1] = 872; TBS[25][1] = 904;
	TBS[26][1] = 936; TBS[27][1] = 968; TBS[28][1] = 1000; TBS[29][1] = 1032; TBS[30][1] = 1064;
	TBS[31][1] = 1128; TBS[32][1] = 1160; TBS[33][1] = 1192; TBS[34][1] = 1224; TBS[35][1] = 1256;
	TBS[36][1] = 1288; TBS[37][1] = 1352; TBS[38][1] = 1384; TBS[39][1] = 1416; TBS[40][1] = 1416;
	TBS[41][1] = 1480; TBS[42][1] = 1544; TBS[43][1] = 1544; TBS[44][1] = 1608; TBS[45][1] = 1608;
	TBS[46][1] = 1672; TBS[47][1] = 1736; TBS[48][1] = 1736; TBS[49][1] = 1800; TBS[50][1] = 1800;

	TBS[1][2] = 32; TBS[2][2] = 72; TBS[3][2] = 120; TBS[4][2] = 160; TBS[5][2] = 200;
	TBS[6][2] = 248; TBS[7][2] = 296; TBS[8][2] = 336; TBS[9][2] = 376; TBS[10][2] = 424;
	TBS[11][2] = 472; TBS[12][2] = 520; TBS[13][2] = 568; TBS[14][2] = 616; TBS[15][2] = 648;
	TBS[16][2] = 696; TBS[17][2] = 744; TBS[18][2] = 776; TBS[19][2] = 840; TBS[20][2] = 872;
	TBS[21][2] = 936; TBS[22][2] = 968; TBS[23][2] = 1000; TBS[24][2] = 1064; TBS[25][2] = 1096;
	TBS[26][2] = 1160; TBS[27][2] = 1192; TBS[28][2] = 1256; TBS[29][2] = 1288; TBS[30][2] = 1320;
	TBS[31][2] = 1384; TBS[32][2] = 1416; TBS[33][2] = 1480; TBS[34][2] = 1544; TBS[35][2] = 1544;
	TBS[36][2] = 1608; TBS[37][2] = 1672; TBS[38][2] = 1672; TBS[39][2] = 1736; TBS[40][2] = 1800;
	TBS[41][2] = 1800; TBS[42][2] = 1864; TBS[43][2] = 1928; TBS[44][2] = 1992; TBS[45][2] = 2024;
	TBS[46][2] = 2088; TBS[47][2] = 2088; TBS[48][2] = 2152; TBS[49][2] = 2216; TBS[50][2] = 2216;

	TBS[1][3] = 40; TBS[2][3] = 104; TBS[3][3] = 152; TBS[4][3] = 208; TBS[5][3] = 272;
	TBS[6][3] = 320; TBS[7][3] = 392; TBS[8][3] = 440; TBS[9][3] = 504; TBS[10][3] = 568;
	TBS[11][3] = 616; TBS[12][3] = 680; TBS[13][3] = 744; TBS[14][3] = 808; TBS[15][3] = 872;
	TBS[16][3] = 904; TBS[17][3] = 968; TBS[18][3] = 1032; TBS[19][3] = 1096; TBS[20][3] = 1160;
	TBS[21][3] = 1224; TBS[22][3] = 1256; TBS[23][3] = 1320; TBS[24][3] = 1384; TBS[25][3] = 1416;
	TBS[26][3] = 1480; TBS[27][3] = 1544; TBS[28][3] = 1608; TBS[29][3] = 1672; TBS[30][3] = 1736;
	TBS[31][3] = 1800; TBS[32][3] = 1864; TBS[33][3] = 1928; TBS[34][3] = 1992; TBS[35][3] = 2024;
	TBS[36][3] = 2088; TBS[37][3] = 2152; TBS[38][3] = 2216; TBS[39][3] = 2280; TBS[40][3] = 2344;
	TBS[41][3] = 2408; TBS[42][3] = 2472; TBS[43][3] = 2536; TBS[44][3] = 2536; TBS[45][3] = 2600;
	TBS[46][3] = 2664; TBS[47][3] = 2728; TBS[48][3] = 2792; TBS[49][3] = 2856; TBS[50][3] = 2856;

	TBS[1][4] = 48; TBS[2][4] = 120; TBS[3][4] = 200; TBS[4][4] = 264; TBS[5][4] = 320;
	TBS[6][4] = 408; TBS[7][4] = 488; TBS[8][4] = 552; TBS[9][4] = 632; TBS[10][4] = 696;
	TBS[11][4] = 776; TBS[12][4] = 840; TBS[13][4] = 904; TBS[14][4] = 1000; TBS[15][4] = 1064;
	TBS[16][4] = 1128; TBS[17][4] = 1192; TBS[18][4] = 1288; TBS[19][4] = 1352; TBS[20][4] = 1416;
	TBS[21][4] = 1480; TBS[22][4] = 1544; TBS[23][4] = 1608; TBS[24][4] = 1736; TBS[25][4] = 1800;
	TBS[26][4] = 1864; TBS[27][4] = 1928; TBS[28][4] = 1992; TBS[29][4] = 2088; TBS[30][4] = 2152;
	TBS[31][4] = 2216; TBS[32][4] = 2280; TBS[33][4] = 2344; TBS[34][4] = 2408; TBS[35][4] = 2472;
	TBS[36][4] = 2600; TBS[37][4] = 2664; TBS[38][4] = 2728; TBS[39][4] = 2792; TBS[40][4] = 2856;
	TBS[41][4] = 2984; TBS[42][4] = 2984; TBS[43][4] = 3112; TBS[44][4] = 3112; TBS[45][4] = 3240;
	TBS[46][4] = 3240; TBS[47][4] = 3368; TBS[48][4] = 3496; TBS[49][4] = 3496; TBS[50][4] = 3624;

	TBS[1][5] = 72; TBS[2][5] = 152; TBS[3][5] = 232; TBS[4][5] = 320; TBS[5][5] = 424;
	TBS[6][5] = 504; TBS[7][5] = 600; TBS[8][5] = 680; TBS[9][5] = 776; TBS[10][5] = 872;
	TBS[11][5] = 968; TBS[12][5] = 1032; TBS[13][5] = 1128; TBS[14][5] = 1224; TBS[15][5] = 1320;
	TBS[16][5] = 1384; TBS[17][5] = 1480; TBS[18][5] = 1544; TBS[19][5] = 1672; TBS[20][5] = 1736;
	TBS[21][5] = 1864; TBS[22][5] = 1928; TBS[23][5] = 2024; TBS[24][5] = 2088; TBS[25][5] = 2216;
	TBS[26][5] = 2280; TBS[27][5] = 2344; TBS[28][5] = 2472; TBS[29][5] = 2536; TBS[30][5] = 2664;
	TBS[31][5] = 2728; TBS[32][5] = 2792; TBS[33][5] = 2856; TBS[34][5] = 2984; TBS[35][5] = 3112;
	TBS[36][5] = 3112; TBS[37][5] = 3240; TBS[38][5] = 3368; TBS[39][5] = 3496; TBS[40][5] = 3496;
	TBS[41][5] = 3624; TBS[42][5] = 3752; TBS[43][5] = 3752; TBS[44][5] = 3880; TBS[45][5] = 4008;
	TBS[46][5] = 4008; TBS[47][5] = 4136; TBS[48][5] = 4264; TBS[49][5] = 4392; TBS[50][5] = 4392;

	TBS[1][6] = 88; TBS[2][6] = 176; TBS[3][6] = 288; TBS[4][6] = 392; TBS[5][6] = 504;
	TBS[6][6] = 600; TBS[7][6] = 712; TBS[8][6] = 808; TBS[9][6] = 936; TBS[10][6] = 1032;
	TBS[11][6] = 1128; TBS[12][6] = 1224; TBS[13][6] = 1352; TBS[14][6] = 1480; TBS[15][6] = 1544;
	TBS[16][6] = 1672; TBS[17][6] = 1736; TBS[18][6] = 1864; TBS[19][6] = 1992; TBS[20][6] = 2088;
	TBS[21][6] = 2216; TBS[22][6] = 2280; TBS[23][6] = 2408; TBS[24][6] = 2472; TBS[25][6] = 2600;
	TBS[26][6] = 2728; TBS[27][6] = 2792; TBS[28][6] = 2984; TBS[29][6] = 2984; TBS[30][6] = 3112;
	TBS[31][6] = 3240; TBS[32][6] = 3368; TBS[33][6] = 3496; TBS[34][6] = 3496; TBS[35][6] = 3624;
	TBS[36][6] = 3752; TBS[37][6] = 3880; TBS[38][6] = 4008; TBS[39][6] = 4136; TBS[40][6] = 4136;
	TBS[41][6] = 4264; TBS[42][6] = 4392; TBS[43][6] = 4584; TBS[44][6] = 4584; TBS[45][6] = 4776;
	TBS[46][6] = 4776; TBS[47][6] = 4968; TBS[48][6] = 4968; TBS[49][6] = 5160; TBS[50][6] = 5160;

	TBS[1][7] = 104; TBS[2][7] = 232; TBS[3][7] = 320; TBS[4][7] = 472; TBS[5][7] = 584;
	TBS[6][7] = 712; TBS[7][7] = 840; TBS[8][7] = 968; TBS[9][7] = 1096; TBS[10][7] = 1224;
	TBS[11][7] = 1320; TBS[12][7] = 1480; TBS[13][7] = 1608; TBS[14][7] = 1672; TBS[15][7] = 1800;
	TBS[16][7] = 1928; TBS[17][7] = 2088; TBS[18][7] = 2216; TBS[19][7] = 2344; TBS[20][7] = 2472;
	TBS[21][7] = 2536; TBS[22][7] = 2664; TBS[23][7] = 2792; TBS[24][7] = 2984; TBS[25][7] = 3112;
	TBS[26][7] = 3240; TBS[27][7] = 3368; TBS[28][7] = 3368; TBS[29][7] = 3496; TBS[30][7] = 3624;
	TBS[31][7] = 3752; TBS[32][7] = 3880; TBS[33][7] = 4008; TBS[34][7] = 4136; TBS[35][7] = 4264;
	TBS[36][7] = 4392; TBS[37][7] = 4584; TBS[38][7] = 4584; TBS[39][7] = 4776; TBS[40][7] = 4968;
	TBS[41][7] = 4968; TBS[42][7] = 5160; TBS[43][7] = 5352; TBS[44][7] = 5352; TBS[45][7] = 5544;
	TBS[46][7] = 5736; TBS[47][7] = 5736; TBS[48][7] = 5992; TBS[49][7] = 5992; TBS[50][7] = 6200;

	TBS[1][8] = 120; TBS[2][8] = 248; TBS[3][8] = 392; TBS[4][8] = 536; TBS[5][8] = 680;
	TBS[6][8] = 808; TBS[7][8] = 968; TBS[8][8] = 1096; TBS[9][8] = 1256; TBS[10][8] = 1384;
	TBS[11][8] = 1544; TBS[12][8] = 1672; TBS[13][8] = 1800; TBS[14][8] = 1928; TBS[15][8] = 2088;
	TBS[16][8] = 2216; TBS[17][8] = 2344; TBS[18][8] = 2536; TBS[19][8] = 2664; TBS[20][8] = 2792;
	TBS[21][8] = 2984; TBS[22][8] = 3112; TBS[23][8] = 3240; TBS[24][8] = 3368; TBS[25][8] = 3496;
	TBS[26][8] = 3624; TBS[27][8] = 3752; TBS[28][8] = 3880; TBS[29][8] = 4008; TBS[30][8] = 4264;
	TBS[31][8] = 4392; TBS[32][8] = 4584; TBS[33][8] = 4584; TBS[34][8] = 4776; TBS[35][8] = 4968;
	TBS[36][8] = 4968; TBS[37][8] = 5160; TBS[38][8] = 5352; TBS[39][8] = 5544; TBS[40][8] = 5544;
	TBS[41][8] = 5736; TBS[42][8] = 5992; TBS[43][8] = 5992; TBS[44][8] = 6200; TBS[45][8] = 6200;
	TBS[46][8] = 6456; TBS[47][8] = 6456; TBS[48][8] = 6712; TBS[49][8] = 6968; TBS[50][8] = 6968;

	TBS[1][9] = 136; TBS[2][9] = 296; TBS[3][9] = 456; TBS[4][9] = 616; TBS[5][9] = 776;
	TBS[6][9] = 936; TBS[7][9] = 1096; TBS[8][9] = 1256; TBS[9][9] = 1416; TBS[10][9] = 1544;
	TBS[11][9] = 1736; TBS[12][9] = 1864; TBS[13][9] = 2024; TBS[14][9] = 2216; TBS[15][9] = 2344;
	TBS[16][9] = 2536; TBS[17][9] = 2664; TBS[18][9] = 2856; TBS[19][9] = 2984; TBS[20][9] = 3112;
	TBS[21][9] = 3368; TBS[22][9] = 3496; TBS[23][9] = 3624; TBS[24][9] = 3752; TBS[25][9] = 4008;
	TBS[26][9] = 4136; TBS[27][9] = 4264; TBS[28][9] = 4392; TBS[29][9] = 4584; TBS[30][9] = 4776;
	TBS[31][9] = 4968; TBS[32][9] = 5160; TBS[33][9] = 5160; TBS[34][9] = 5352; TBS[35][9] = 5544;
	TBS[36][9] = 5736; TBS[37][9] = 5736; TBS[38][9] = 5992; TBS[39][9] = 6200; TBS[40][9] = 6200;
	TBS[41][9] = 6456; TBS[42][9] = 6712; TBS[43][9] = 6712; TBS[44][9] = 6968; TBS[45][9] = 6968;
	TBS[46][9] = 7224; TBS[47][9] = 7480; TBS[48][9] = 7480; TBS[49][9] = 7736; TBS[50][9] = 7992;

	TBS[1][10] = 152; TBS[2][10] = 320; TBS[3][10] = 504; TBS[4][10] = 680; TBS[5][10] = 872;
	TBS[6][10] = 1032; TBS[7][10] = 1224; TBS[8][10] = 1384; TBS[9][10] = 1544; TBS[10][10] = 1736;
	TBS[11][10] = 1928; TBS[12][10] = 2088; TBS[13][10] = 2280; TBS[14][10] = 2472; TBS[15][10] = 2664;
	TBS[16][10] = 2792; TBS[17][10] = 2984; TBS[18][10] = 3112; TBS[19][10] = 3368; TBS[20][10] = 3496;
	TBS[21][10] = 3752; TBS[22][10] = 3880; TBS[23][10] = 4008; TBS[24][10] = 4264; TBS[25][10] = 4392;
	TBS[26][10] = 4584; TBS[27][10] = 4776; TBS[28][10] = 4968; TBS[29][10] = 5160; TBS[30][10] = 5352;
	TBS[31][10] = 5544; TBS[32][10] = 5736; TBS[33][10] = 5736; TBS[34][10] = 5992; TBS[35][10] = 6200;
	TBS[36][10] = 6200; TBS[37][10] = 6456; TBS[38][10] = 6712; TBS[39][10] = 6712; TBS[40][10] = 6968;
	TBS[41][10] = 7224; TBS[42][10] = 7480; TBS[43][10] = 7480; TBS[44][10] = 7736; TBS[45][10] = 7992;
	TBS[46][10] = 7992; TBS[47][10] = 8248; TBS[48][10] = 8504; TBS[49][10] = 8504; TBS[50][10] = 8760;

	TBS[1][11] = 176; TBS[2][11] = 376; TBS[3][11] = 584; TBS[4][11] = 776; TBS[5][11] = 1000;
	TBS[6][11] = 1192; TBS[7][11] = 1384; TBS[8][11] = 1608; TBS[9][11] = 1800; TBS[10][11] = 2024;
	TBS[11][11] = 2216; TBS[12][11] = 2408; TBS[13][11] = 2600; TBS[14][11] = 2792; TBS[15][11] = 2984;
	TBS[16][11] = 3240; TBS[17][11] = 3496; TBS[18][11] = 3624; TBS[19][11] = 3880; TBS[20][11] = 4008;
	TBS[21][11] = 4264; TBS[22][11] = 4392; TBS[23][11] = 4584; TBS[24][11] = 4776; TBS[25][11] = 4968;
	TBS[26][11] = 5352; TBS[27][11] = 5544; TBS[28][11] = 5736; TBS[29][11] = 5992; TBS[30][11] = 5992;
	TBS[31][11] = 6200; TBS[32][11] = 6456; TBS[33][11] = 6712; TBS[34][11] = 6968; TBS[35][11] = 6968;
	TBS[36][11] = 7224; TBS[37][11] = 7480; TBS[38][11] = 7736; TBS[39][11] = 7736; TBS[40][11] = 7992;
	TBS[41][11] = 8248; TBS[42][11] = 8504; TBS[43][11] = 8760; TBS[44][11] = 8760; TBS[45][11] = 9144;
	TBS[46][11] = 9144; TBS[47][11] = 9528; TBS[48][11] = 9528; TBS[49][11] = 9912; TBS[50][11] = 9912;

	TBS[1][12] = 208; TBS[2][12] = 440; TBS[3][12] = 680; TBS[4][12] = 904; TBS[5][12] = 1128;
	TBS[6][12] = 1352; TBS[7][12] = 1608; TBS[8][12] = 1800; TBS[9][12] = 2024; TBS[10][12] = 2280;
	TBS[11][12] = 2472; TBS[12][12] = 2728; TBS[13][12] = 2984; TBS[14][12] = 3240; TBS[15][12] = 3368;
	TBS[16][12] = 3624; TBS[17][12] = 3880; TBS[18][12] = 4136; TBS[19][12] = 4392; TBS[20][12] = 4584;
	TBS[21][12] = 4776; TBS[22][12] = 4968; TBS[23][12] = 5352; TBS[24][12] = 5544; TBS[25][12] = 5736;
	TBS[26][12] = 5992; TBS[27][12] = 6200; TBS[28][12] = 6456; TBS[29][12] = 6712; TBS[30][12] = 6712;
	TBS[31][12] = 6968; TBS[32][12] = 7224; TBS[33][12] = 7480; TBS[34][12] = 7736; TBS[35][12] = 7992;
	TBS[36][12] = 8248; TBS[37][12] = 8504; TBS[38][12] = 8760; TBS[39][12] = 8760; TBS[40][12] = 9144;
	TBS[41][12] = 9528; TBS[42][12] = 9528; TBS[43][12] = 9912; TBS[44][12] = 9912; TBS[45][12] = 10296;
	TBS[46][12] = 10680; TBS[47][12] = 10680; TBS[48][12] = 11064; TBS[49][12] = 11064; TBS[50][12] = 11448;

	TBS[1][13] = 232; TBS[2][13] = 488; TBS[3][13] = 744; TBS[4][13] = 1000; TBS[5][13] = 1256;
	TBS[6][13] = 1544; TBS[7][13] = 1800; TBS[8][13] = 2024; TBS[9][13] = 2280; TBS[10][13] = 2536;
	TBS[11][13] = 2856; TBS[12][13] = 3112; TBS[13][13] = 3368; TBS[14][13] = 3624; TBS[15][13] = 3880;
	TBS[16][13] = 4136; TBS[17][13] = 4392; TBS[18][13] = 4584; TBS[19][13] = 4968; TBS[20][13] = 5160;
	TBS[21][13] = 5352; TBS[22][13] = 5736; TBS[23][13] = 5992; TBS[24][13] = 6200; TBS[25][13] = 6456;
	TBS[26][13] = 6712; TBS[27][13] = 6968; TBS[28][13] = 7224; TBS[29][13] = 7480; TBS[30][13] = 7736;
	TBS[31][13] = 7992; TBS[32][13] = 8248; TBS[33][13] = 8504; TBS[34][13] = 8760; TBS[35][13] = 9144;
	TBS[36][13] = 9144; TBS[37][13] = 9528; TBS[38][13] = 9912; TBS[39][13] = 9912; TBS[40][13] = 10296;
	TBS[41][13] = 10680; TBS[42][13] = 10680; TBS[43][13] = 11064; TBS[44][13] = 11448; TBS[45][13] = 11448;
	TBS[46][13] = 11832; TBS[47][13] = 12216; TBS[48][13] = 12216; TBS[49][13] = 12576; TBS[50][13] = 12960;

	TBS[1][14] = 264; TBS[2][14] = 552; TBS[3][14] = 840; TBS[4][14] = 1128; TBS[5][14] = 1416;
	TBS[6][14] = 1736; TBS[7][14] = 1992; TBS[8][14] = 2280; TBS[9][14] = 2600; TBS[10][14] = 2856;
	TBS[11][14] = 3112; TBS[12][14] = 3496; TBS[13][14] = 3752; TBS[14][14] = 4008; TBS[15][14] = 4264;
	TBS[16][14] = 4584; TBS[17][14] = 4968; TBS[18][14] = 5160; TBS[19][14] = 5544; TBS[20][14] = 5736;
	TBS[21][14] = 5992; TBS[22][14] = 6200; TBS[23][14] = 6456; TBS[24][14] = 6968; TBS[25][14] = 7224;
	TBS[26][14] = 7480; TBS[27][14] = 7736; TBS[28][14] = 7992; TBS[29][14] = 8248; TBS[30][14] = 8504;
	TBS[31][14] = 8760; TBS[32][14] = 9144; TBS[33][14] = 9528; TBS[34][14] = 9912; TBS[35][14] = 9912;
	TBS[36][14] = 10296; TBS[37][14] = 10680; TBS[38][14] = 11064; TBS[39][14] = 11064; TBS[40][14] = 11448;
	TBS[41][14] = 11832; TBS[42][14] = 12216; TBS[43][14] = 12216; TBS[44][14] = 12576; TBS[45][14] = 12960;
	TBS[46][14] = 12960; TBS[47][14] = 13536; TBS[48][14] = 13536; TBS[49][14] = 14112; TBS[50][14] = 14112;

	TBS[1][15] = 280; TBS[2][15] = 600; TBS[3][15] = 904; TBS[4][15] = 1224; TBS[5][15] = 1544;
	TBS[6][15] = 1800; TBS[7][15] = 2152; TBS[8][15] = 2472; TBS[9][15] = 2728; TBS[10][15] = 3112;
	TBS[11][15] = 3368; TBS[12][15] = 3624; TBS[13][15] = 4008; TBS[14][15] = 4264; TBS[15][15] = 4584;
	TBS[16][15] = 4968; TBS[17][15] = 5160; TBS[18][15] = 5544; TBS[19][15] = 5736; TBS[20][15] = 6200;
	TBS[21][15] = 6456; TBS[22][15] = 6712; TBS[23][15] = 6968; TBS[24][15] = 7224; TBS[25][15] = 7736;
	TBS[26][15] = 7992; TBS[27][15] = 8248; TBS[28][15] = 8504; TBS[29][15] = 8760; TBS[30][15] = 9144;
	TBS[31][15] = 9528; TBS[32][15] = 9912; TBS[33][15] = 10296; TBS[34][15] = 10296; TBS[35][15] = 10680;
	TBS[36][15] = 11064; TBS[37][15] = 11448; TBS[38][15] = 11832; TBS[39][15] = 11832; TBS[40][15] = 12216;
	TBS[41][15] = 12576; TBS[42][15] = 12960; TBS[43][15] = 12960; TBS[44][15] = 13536; TBS[45][15] = 13536;
	TBS[46][15] = 14112; TBS[47][15] = 14688; TBS[48][15] = 14688; TBS[49][15] = 15264; TBS[50][15] = 15264;

	TBS[1][16] = 320; TBS[2][16] = 632; TBS[3][16] = 968; TBS[4][16] = 1288; TBS[5][16] = 1608;
	TBS[6][16] = 1928; TBS[7][16] = 2280; TBS[8][16] = 2600; TBS[9][16] = 2984; TBS[10][16] = 3240;
	TBS[11][16] = 3624; TBS[12][16] = 3880; TBS[13][16] = 4264; TBS[14][16] = 4584; TBS[15][16] = 4968;
	TBS[16][16] = 5160; TBS[17][16] = 5544; TBS[18][16] = 5992; TBS[19][16] = 6200; TBS[20][16] = 6456;
	TBS[21][16] = 6712; TBS[22][16] = 7224; TBS[23][16] = 7480; TBS[24][16] = 7736; TBS[25][16] = 7992;
	TBS[26][16] = 8504; TBS[27][16] = 8760; TBS[28][16] = 9144; TBS[29][16] = 9528; TBS[30][16] = 9912;
	TBS[31][16] = 9912; TBS[32][16] = 10296; TBS[33][16] = 10680; TBS[34][16] = 11064; TBS[35][16] = 11448;
	TBS[36][16] = 11832; TBS[37][16] = 12216; TBS[38][16] = 12216; TBS[39][16] = 12576; TBS[40][16] = 12960;
	TBS[41][16] = 13536; TBS[42][16] = 13536; TBS[43][16] = 14112; TBS[44][16] = 14112; TBS[45][16] = 14688;
	TBS[46][16] = 14688; TBS[47][16] = 15264; TBS[48][16] = 15840; TBS[49][16] = 15840; TBS[50][16] = 16416;

	TBS[1][17] = 336; TBS[2][17] = 696; TBS[3][17] = 1064; TBS[4][17] = 1416; TBS[5][17] = 1800;
	TBS[6][17] = 2152; TBS[7][17] = 2536; TBS[8][17] = 2856; TBS[9][17] = 3240; TBS[10][17] = 3624;
	TBS[11][17] = 4008; TBS[12][17] = 4392; TBS[13][17] = 4776; TBS[14][17] = 5160; TBS[15][17] = 5352;
	TBS[16][17] = 5736; TBS[17][17] = 6200; TBS[18][17] = 6456; TBS[19][17] = 6712; TBS[20][17] = 7224;
	TBS[21][17] = 7480; TBS[22][17] = 7992; TBS[23][17] = 8248; TBS[24][17] = 8760; TBS[25][17] = 9144;
	TBS[26][17] = 9528; TBS[27][17] = 9912; TBS[28][17] = 10296; TBS[29][17] = 10296; TBS[30][17] = 10680;
	TBS[31][17] = 11064; TBS[32][17] = 11448; TBS[33][17] = 11832; TBS[34][17] = 12216; TBS[35][17] = 12576;
	TBS[36][17] = 12960; TBS[37][17] = 13536; TBS[38][17] = 13536; TBS[39][17] = 14112; TBS[40][17] = 14688;
	TBS[41][17] = 14688; TBS[42][17] = 15264; TBS[43][17] = 15264; TBS[44][17] = 15840; TBS[45][17] = 16416;
	TBS[46][17] = 16416; TBS[47][17] = 16992; TBS[48][17] = 17568; TBS[49][17] = 17568; TBS[50][17] = 18336;

	TBS[1][18] = 376; TBS[2][18] = 776; TBS[3][18] = 1160; TBS[4][18] = 1544; TBS[5][18] = 1992;
	TBS[6][18] = 2344; TBS[7][18] = 2792; TBS[8][18] = 3112; TBS[9][18] = 3624; TBS[10][18] = 4008;
	TBS[11][18] = 4392; TBS[12][18] = 4776; TBS[13][18] = 5160; TBS[14][18] = 5544; TBS[15][18] = 5992;
	TBS[16][18] = 6200; TBS[17][18] = 6712; TBS[18][18] = 7224; TBS[19][18] = 7480; TBS[20][18] = 7992;
	TBS[21][18] = 8248; TBS[22][18] = 8760; TBS[23][18] = 9144; TBS[24][18] = 9528; TBS[25][18] = 9912;
	TBS[26][18] = 10296; TBS[27][18] = 10680; TBS[28][18] = 11064; TBS[29][18] = 11448; TBS[30][18] = 11832;
	TBS[31][18] = 12216; TBS[32][18] = 12576; TBS[33][18] = 12960; TBS[34][18] = 13536; TBS[35][18] = 14112;
	TBS[36][18] = 14112; TBS[37][18] = 14688; TBS[38][18] = 15264; TBS[39][18] = 15264; TBS[40][18] = 15840;
	TBS[41][18] = 16416; TBS[42][18] = 16416; TBS[43][18] = 16992; TBS[44][18] = 17568; TBS[45][18] = 17568;
	TBS[46][18] = 18336; TBS[47][18] = 18336; TBS[48][18] = 19080; TBS[49][18] = 19080; TBS[50][18] = 19848;

	TBS[1][19] = 408; TBS[2][19] = 840; TBS[3][19] = 1288; TBS[4][19] = 1736; TBS[5][19] = 2152;
	TBS[6][19] = 2600; TBS[7][19] = 2984; TBS[8][19] = 3496; TBS[9][19] = 3880; TBS[10][19] = 4264;
	TBS[11][19] = 4776; TBS[12][19] = 5160; TBS[13][19] = 5544; TBS[14][19] = 5992; TBS[15][19] = 6456;
	TBS[16][19] = 6968; TBS[17][19] = 7224; TBS[18][19] = 7736; TBS[19][19] = 8248; TBS[20][19] = 8504;
	TBS[21][19] = 9144; TBS[22][19] = 9528; TBS[23][19] = 9912; TBS[24][19] = 10296; TBS[25][19] = 10680;
	TBS[26][19] = 11064; TBS[27][19] = 11448; TBS[28][19] = 12216; TBS[29][19] = 12576; TBS[30][19] = 12960;
	TBS[31][19] = 13536; TBS[32][19] = 13536; TBS[33][19] = 14112; TBS[34][19] = 14688; TBS[35][19] = 15264;
	TBS[36][19] = 15264; TBS[37][19] = 15840; TBS[38][19] = 16416; TBS[39][19] = 16992; TBS[40][19] = 16992;
	TBS[41][19] = 17568; TBS[42][19] = 18336; TBS[43][19] = 18336; TBS[44][19] = 19080; TBS[45][19] = 19080;
	TBS[46][19] = 19848; TBS[47][19] = 20616; TBS[48][19] = 20616; TBS[49][19] = 21384; TBS[50][19] = 21384;

	TBS[1][20] = 440; TBS[2][20] = 904; TBS[3][20] = 1384; TBS[4][20] = 1864; TBS[5][20] = 2344;
	TBS[6][20] = 2792; TBS[7][20] = 3240; TBS[8][20] = 3752; TBS[9][20] = 4136; TBS[10][20] = 4584;
	TBS[11][20] = 5160; TBS[12][20] = 5544; TBS[13][20] = 5992; TBS[14][20] = 6456; TBS[15][20] = 6968;
	TBS[16][20] = 7480; TBS[17][20] = 7992; TBS[18][20] = 8248; TBS[19][20] = 8760; TBS[20][20] = 9144;
	TBS[21][20] = 9912; TBS[22][20] = 10296; TBS[23][20] = 10680; TBS[24][20] = 11064; TBS[25][20] = 11448;
	TBS[26][20] = 12216; TBS[27][20] = 12576; TBS[28][20] = 12960; TBS[29][20] = 13536; TBS[30][20] = 14112;
	TBS[31][20] = 14688; TBS[32][20] = 14688; TBS[33][20] = 15264; TBS[34][20] = 15840; TBS[35][20] = 16416;
	TBS[36][20] = 16992; TBS[37][20] = 16992; TBS[38][20] = 17568; TBS[39][20] = 18336; TBS[40][20] = 18336;
	TBS[41][20] = 19080; TBS[42][20] = 19848; TBS[43][20] = 19848; TBS[44][20] = 20616; TBS[45][20] = 20616;
	TBS[46][20] = 21384; TBS[47][20] = 22152; TBS[48][20] = 22152; TBS[49][20] = 22920; TBS[50][20] = 22920;

	TBS[1][21] = 488; TBS[2][21] = 1000; TBS[3][21] = 1480; TBS[4][21] = 1992; TBS[5][21] = 2472;
	TBS[6][21] = 2984; TBS[7][21] = 3496; TBS[8][21] = 4008; TBS[9][21] = 4584; TBS[10][21] = 4968;
	TBS[11][21] = 5544; TBS[12][21] = 5992; TBS[13][21] = 6456; TBS[14][21] = 6968; TBS[15][21] = 7480;
	TBS[16][21] = 7992; TBS[17][21] = 8504; TBS[18][21] = 9144; TBS[19][21] = 9528; TBS[20][21] = 9912;
	TBS[21][21] = 10680; TBS[22][21] = 11064; TBS[23][21] = 11448; TBS[24][21] = 12216; TBS[25][21] = 12576;
	TBS[26][21] = 12960; TBS[27][21] = 13536; TBS[28][21] = 14112; TBS[29][21] = 14688; TBS[30][21] = 15264;
	TBS[31][21] = 15840; TBS[32][21] = 15840; TBS[33][21] = 16416; TBS[34][21] = 16992; TBS[35][21] = 17568;
	TBS[36][21] = 18336; TBS[37][21] = 18336; TBS[38][21] = 19080; TBS[39][21] = 19848; TBS[40][21] = 19848;
	TBS[41][21] = 20616; TBS[42][21] = 21384; TBS[43][21] = 21384; TBS[44][21] = 22152; TBS[45][21] = 22920;
	TBS[46][21] = 22920; TBS[47][21] = 23688; TBS[48][21] = 24496; TBS[49][21] = 24496; TBS[50][21] = 25456;

	TBS[1][22] = 520; TBS[2][22] = 1064; TBS[3][22] = 1608; TBS[4][22] = 2152; TBS[5][22] = 2664;
	TBS[6][22] = 3240; TBS[7][22] = 3752; TBS[8][22] = 4264; TBS[9][22] = 4776; TBS[10][22] = 5352;
	TBS[11][22] = 5992; TBS[12][22] = 6456; TBS[13][22] = 6968; TBS[14][22] = 7480; TBS[15][22] = 7992;
	TBS[16][22] = 8504; TBS[17][22] = 9144; TBS[18][22] = 9528; TBS[19][22] = 10296; TBS[20][22] = 10680;
	TBS[21][22] = 11448; TBS[22][22] = 11832; TBS[23][22] = 12576; TBS[24][22] = 12960; TBS[25][22] = 13536;
	TBS[26][22] = 14112; TBS[27][22] = 14688; TBS[28][22] = 15264; TBS[29][22] = 15840; TBS[30][22] = 16416;
	TBS[31][22] = 16992; TBS[32][22] = 16992; TBS[33][22] = 17568; TBS[34][22] = 18336; TBS[35][22] = 19080;
	TBS[36][22] = 19080; TBS[37][22] = 19848; TBS[38][22] = 20616; TBS[39][22] = 21384; TBS[40][22] = 21384;
	TBS[41][22] = 22152; TBS[42][22] = 22920; TBS[43][22] = 22920; TBS[44][22] = 23688; TBS[45][22] = 24496;
	TBS[46][22] = 24496; TBS[47][22] = 25456; TBS[48][22] = 25456; TBS[49][22] = 26416; TBS[50][22] = 27376;

	TBS[1][23] = 552; TBS[2][23] = 1128; TBS[3][23] = 1736; TBS[4][23] = 2280; TBS[5][23] = 2856;
	TBS[6][23] = 3496; TBS[7][23] = 4008; TBS[8][23] = 4584; TBS[9][23] = 5160; TBS[10][23] = 5736;
	TBS[11][23] = 6200; TBS[12][23] = 6968; TBS[13][23] = 7480; TBS[14][23] = 7992; TBS[15][23] = 8504;
	TBS[16][23] = 9144; TBS[17][23] = 9912; TBS[18][23] = 10296; TBS[19][23] = 11064; TBS[20][23] = 11448;
	TBS[21][23] = 12216; TBS[22][23] = 12576; TBS[23][23] = 12960; TBS[24][23] = 13536; TBS[25][23] = 14112;
	TBS[26][23] = 14688; TBS[27][23] = 15264; TBS[28][23] = 15840; TBS[29][23] = 16416; TBS[30][23] = 16992;
	TBS[31][23] = 17568; TBS[32][23] = 18336; TBS[33][23] = 19080; TBS[34][23] = 19848; TBS[35][23] = 19848;
	TBS[36][23] = 20616; TBS[37][23] = 21384; TBS[38][23] = 22152; TBS[39][23] = 22152; TBS[40][23] = 22920;
	TBS[41][23] = 23688; TBS[42][23] = 24496; TBS[43][23] = 24496; TBS[44][23] = 25456; TBS[45][23] = 25456;
	TBS[46][23] = 26416; TBS[47][23] = 27376; TBS[48][23] = 27376; TBS[49][23] = 28336; TBS[50][23] = 28336;

	TBS[1][24] = 584; TBS[2][24] = 1192; TBS[3][24] = 1800; TBS[4][24] = 2408; TBS[5][24] = 2984;
	TBS[6][24] = 3624; TBS[7][24] = 4264; TBS[8][24] = 4968; TBS[9][24] = 5544; TBS[10][24] = 5992;
	TBS[11][24] = 6712; TBS[12][24] = 7224; TBS[13][24] = 7992; TBS[14][24] = 8504; TBS[15][24] = 9144;
	TBS[16][24] = 9912; TBS[17][24] = 10296; TBS[18][24] = 11064; TBS[19][24] = 11448; TBS[20][24] = 12216;
	TBS[21][24] = 12960; TBS[22][24] = 13536; TBS[23][24] = 14112; TBS[24][24] = 14688; TBS[25][24] = 15264;
	TBS[26][24] = 15840; TBS[27][24] = 16416; TBS[28][24] = 16992; TBS[29][24] = 17568; TBS[30][24] = 18336;
	TBS[31][24] = 19080; TBS[32][24] = 19848; TBS[33][24] = 19848; TBS[34][24] = 20616; TBS[35][24] = 21384;
	TBS[36][24] = 22152; TBS[37][24] = 22920; TBS[38][24] = 22920; TBS[39][24] = 23688; TBS[40][24] = 24496;
	TBS[41][24] = 25456; TBS[42][24] = 25456; TBS[43][24] = 26416; TBS[44][24] = 26416; TBS[45][24] = 27376;
	TBS[46][24] = 28336; TBS[47][24] = 28336; TBS[48][24] = 29296; TBS[49][24] = 29296; TBS[50][24] = 30576;

	TBS[1][25] = 616; TBS[2][25] = 1256; TBS[3][25] = 1864; TBS[4][25] = 2536; TBS[5][25] = 3112;
	TBS[6][25] = 3752; TBS[7][25] = 4392; TBS[8][25] = 5160; TBS[9][25] = 5736; TBS[10][25] = 6200;
	TBS[11][25] = 6968; TBS[12][25] = 7480; TBS[13][25] = 8248; TBS[14][25] = 8760; TBS[15][25] = 9528;
	TBS[16][25] = 10296; TBS[17][25] = 10680; TBS[18][25] = 11448; TBS[19][25] = 12216; TBS[20][25] = 12576;
	TBS[21][25] = 13536; TBS[22][25] = 14112; TBS[23][25] = 14688; TBS[24][25] = 15264; TBS[25][25] = 15840;
	TBS[26][25] = 16416; TBS[27][25] = 16992; TBS[28][25] = 17568; TBS[29][25] = 18336; TBS[30][25] = 19080;
	TBS[31][25] = 19848; TBS[32][25] = 20616; TBS[33][25] = 20616; TBS[34][25] = 21384; TBS[35][25] = 22152;
	TBS[36][25] = 22920; TBS[37][25] = 23688; TBS[38][25] = 24496; TBS[39][25] = 24496; TBS[40][25] = 25456;
	TBS[41][25] = 26416; TBS[42][25] = 26416; TBS[43][25] = 27376; TBS[44][25] = 28336; TBS[45][25] = 28336;
	TBS[46][25] = 29296; TBS[47][25] = 29296; TBS[48][25] = 30576; TBS[49][25] = 31704; TBS[50][25] = 31704;

	TBS[1][26] = 712;   TBS[2][26] = 1480; TBS[3][26] = 2216; TBS[4][26] = 2984; TBS[5][26] = 3752;
	TBS[6][26] = 4392;  TBS[7][26] = 5160; TBS[8][26] = 5992; TBS[9][26] = 6712; TBS[10][26] = 7480;
	TBS[11][26] = 8248;  TBS[12][26] = 8760; TBS[13][26] = 9528; TBS[14][26] = 10296; TBS[15][26] = 11064;
	TBS[16][26] = 11832; TBS[17][26] = 12576; TBS[18][26] = 13536; TBS[19][26] = 14112; TBS[20][26] = 14688;
	TBS[21][26] = 15264; TBS[22][26] = 16416; TBS[23][26] = 16992; TBS[24][26] = 17568; TBS[25][26] = 18336;
	TBS[26][26] = 19080; TBS[27][26] = 19848; TBS[28][26] = 20616; TBS[29][26] = 21384; TBS[30][26] = 22152;
	TBS[31][26] = 22920; TBS[32][26] = 23688; TBS[33][26] = 24496; TBS[34][26] = 25456; TBS[35][26] = 25456;
	TBS[36][26] = 26416; TBS[37][26] = 27376; TBS[38][26] = 28336; TBS[39][26] = 29296; TBS[40][26] = 29296;
	TBS[41][26] = 30576; TBS[42][26] = 30576; TBS[43][26] = 31704; TBS[44][26] = 32856; TBS[45][26] = 32856;
	TBS[46][26] = 34008; TBS[47][26] = 35160; TBS[48][26] = 35160; TBS[49][26] = 36696; TBS[50][26] = 36696;


	/*
	for( int i = 0 ; i < 51 ; i++ )
	{
	for( int j = 0 ; j < 27 ; j++)
	{
	TBS[i][j] = TBS[i][j] * 2 ;
	}
	}
	*/
}// end of TBS table setting


void table_setting_RBIR(void)
{
	RBIR_value[0][0] = 0.0072;   RBIR_value[0][1] = 0.008;     RBIR_value[0][2] = 0.009;
	RBIR_value[0][3] = 0.0101;   RBIR_value[0][4] = 0.0114;    RBIR_value[0][5] = 0.0127;
	RBIR_value[0][6] = 0.0143;   RBIR_value[0][7] = 0.0159;    RBIR_value[0][8] = 0.0179;
	RBIR_value[0][9] = 0.02;     RBIR_value[0][10] = 0.0225;   RBIR_value[0][11] = 0.0251;
	RBIR_value[0][12] = 0.0282;   RBIR_value[0][13] = 0.0315;   RBIR_value[0][14] = 0.0352;
	RBIR_value[0][15] = 0.0394;   RBIR_value[0][16] = 0.0442;   RBIR_value[0][17] = 0.0493;
	RBIR_value[0][18] = 0.0551;   RBIR_value[0][19] = 0.0616;   RBIR_value[0][20] = 0.0688;
	RBIR_value[0][21] = 0.0767;   RBIR_value[0][22] = 0.0855;   RBIR_value[0][23] = 0.0953;
	RBIR_value[0][24] = 0.1061;   RBIR_value[0][25] = 0.118;    RBIR_value[0][26] = 0.1311;
	RBIR_value[0][27] = 0.1456;   RBIR_value[0][28] = 0.1615;   RBIR_value[0][29] = 0.1788;
	RBIR_value[0][30] = 0.1978;   RBIR_value[0][31] = 0.2184;   RBIR_value[0][32] = 0.2407;
	RBIR_value[0][33] = 0.265;    RBIR_value[0][34] = 0.291;    RBIR_value[0][35] = 0.319;
	RBIR_value[0][36] = 0.3489;   RBIR_value[0][37] = 0.3806;   RBIR_value[0][38] = 0.4141;
	RBIR_value[0][39] = 0.4493;   RBIR_value[0][40] = 0.4859;   RBIR_value[0][41] = 0.5239;
	RBIR_value[0][42] = 0.5628;   RBIR_value[0][43] = 0.6024;   RBIR_value[0][44] = 0.6422;
	RBIR_value[0][45] = 0.6817;   RBIR_value[0][46] = 0.7207;   RBIR_value[0][47] = 0.7584;
	RBIR_value[0][48] = 0.7944;   RBIR_value[0][49] = 0.8281;   RBIR_value[0][50] = 0.8592;
	RBIR_value[0][51] = 0.8872;   RBIR_value[0][52] = 0.9119;   RBIR_value[0][53] = 0.9331;
	RBIR_value[0][54] = 0.9507;   RBIR_value[0][55] = 0.9649;   RBIR_value[0][56] = 0.976;
	RBIR_value[0][57] = 0.9842;   RBIR_value[0][58] = 0.9901;   RBIR_value[0][59] = 0.9942;
	RBIR_value[0][60] = 0.9968;   RBIR_value[0][61] = 0.9983;   RBIR_value[0][62] = 0.9992;
	RBIR_value[0][63] = 0.9997;   RBIR_value[0][64] = 0.9999;   RBIR_value[0][65] = 1;
	RBIR_value[0][66] = 1;        RBIR_value[0][67] = 1;        RBIR_value[0][68] = 1;
	RBIR_value[0][69] = 1;        RBIR_value[0][70] = 1;        RBIR_value[0][71] = 1;
	RBIR_value[0][72] = 1;        RBIR_value[0][73] = 1;        RBIR_value[0][74] = 1;
	RBIR_value[0][75] = 1;        RBIR_value[0][76] = 1;        RBIR_value[0][77] = 1;
	RBIR_value[0][78] = 1;        RBIR_value[0][79] = 1;        RBIR_value[0][80] = 1;
	RBIR_value[0][81] = 1;        RBIR_value[0][82] = 1;        RBIR_value[0][83] = 1;
	RBIR_value[0][84] = 1;        RBIR_value[0][85] = 1;        RBIR_value[0][86] = 1;
	RBIR_value[0][87] = 1;        RBIR_value[0][88] = 1;        RBIR_value[0][89] = 1;
	RBIR_value[0][90] = 1;        RBIR_value[0][91] = 1;        RBIR_value[0][92] = 1;
	RBIR_value[0][93] = 1;        RBIR_value[0][94] = 1;

	RBIR_value[1][0] = 0.0036;   RBIR_value[1][1] = 0.004;    RBIR_value[1][2] = 0.0045;
	RBIR_value[1][3] = 0.005;    RBIR_value[1][4] = 0.0057;   RBIR_value[1][5] = 0.0063;
	RBIR_value[1][6] = 0.0071;   RBIR_value[1][7] = 0.008;    RBIR_value[1][8] = 0.0089;
	RBIR_value[1][9] = 0.01;     RBIR_value[1][10] = 0.0112;   RBIR_value[1][11] = 0.0126;
	RBIR_value[1][12] = 0.0141;   RBIR_value[1][13] = 0.0158;   RBIR_value[1][14] = 0.0176;
	RBIR_value[1][15] = 0.0197;   RBIR_value[1][16] = 0.0221;   RBIR_value[1][17] = 0.0247;
	RBIR_value[1][18] = 0.0276;   RBIR_value[1][19] = 0.0308;   RBIR_value[1][20] = 0.0344;
	RBIR_value[1][21] = 0.0384;   RBIR_value[1][22] = 0.0428;   RBIR_value[1][23] = 0.0476;
	RBIR_value[1][24] = 0.0531;   RBIR_value[1][25] = 0.059;    RBIR_value[1][26] = 0.0656;
	RBIR_value[1][27] = 0.0728;   RBIR_value[1][28] = 0.0808;   RBIR_value[1][29] = 0.0895;
	RBIR_value[1][30] = 0.099;    RBIR_value[1][31] = 0.1094;   RBIR_value[1][32] = 0.1206;
	RBIR_value[1][33] = 0.1329;   RBIR_value[1][34] = 0.1461;   RBIR_value[1][35] = 0.1603;
	RBIR_value[1][36] = 0.1756;   RBIR_value[1][37] = 0.192;    RBIR_value[1][38] = 0.2094;
	RBIR_value[1][39] = 0.2279;   RBIR_value[1][40] = 0.2474;   RBIR_value[1][41] = 0.268;
	RBIR_value[1][42] = 0.2896;   RBIR_value[1][43] = 0.3122;   RBIR_value[1][44] = 0.3357;
	RBIR_value[1][45] = 0.36;     RBIR_value[1][46] = 0.3852;   RBIR_value[1][47] = 0.4112;
	RBIR_value[1][48] = 0.4379;   RBIR_value[1][49] = 0.4653;   RBIR_value[1][50] = 0.4933;
	RBIR_value[1][51] = 0.5219;   RBIR_value[1][52] = 0.5509;   RBIR_value[1][53] = 0.5804;
	RBIR_value[1][54] = 0.6103;   RBIR_value[1][55] = 0.6403;   RBIR_value[1][56] = 0.6709;
	RBIR_value[1][57] = 0.7014;   RBIR_value[1][58] = 0.7317;   RBIR_value[1][59] = 0.7617;
	RBIR_value[1][60] = 0.791;    RBIR_value[1][61] = 0.8193;   RBIR_value[1][62] = 0.8463;
	RBIR_value[1][63] = 0.8716;   RBIR_value[1][64] = 0.8949;   RBIR_value[1][65] = 0.9158;
	RBIR_value[1][66] = 0.9343;   RBIR_value[1][67] = 0.9501;   RBIR_value[1][68] = 0.9633;
	RBIR_value[1][69] = 0.9739;   RBIR_value[1][70] = 0.9821;   RBIR_value[1][71] = 0.9883;
	RBIR_value[1][72] = 0.9927;   RBIR_value[1][73] = 0.9957;   RBIR_value[1][74] = 0.9976;
	RBIR_value[1][75] = 0.9988;   RBIR_value[1][76] = 0.9994;   RBIR_value[1][77] = 0.9997;
	RBIR_value[1][78] = 0.9999;   RBIR_value[1][79] = 1;        RBIR_value[1][80] = 1;
	RBIR_value[1][81] = 1;        RBIR_value[1][82] = 1;        RBIR_value[1][83] = 1;
	RBIR_value[1][84] = 1;        RBIR_value[1][85] = 1;        RBIR_value[1][86] = 1;
	RBIR_value[1][87] = 1;        RBIR_value[1][88] = 1;        RBIR_value[1][89] = 1;
	RBIR_value[1][90] = 1;        RBIR_value[1][91] = 1;        RBIR_value[1][92] = 1;
	RBIR_value[1][93] = 1;        RBIR_value[1][94] = 1;


	RBIR_value[2][0] = 0.0024;   RBIR_value[2][1] = 0.0027;    RBIR_value[2][2] = 0.003;
	RBIR_value[2][3] = 0.0034;   RBIR_value[2][4] = 0.0038;    RBIR_value[2][5] = 0.0043;
	RBIR_value[2][6] = 0.0047;   RBIR_value[2][7] = 0.0054;    RBIR_value[2][8] = 0.006;
	RBIR_value[2][9] = 0.0067;   RBIR_value[2][10] = 0.0075;   RBIR_value[2][11] = 0.0084;
	RBIR_value[2][12] = 0.0094;   RBIR_value[2][13] = 0.0106;   RBIR_value[2][14] = 0.0117;
	RBIR_value[2][15] = 0.0132;   RBIR_value[2][16] = 0.0147;   RBIR_value[2][17] = 0.0165;
	RBIR_value[2][18] = 0.0184;   RBIR_value[2][19] = 0.0207;   RBIR_value[2][20] = 0.0229;
	RBIR_value[2][21] = 0.0257;   RBIR_value[2][22] = 0.0285;   RBIR_value[2][23] = 0.0319;
	RBIR_value[2][24] = 0.0354;   RBIR_value[2][25] = 0.0396;   RBIR_value[2][26] = 0.0437;
	RBIR_value[2][27] = 0.0488;   RBIR_value[2][28] = 0.0539;   RBIR_value[2][29] = 0.0599;
	RBIR_value[2][30] = 0.066;    RBIR_value[2][31] = 0.0732;   RBIR_value[2][32] = 0.0805;
	RBIR_value[2][33] = 0.089;    RBIR_value[2][34] = 0.0974;   RBIR_value[2][35] = 0.1073;
	RBIR_value[2][36] = 0.1172;   RBIR_value[2][37] = 0.1285;   RBIR_value[2][38] = 0.1398;
	RBIR_value[2][39] = 0.1525;   RBIR_value[2][40] = 0.1653;   RBIR_value[2][41] = 0.1795;
	RBIR_value[2][42] = 0.1937;   RBIR_value[2][43] = 0.2092;   RBIR_value[2][44] = 0.2247;
	RBIR_value[2][45] = 0.2415;   RBIR_value[2][46] = 0.2583;   RBIR_value[2][47] = 0.2763;
	RBIR_value[2][48] = 0.2942;   RBIR_value[2][49] = 0.3132;   RBIR_value[2][50] = 0.3321;
	RBIR_value[2][51] = 0.3519;   RBIR_value[2][52] = 0.3718;   RBIR_value[2][53] = 0.3924;
	RBIR_value[2][54] = 0.4131;   RBIR_value[2][55] = 0.4345;   RBIR_value[2][56] = 0.4558;
	RBIR_value[2][57] = 0.4778;   RBIR_value[2][58] = 0.4997;   RBIR_value[2][59] = 0.5223;
	RBIR_value[2][60] = 0.5448;   RBIR_value[2][61] = 0.5677;   RBIR_value[2][62] = 0.5907;
	RBIR_value[2][63] = 0.6141;   RBIR_value[2][64] = 0.6374;   RBIR_value[2][65] = 0.6611;
	RBIR_value[2][66] = 0.6848;   RBIR_value[2][67] = 0.7087;   RBIR_value[2][68] = 0.7325;
	RBIR_value[2][69] = 0.7564;   RBIR_value[2][70] = 0.7802;   RBIR_value[2][71] = 0.8036;
	RBIR_value[2][72] = 0.8269;   RBIR_value[2][73] = 0.8489;   RBIR_value[2][74] = 0.8708;
	RBIR_value[2][75] = 0.8904;   RBIR_value[2][76] = 0.91;     RBIR_value[2][77] = 0.9262;
	RBIR_value[2][78] = 0.9425;   RBIR_value[2][79] = 0.9547;   RBIR_value[2][80] = 0.9668;
	RBIR_value[2][81] = 0.9732;   RBIR_value[2][82] = 0.9796;   RBIR_value[2][83] = 0.984;
	RBIR_value[2][84] = 0.9883;   RBIR_value[2][85] = 0.991;    RBIR_value[2][86] = 0.9937;
	RBIR_value[2][87] = 0.9954;   RBIR_value[2][88] = 0.9971;   RBIR_value[2][89] = 0.9983;
	RBIR_value[2][90] = 0.9995;   RBIR_value[2][91] = 0.9998;   RBIR_value[2][92] = 1;
	RBIR_value[2][93] = 1;        RBIR_value[2][94] = 1;

}



#if 0  // unused: Load_BLER_vienna_file — never called (commented out in StandardInitialization)
void Load_BLER_vienna_file()
{
	//string trash;

	int SINR_dB_idx = 0;
	int CQI_idx = 0;
	Real BLER_linear = 0.;
	int loop = 0;


	ifstream BLER_vienna_input;
	BLER_vienna_input.open("BLER_vienna.dat");

	if (BLER_vienna_input.is_open())
	{
		cout << "=== BLER_vienna.DAT file open ===" << endl;
	}
	else
	{
		cout << "=== NO FILE!!!!!! Please CHECK ===" << endl;
		getchar();
	}

	while (!BLER_vienna_input.eof())
	{
		for (int cqi_index = 0; cqi_index < 15; cqi_index++)
		{
			BLER_vienna_input >> BLER_linear;

			/// vienna BLER [ SINR_dB ][ CQI_idx ]
			BLER_vienna[loop][cqi_index] = BLER_linear;



		}
		//getchar();
		//cout << "loop = " << loop << endl;

		loop++;
	}

	cout << "=== BLER_vienna.DAT load Complete ===" << endl << endl;
	BLER_vienna_input.close();

}
#endif  // unused: Load_BLER_vienna_file



void Load_5G_BLER_vienna_file()
{
	//string trash;

	int SINR_dB_idx = 0;
	int CQI_idx = 0;
	Real BLER_linear = 0.;
	int loop = 0;


	ifstream NEW5GBLER_input;
	NEW5GBLER_input.open("BLER_5G.dat");

	if (NEW5GBLER_input.is_open())
	{
		cout << "=== NEW5GBLER.dat file open ===" << endl;
	}
	else
	{
		cout << "=== NO NEW5GBLER.dat FILE!!!!!! Please CHECK ===" << endl;
		getchar();
	}

	while (!NEW5GBLER_input.eof())
	{
		for (int cqi_index = 0; cqi_index < 15; cqi_index++)
		{
			NEW5GBLER_input >> BLER_linear;

			/// 5GNEWBLER.dat [ SINR_dB ][ CQI_idx ]
			NEW5GBLER[loop][cqi_index] = BLER_linear;

		}


		//getchar();
		//cout << "loop = " << loop << endl;

		loop++;
	}

	cout << "=== NEW5GBLER.dat load Complete ===" << endl;
	NEW5GBLER_input.close();

}