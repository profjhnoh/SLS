#include <iostream>
#include <fstream>

using namespace std;

#include <cstring>
#include <cstdlib>

#include "common.h"
#include "const.h"
/*===================================================================
FUNCTION: PerDropStatistics( void )

AUTHOR: SJ

DESCRIPTION: Per-Drop Statistics are output

NOTES:

===================================================================*/
double Total_Avr;
double Total_Min;
double SF_avr;

/*===================================================================
FUNCTION: int getIndex(vector<int> v, int K)
===================================================================*/
int getIndex(vector<int> v, int K);

void PerDropStatistics(int drop_idx)
{

	//// For calibration test
	/*
	ofstream TESToutputfile("TEST_output_file_MS.dat");
	TESToutputfile << "--------------------------- Simulation Configure ---------------------------" << endl;
	TESToutputfile << " MS_index " << '\t' << "MS_StaticGain_dB" << '\t' << "MS_interference_dB" << '\t' << "MS_Geometry_dB" << '\t' << "MS_throughput" << '\t' << "MS_sum_of_ch_coef_mean" << '\t' << "MS_BLER_mean" << '\t' << "MS_decision_BLER_mean" << endl;
	TESToutputfile << "----------------------------------------------------------------------------" << endl;
	TESToutputfile << "--------------------------- drop_idx = " << drop_idx + 1 << " --------------------------------- " << endl;

	for (int ms_idx = 0; ms_idx < num_MS; ms_idx++)
	{
		double linear_signal = dBm2linear(link[ms_idx].str_signal);
		double linear_interference = dBm2linear(link[ms_idx].interference);

		linear_interference = linear_interference + noise;

		double sum_of_ch_coef = 0.;

		for (int rb_idx = 0; rb_idx < num_rb; rb_idx++)
		{
			
			for (int r_idx = 0; r_idx < NUM_RX; r_idx++) // Received    Channel
			{
				for (int t_idx = 0; t_idx < NUM_TX; t_idx++) // Transmitted Channel
				{
					sum_of_ch_coef += pow(abs(ms[ms_idx].H_m[rb_idx](r_idx, t_idx)), 2.);
				}
			}
		}

		TESToutputfile << ms_idx << '\t' << link[ms_idx].str_signal << '\t' << link[ms_idx].interference << '\t' << linear2dB(linear_signal / linear_interference) << '\t' << ms[ms_idx].Return_Throughput(ms_idx) * 1000. / run_times * pow(10., -7.) << '\t' << sum_of_ch_coef / num_rb << '\t' << ms[ms_idx].BLER_SUM_Value / (double)ms[ms_idx].BLER_SUM_count << '\t' << ms[ms_idx].decision_BLER_SUM_Value / (double)ms[ms_idx].decision_BLER_SUM_count << endl;
	}
	TESToutputfile.close();
	*/


	if (drop_idx == 0)
	{
		Total_Avr = 0.;
		Total_Min = 0.;
		SF_avr = 0.;
	}


	// For Calibration
	/*
	ofstream TESToutputSECTOR("TEST_output_file_SECTOR.dat");
	TESToutputSECTOR << "--------------------------- Simulation Configure ---------------------------" << endl;
	TESToutputSECTOR << " Sector_index " << '\t' << "link_MS_Size" << '\t' << "link_MS_index" << '\t' << "Sector_throughput" << endl;
	TESToutputSECTOR << "----------------------------------------------------------------------------" << endl;
	TESToutputSECTOR << "--------------------------- drop_idx = " << drop_idx + 1 << " --------------------------------- " << endl;
	*/

	// Determine sector range based on single_cell_mode
	// In single_cell_mode, only analyze center BS (BS 0) = sectors 0, 1, 2
	int start_sector = 0;
	int end_sector = (single_cell_mode == 1) ? 3 : num_SECTORS;
	int num_active_sectors = end_sector - start_sector;

	double *per_sector_thru = new double[num_SECTORS];
	//double per_sector_thru[num_SECTORS];

	for (int sector_idx = start_sector; sector_idx < end_sector; sector_idx++)
	{
		per_sector_thru[sector_idx] = 0;
		for (int ue_idx = 0; ue_idx < (int)sector[sector_idx].ue_in_control.size(); ue_idx++)
		{
			int ue2consider = sector[sector_idx].ue_in_control[ue_idx];
			per_sector_thru[sector_idx] += ms[ue2consider].Return_Throughput(ue2consider);
			//TESToutputSECTOR << sector_idx << '\t' << (int)sector[sector_idx].ue_in_control.size() << '\t' << ue2consider << '\t' << per_sector_thru[sector_idx] * 1000. / run_times * pow(10., -7.) << endl;
		}
	}

	//TESToutputSECTOR.close();
	double avg_thru = 0;
	double min_thru = per_sector_thru[start_sector];
	double max_thru = per_sector_thru[start_sector];


	for (int sector_idx = start_sector; sector_idx < end_sector; sector_idx++)
	{
		avg_thru += per_sector_thru[sector_idx];
		if (min_thru > per_sector_thru[sector_idx])
		{
			min_thru = per_sector_thru[sector_idx];
		}

		if (max_thru < per_sector_thru[sector_idx])
		{
			max_thru = per_sector_thru[sector_idx];
		}
	}
	
	
	ofstream  OutputFile;
	if (drop_idx == 0)
	{
		/// time output
		//time_t t = time(NULL);
		//struct tm tm = *localtime(&t);

		char c_year[10];
		char c_mon [10];
		char c_day [10];
		char c_hour[10];
		char c_min [10];


		//int year = tm.tm_year + 1900;
		//int mon = tm.tm_mon + 1;
		//int day = tm.tm_mday;
		//int hour = tm.tm_hour;
		//int min = tm.tm_min;

		sprintf(c_year, "%04d", g_year);
		sprintf(c_mon, "%02d",  g_mon);
		sprintf(c_day, "%02d",  g_day);
		sprintf(c_hour, "%02d", g_hour);
		sprintf(c_min, "%02d",  g_min);

		OutputFile.open(file_name, ios::out);
		
		OutputFile << "configuration_file : " << cfg_name << endl;
		OutputFile << "time : " << c_year << "-"<< c_mon << "-" << c_day << "_" << c_hour << "-" << c_min << endl;
		OutputFile << endl;

		OutputFile << "Configuration :" << endl;
		OutputFile << "    " << "_seed                   : " << _seed << endl;
		OutputFile << "    " << "num_drops               : " << num_drop << endl;
		OutputFile << "    " << "run_times               : " << run_times << endl;
		OutputFile << "    " << "num_user_per_sector     : " << num_MS_persector << endl;
		OutputFile << "    " << "Indoor_TRxP             : " << num_Indoor_TRxP << endl;
		OutputFile << "    " << "simple_num_BS           : " << num_BS << endl;
		OutputFile << "    " << "scenario                : " << TYPE << endl;
		OutputFile << "    " << "Scheduling_Type         : " << Scheduling_Type << endl;
		OutputFile << "    " << "Num_RBs                 : " << num_rb                  << endl;
		OutputFile << "    " << "fft_size                : " << fft_size                << endl;
		OutputFile << "    " << "subcarrier_spacing      : " << subcarrier_spacing      << endl;
		OutputFile << "    " << "bandwidth               : " << bandwidth               << endl;
		OutputFile << "    " << "numerology              : " << numerology              << endl;
		OutputFile << "    " << "carrier_freq            : " << carrier_freq            << endl;
		OutputFile << "    " << "Calibration_mode        : " << Calibration_mode        << endl;
		OutputFile << "    " << "BS_M                    : " << BS_M                    << endl;
		OutputFile << "    " << "BS_N                    : " << BS_N                    << endl;
		OutputFile << "    " << "BS_P                    : " << BS_P                    << endl;
		OutputFile << "    " << "BS_Mg                   : " << BS_Mg                   << endl;
		OutputFile << "    " << "BS_Ng                   : " << BS_Ng                   << endl;
		OutputFile << "    " << "BS_Mp                   : " << BS_Mp                   << endl;
		OutputFile << "    " << "BS_Np                   : " << BS_Np                   << endl;
		OutputFile << "    " << "BS_dH                   : " << BS_dH                   << endl;
		OutputFile << "    " << "BS_dV                   : " << BS_dV                   << endl;
		OutputFile << "    " << "BS_dgH                  : " << BS_dgH                  << endl;
		OutputFile << "    " << "BS_dgV                  : " << BS_dgV                  << endl;
		OutputFile << "    " << "MS_M                    : " << MS_M                    << endl;
		OutputFile << "    " << "MS_N                    : " << MS_N                    << endl;
		OutputFile << "    " << "MS_P                    : " << MS_P                    << endl;
		OutputFile << "    " << "MS_Mg                   : " << MS_Mg                   << endl;
		OutputFile << "    " << "MS_Ng                   : " << MS_Ng                   << endl;
		OutputFile << "    " << "MS_Mp                   : " << MS_Mp                   << endl;
		OutputFile << "    " << "MS_Np                   : " << MS_Np                   << endl;
		OutputFile << "    " << "MS_dH                   : " << MS_dH                   << endl;
		OutputFile << "    " << "MS_dV                   : " << MS_dV                   << endl;
		OutputFile << "    " << "MS_dgH                  : " << MS_dgH                  << endl;
		OutputFile << "    " << "MS_dgV                  : " << MS_dgV                  << endl;
		OutputFile << "    " << "Mechanic_tilt           : " << Mechanic_tilt           << endl;
		OutputFile << "    " << "num_compute_coef        : " << num_compute_coef        << endl;
		OutputFile << "    " << "num_neighbor            : " << num_neighbor            << endl;
		OutputFile << "    " << "grid_interval           : " << grid_interval           << endl;
		OutputFile << "    " << "cqi_history_length      : " << cqi_history_length      << endl;
		OutputFile << "    " << "mx_ue_mumimo            : " << mx_ue_mumimo            << endl;
		OutputFile << "    " << "NUM_UE_Layer            : " << NUM_UE_Layer            << endl;
		OutputFile << "    " << "N_pf                    : " << N_pf                    << endl;
		OutputFile << "    " << "SCHEDULE_DELAY          : " << SCHEDULE_DELAY          << endl;
		OutputFile << "    " << "num_of_threads          : " << num_of_threads          << endl;
		OutputFile << "    " << "g_comp_mode             : " << g_comp_mode             << endl;
		OutputFile << "    " << "comp_ue_pct             : " << comp_ue_pct             << endl;
		OutputFile << "    " << "g_static_gain_ratio_comp: " << g_static_gain_ratio_comp<< endl;
		OutputFile << endl;

		OutputFile << "Results : " << drop_idx + 1 << endl;

	}
	else { OutputFile.open(file_name, ios::app); }

	// 3GPP standard: slot duration = 1ms / 2^μ
	double slot_duration = 1.0e-3 / pow(2.0, numerology);  // seconds

	OutputFile << "    " << endl;
	OutputFile << "  - " << "Sector_Thruput:" << endl;
    OutputFile << "    " << "    " << "Avg : " << avg_thru / num_active_sectors / (run_times*slot_duration) * 1e-6 << endl;
    OutputFile << "    " << "    " << "Max : " << max_thru                      / (run_times*slot_duration) * 1e-6 << endl;
    OutputFile << "    " << "    " << "Min : " << min_thru                      / (run_times*slot_duration) * 1e-6 << endl;
	OutputFile << "    " << endl;

	Total_Avr += avg_thru / num_active_sectors / (run_times*slot_duration) * 1e-6;

	//compute ue thru 
	unsigned long int  *per_ue_thru = new unsigned long int [num_MS];
	unsigned long int  *ue_statistics = new unsigned long int [num_MS]; 

	//long double per_ue_thru[num_MS];   //
	//long double ue_statistics[num_MS]; //
	int cell_edge_ue_idx = (int)(num_MS * 0.05); // cell edge user

	ofstream  Througput;

	char ue_throughput_file_name[100];
	sprintf(ue_throughput_file_name,"./%s/UE_througpht_drop_%d.csv",folder_name,drop_idx);
	Througput.open(ue_throughput_file_name, ios::out);

	/*
	if (drop_idx == 0)
		Througput.open("UE_througpht_drop_0.csv", ios::out);
	else if (drop_idx == 1)
		Througput.open("UE_througpht_drop_1.csv", ios::out);		
	else if (drop_idx == 2)
		Througput.open("UE_througpht_drop_2.csv", ios::out);		
	else if (drop_idx == 3)
		Througput.open("UE_througpht_drop_3.csv", ios::out);		
	else if (drop_idx == 4)
		Througput.open("UE_througpht_drop_4.csv", ios::out);	
	*/

	Througput
	<< "ue_idx" <<","                  //  1
	<< "per_ue_thru" <<","             //  2
	<< "avr_cqi[scheduling/(log(1+SNR)]" <<","                 //  3
	<< "total_num_scheduled" <<","     //  4
	<< "total_num_tx" <<","            //  5
	<< "failed" <<","                  //  6
	<< "total_num_re_tx" <<","         //  7
	<< "re_tx_failed" <<","            //  8
	<< "avr cqi idx[ReceiveDL]" <<","                 //  9
	<< "avr mcs" <<","                 // 10
	<< "avr effective sinr"<<","                // 11
	<< "avr sinr estimate[After scheduling]"<<","        // 12
	<< "num_comp_tx"<<","              // 13
	<< "bs_idx"<<","                   // 14
	<< "sector_idx"<<","               // 15
	<< "comp_sector_idx"<<","          // 16
	<< "ue_x"<<","                     // 17
	<< "ue_y"                          // 18
	<< endl;

	for (int ue_idx = 0; ue_idx < num_MS; ue_idx++)
	{
		per_ue_thru[ue_idx] = ms[ue_idx].Return_Throughput(ue_idx);
		int ue_idx_in_sector = getIndex( sector[links[ue_idx]._sector_in_control].ue_in_control , ue_idx );

		double avr_cqi = 0;
		for (int rb_idx = 0; rb_idx < num_rb; rb_idx++)
		{
			avr_cqi += sector[links[ue_idx]._sector_in_control].CQI_AVR[ue_idx_in_sector][rb_idx]/num_rb;
		}

		Througput
		<< ue_idx <<","                                                               //  1
		<< per_ue_thru[ue_idx]/(run_times*slot_duration) * 1e-6 <<","                 //  2
		<< avr_cqi <<","                                                              //  3
		<< links[ue_idx].total_num_scheduled <<","                                    //  4
		<< links[ue_idx].total_num_tx <<","                                           //  5
		<< links[ue_idx].failed <<","                                                 //  6
		<< links[ue_idx].total_num_re_tx <<","                                        //  7
		<< links[ue_idx].re_tx_failed <<","                                           //  8
		<< (double) ms[ue_idx].sum_cqi/ (double)links[ue_idx].total_num_tx <<","      //  9
		<< (double) ms[ue_idx].sum_mcs_type/ (double)links[ue_idx].total_num_tx <<"," // 10
		<< ms[ue_idx].total_esinr/links[ue_idx].total_num_tx <<","                    // 11
		<< ms[ue_idx].total_estimate_sinr/links[ue_idx].total_num_tx <<","            // 12
		<< ms[ue_idx].num_comp_tx <<","                                               // 13
		<< links[ue_idx].self_bs_idx <<","                                            // 14
		<< links[ue_idx]._sector_in_control <<","                                     // 15
		<< links[ue_idx].comp_sector_idx <<","                                        // 16
		<< ms[ue_idx].loc.x <<","                                                     // 17
		<< ms[ue_idx].loc.y                                                           // 18
		<< endl;
		ue_statistics[ue_idx] = per_ue_thru[ue_idx];
	}


	//getchar();

	long double tmp;
	bool swap;
	do
	{
		swap = false; //sorting frequence in order of cqi
		for (int ue_idx = 0; ue_idx < num_MS - 1; ue_idx++)
		{
			if (ue_statistics[ue_idx] > ue_statistics[ue_idx + 1])
			{
				tmp = ue_statistics[ue_idx];
				ue_statistics[ue_idx] = ue_statistics[ue_idx + 1];
				ue_statistics[ue_idx + 1] = tmp;
				swap = true;
			}
		}
	} while (swap);

	long double avg_ue_thru = 0;
	long double min_ue_thru = ue_statistics[cell_edge_ue_idx];
	long double max_ue_thru = ue_statistics[num_MS - 1];

	for (int ue_idx = 0; ue_idx < num_MS; ue_idx++)
	{
		avg_ue_thru += ue_statistics[ue_idx];
	}
	//getchar();

	double avg_initial_pe = 0;
	for (int ue_idx = 0; ue_idx < num_MS; ue_idx++)
	{
		avg_initial_pe += ms[ue_idx].Return_initial_packet_error() / num_MS;
	}

	/////////////////////// Average Capacity /////////////////
	double capacity_avr_drop = 0;
	for (int idx = 0; idx < (run_times - 1); idx++)
	{
		//capacity_avr_drop += Capacity_avr_per_drop[drop_idx][idx];
		capacity_avr_drop += 0.;
	}

	//////////////////////////////////////////////////////////
	OutputFile << "    " << "UE_Thruput:" << endl;
	OutputFile << "    " << "    " << "Avg : " << avg_ue_thru / num_MS / (run_times*slot_duration) * 1e-6 << endl;
	OutputFile << "    " << "    " << "Max : " << max_ue_thru          / (run_times*slot_duration) * 1e-6 << endl;
	OutputFile << "    " << "    " << "Min : " << min_ue_thru          / (run_times*slot_duration) * 1e-6 << endl;
	OutputFile << endl;

	//OutputFile << " Initial Packet Erro = " << avg_initial_pe << endl;	
	//OutputFile << " Average Spectral Efficiency (from MIMO capacity) = " << (capacity_avr_drop / (run_times - 1)) / bandwidth << "bits/s/Hz" << endl;
	//OutputFile << " Shannon Capacity (from MIMO capacity) = " << (capacity_avr_drop / (run_times - 1)) / (pow(10, 6)) << "Mbits/s" << endl;

	// 검증용 인석 추가
	OutputFile << "    " << "Schedule_map:" << endl;
	OutputFile << "    " << "    " << "sector_selected_ue : "        << sector_selected_ue                  / (run_times - SCHEDULE_DELAY) << endl;
	OutputFile << "    " << "    " << "scheduled_ue_mcs : "          << scheduled_ue_mcs                    / (run_times - SCHEDULE_DELAY) << endl;
	OutputFile << "    " << "    " << "scheduled_ue_cqi : "          << scheduled_ue_cqi                    / (run_times - SCHEDULE_DELAY) << endl;
	OutputFile << "    " << "    " << "scheduled_ue_widebandSINR : " << linear2dB(scheduled_ue_widebandSINR / (run_times - SCHEDULE_DELAY)) << endl;
	OutputFile  << endl;

	OutputFile << "    " << "Scheduling:" << endl;
	OutputFile << "    " << "    " << "METRIC : " << sector_metric     / (N_pf + run_times - SCHEDULE_DELAY) << endl;
	OutputFile << "    " << "    " << "CQI_read : " << sector_cqi_read / (N_pf + run_times - SCHEDULE_DELAY) << endl;
	OutputFile << "    " << "    " << "CQI_AVR : " << sector_cqi_avr   / (N_pf + run_times - SCHEDULE_DELAY) << endl;
	OutputFile  << endl;

	OutputFile << "    " << "Receive:" << endl;
	OutputFile << "    " << "    " << "ue_effective_sinr : " << linear2dB(ue_effective_sinr / (N_pf + run_times - SCHEDULE_DELAY)) << endl;
	OutputFile << "    " << "    " << "ue_info_bits : "      << ue_info_bits                / (N_pf + run_times - SCHEDULE_DELAY) << endl;
	OutputFile << "    " << "    " << "ue_bler : "           << ue_bler                     / (N_pf + run_times - SCHEDULE_DELAY) << endl;
	OutputFile  << endl;

	// 3GPP standard: use proper slot_duration instead of "1000 * numerology"
	// Original formula was incorrect approximation
	Total_Min += min_ue_thru / (run_times * slot_duration) * 1e-6;

	SF_avr += (capacity_avr_drop / (run_times - 1)) / bandwidth;

	/*
	if (drop_idx == num_drop - 1)
	{
		OutputFile << "Configuration:" << endl;
		OutputFile << "    " << "_seed : " << _seed << endl;
		OutputFile << "    " << "num_drops : " << num_drop << endl;
		OutputFile << "    " << "run_times : " << run_times << endl;
		OutputFile << "    " << "num_user_per_sector : " << num_MS_persector << endl;
		OutputFile << "    " << "Indoor_TRxP : " << num_Indoor_TRxP << endl;
		OutputFile << "    " << "simple_num_BS : " << num_BS << endl;
		OutputFile << "    " << "scenario : " << TYPE << endl;
		OutputFile << "    " << "Scheduling_Type : " << Scheduling_Type << endl;
		OutputFile << endl;
		//OutputFile << " Sector  AVR = " << Total_Avr / num_drop << endl;
		//OutputFile << " Minimum AVR = " << Total_Min / num_drop << endl;
		//OutputFile << " Spectral Efficiency AVR (from MIMO capacity) = " << SF_avr / num_drop << endl;
		OutputFile.close();
	}
	*/

	

}

void scheduling_statistics()
{
	double sector_selected_ue_per_time = 0;
	double scheduled_ue_mcs_per_time = 0;
	double scheduled_ue_cqi_per_time = 0;
	double scheduled_ue_widebandSINR_per_time = 0;

	double sector_metric_per_time = 0;
	double sector_cqi_read_per_time = 0;
	double sector_cqi_avr_per_time = 0;

	// In single_cell_mode, only collect statistics for center BS (sectors 0-2)
	int end_sector = (single_cell_mode == 1) ? 3 : num_SECTORS;

	for(int sector_idx = 0; sector_idx < end_sector; sector_idx++)
	{
		for (int rb_idx = 0; rb_idx < num_rb; rb_idx++)
		{
			for (int stream = 0; stream < mx_ue_mumimo; stream++)
			{
				sector_selected_ue_per_time += sector[sector_idx].ppschedulewrite[rb_idx][stream].ue_selected; // 스케줄된 ue의 index, 의미없지만 어차피 검증용이니 뽑아봄
				scheduled_ue_mcs_per_time += sector[sector_idx].ppschedulewrite[rb_idx][stream].mcs_selected; // 스케줄된 ue의 mcs index
				scheduled_ue_cqi_per_time += sector[sector_idx].ppschedulewrite[rb_idx][stream].cqi_selected; // 스케줄된 ue의 cqi 
				
				double linear_signal = dBm2linear(links[sector[sector_idx].ppschedulewrite[rb_idx][stream].ue_selected].str_signal);
				double linear_interference = dBm2linear(links[sector[sector_idx].ppschedulewrite[rb_idx][stream].ue_selected].interference);
				linear_interference = linear_interference + noise;

				scheduled_ue_widebandSINR_per_time += (linear_signal / linear_interference); // 스케줄된 ue의 wideband SINR을 계산
			}
		}
	}
	for(int sector_idx = 0; sector_idx < end_sector; sector_idx++)
	{
		for (int ue_idx = 0; ue_idx < (int)sector[sector_idx].ue_in_control.size(); ue_idx++)
		{
			for (int rb_idx = 0; rb_idx < num_rb; rb_idx++)
			{
				sector_metric_per_time   += sector[sector_idx].METRIC[ue_idx][rb_idx];
				sector_cqi_read_per_time += sector[sector_idx].CQI_read[ue_idx][rb_idx];
				sector_cqi_avr_per_time  += sector[sector_idx].CQI_AVR[ue_idx][rb_idx];
			}
		}

	}

	sector_selected_ue += sector_selected_ue_per_time / (num_SECTORS * num_rb * mx_ue_mumimo);
	scheduled_ue_mcs   += scheduled_ue_mcs_per_time / (num_SECTORS * num_rb * mx_ue_mumimo);
	scheduled_ue_cqi   += scheduled_ue_cqi_per_time / (num_SECTORS * num_rb * mx_ue_mumimo);
	scheduled_ue_widebandSINR += scheduled_ue_widebandSINR_per_time / (num_SECTORS * num_rb * mx_ue_mumimo);

	//cout << "scheduled_ue_mcs : " << scheduled_ue_mcs / (t - SCHEDULE_DELAY) << endl;
	//cout << "scheduled_ue_cqi : " << scheduled_ue_cqi / (t - SCHEDULE_DELAY) << endl;
    //cout << "scheduled_ue_widebandSINR : " << linear2dB(scheduled_ue_widebandSINR / (t - SCHEDULE_DELAY)) << endl;

	sector_metric   += sector_metric_per_time   / (num_SECTORS * num_MS * num_rb);
	sector_cqi_read += sector_cqi_read_per_time / (num_SECTORS * num_MS * num_rb);
	sector_cqi_avr  += sector_cqi_avr_per_time  / (num_SECTORS * num_MS * num_rb);
}

void measure_statistics()
{
	
	double ue_effective_sinr_per_time = 0;
	double ue_num_traffic_per_time = 0;
	double ue_info_bits_per_time = 0;
	double ue_bler_per_time = 0;


	for(int ms_idx = 0; ms_idx < num_MS; ms_idx++)
	{
		// 스케줄된 ue 기준
		if (ms[ms_idx].rb_indices.size() > 0)
		{
			ue_effective_sinr_per_time += ms[ms_idx].ESINR_linear;
			ue_info_bits_per_time += ms[ms_idx]._info_bits;
			ue_bler_per_time += ms[ms_idx].ue_BLER_Value;
		}

	}

	ue_effective_sinr += ue_effective_sinr_per_time / num_MS;
	//ue_num_traffic += ue_num_traffic_per_time / num_MS;
	ue_info_bits += ue_info_bits_per_time / num_MS;
	ue_bler += ue_bler_per_time / num_MS;
}