#include "Initiallization.h"


#include <cstring>
#include <cstdlib>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

void Set_BS_Location();
void Initialize_CHANNEL();
Real Get_distance(LOCATION, LOCATION);

void CQI_Map_creation_Map_connection();
void Schedule_Map_creation_Connection();
void Type1_Codebook_Gen();
void Type2_Codebook_Gen();
void ClearMAPs2StartADrop();
void Link_configuration();
void Get_CouplingLoss();

void Set_antenna_location_vector();
LOCATION3D Transpose_LCS_to_GCS_location(Real, Real, Real, LOCATION3D);
ComplexReal  Generate_vertical_virtualization_weight(Real theta, Real dv, Real k, Real K);
ComplexReal  Generate_horizontal_virtualization_weight(Real theta, Real phi, Real dH, Real l, Real L);
void Generate_bs_2D_DFT_beam_precoder();
void Declare_Precalc_Parameter();

void Do_precalculation();
void Load_precalculation_file();
void Print_DAT_FILE();
void Do_FT_precalc();
void Do_PMI_precalc(int );
void Do_SINR_estimate_precalc(int);
void Do_CH_COEF_SUM_precalc(int);

void Load_CHIR_precalc_file();
void Load_PMI_precalc_file();
void Load_SINR_estimate_precalc_file();
void Load_CH_COEF_SUM_precalc_file();

void Set_Parameter(int scenario)
{

	/*=====================================================================================
	INPUT: scenario index (1=InH, 2=UMi, 3=UMa, 4=RMa, 5=SMa, 6=3D-UMi, 7=3D-UMa, 8=3D-UMa-H)
	OUTPUT: scenario parameter
	====================================================================================*/
	TYPE = scenario;
	
	if (TYPE == 11) // InH
	{
		if (Configuration_Type == 0) ///// InH config A
		{
			///// IMT-2020 EVal, Page 19
			num_LINK = num_MS;

			//W = 20;//street width
			//h = 20;//avg. building height\
			user_speed = 3;  //[3km/h]
			Total_BS_Tx_power = 21; //[dBm for 20MHz]
			//carrier_freq = 4000000000.; //[4GHz] 4 * 10^9
			Wavelength = light_speed / carrier_freq;
			bs_maxpower = 21; //[dBm for 20MHz] =Total BS TX power
			ms_maxpower = 23; //[dBm UE power class]
			BS_noisefig = 5;
			MS_noisefig = 7;
			ue_antenna_element_gain = 0;
			max_antgain = 5; //[dBi]  BS antenna element gain

			//bandwidth = 10000000.; //[10MHz]  10 * 10^6
			//noise = dBm2linear(thermal_noise) * 15000; //[linear]  // 15000 = UE bandwidth (15kHz)
			//noise = dBm2linear(thermal_noise + MS_noisefig) * bandwidth; //[linear]  
			noise = dBm2linear(thermal_noise + (10. * log10(bandwidth)) + MS_noisefig); //[linear]  

			bs_height = 3.; //[m]
			ms_height = 1.5;
			inter_site_distance = 20;
			min_distance = 0.; // not decided in M.2412
			ANGLE_tilt = 0.; // change by calibration config

			num_floor = 1;
			num_propagation_condition = 2; // LOS / NLOS
			///////////////////////////////////////////////////
		}
		else if (Configuration_Type == 1) ///// InH config B
		{
			///// IMT-2020 EVal, Page 19

			num_LINK = num_MS;

			//W = 20;//street width
			//h = 20;//avg. building height

			user_speed = 3;  //[3km/h]  
			Total_BS_Tx_power = 20; //[dBm for 40MHz]
			//carrier_freq = 30000000000.; //[30GHz]
			Wavelength = light_speed / carrier_freq;
			bs_maxpower = 20; //[dBm for 40MHz] =Total BS TX power
			ms_maxpower = 23; //[dBm UE power class]
			BS_noisefig = 7;
			MS_noisefig = 10;
			ue_antenna_element_gain = 5;
			max_antgain = 5; //[dBi]  BS antenna element gain

			//bandwidth = 40000000.; //[40MHz]  40 * 10 ^ 6
			//noise = dBm2linear(thermal_noise) * 15000; //[linear]  // 15000 = UE bandwidth (15kHz)
			//noise = dBm2linear(thermal_noise + MS_noisefig) * bandwidth; //[linear]  
			noise = dBm2linear(thermal_noise + (10. * log10(bandwidth)) + MS_noisefig); //[linear]  

			bs_height = 3.; //[m]
			ms_height = 1.5;
			inter_site_distance = 20;
			min_distance = 0.; // not decided in M.2412
			ANGLE_tilt = 0.; // change by calibration config

			num_floor = 1;
			num_propagation_condition = 2; // LOS / NLOS

			///////////////////////////////////////////////////
		}
		else if (Configuration_Type == 2) ///// InH config C
		{
			///// IMT-2020 EVal, Page 19

			num_LINK = num_MS;

			//W = 20;//street width
			//h = 20;//avg. building height

			user_speed = 3;  //[3km/h]  

			// jhnoh 220902
			Total_BS_Tx_power = 30; //[dBm for 40MHz]
			
			//carrier_freq = 70000000000.; //[70GHz]
			Wavelength = light_speed / carrier_freq;

			// jhnoh 220902 ----------------------------------------------
			// BS POWER CONTROL
			bs_maxpower = 28; //18; //[dBm for 40MHz] =Total BS TX power
			ms_maxpower = 21; //[dBm UE power class]
			BS_noisefig = 7;
			MS_noisefig = 10;
			ue_antenna_element_gain = 5;
			max_antgain = 5; //[dBi]  BS antenna element gain

			//bandwidth = 40000000.; //[40MHz]  40 * 10 ^ 6
			//noise = dBm2linear(thermal_noise) * 15000; //[linear]  // 15000 = UE bandwidth (15kHz)
			//noise = dBm2linear(thermal_noise + MS_noisefig) * bandwidth; //[linear]  
			noise = dBm2linear(thermal_noise + (10. * log10(bandwidth)) + MS_noisefig); //[linear]  
			
			bs_height = 3.; //[m]
			ms_height = 1.5;

			inter_site_distance = 5; //20;

			min_distance = 0.; // not decided in M.2412
			ANGLE_tilt = 0.; // change by calibration config

			num_floor = 1;
			num_propagation_condition = 2; // LOS / NLOS
			///////////////////////////////////////////////////
		}
		else
		{
			cout << "Configuration Type Error !!" << endl;
			getchar();
		}

	}
	else if (TYPE == 12) // Dense_Urban
	{
		if (Configuration_Type == 0) ///// Dense Urban config A
		{
			///// IMT-2020 EVal, Page 20
			num_LINK = num_MS;

			user_speed = 3;  //[3km/h] Indoor user // outdoor user speed -> decided later
			Total_BS_Tx_power = 41; //[dBm for 10MHz]
			//carrier_freq = 4000000000.; //[4GHz]
			Wavelength = light_speed / carrier_freq;

			bs_maxpower = 41; //[dBm for 10MHz] =Total BS TX power
			ms_maxpower = 23; //[dBm UE power class]
			BS_noisefig = 5;
			MS_noisefig = 7;
			ue_antenna_element_gain = 0;
			max_antgain = 8; //[dBi]  BS antenna element gain

			//bandwidth = 10000000.; //[10MHz]  10 * 10 ^ 6
			//noise = dBm2linear(thermal_noise) * 15000; //[linear]  // 15000 = UE bandwidth (15kHz)
			//noise = dBm2linear(thermal_noise + MS_noisefig) * bandwidth; //[linear]  
			noise               = dBm2linear(thermal_noise + (10. * log10(bandwidth)) + MS_noisefig); //[linear]  
			bs_height           = 25.;    /// [m]
			macro_bs_height     = 25.;
            micro_bs_height     = 10.;			

			ms_height_in        = 0;   /// in LOS prob part
			ms_height_out       = 1.5;
			inter_site_distance = 500;

			min_distance        = 35.; // not decided in M.2412
			ANGLE_tilt          = 0.;  // change by calibration config

			num_floor                 = 8;
			num_propagation_condition = 3;     /// LOS, NLOS, OUT2IN

			//W = 20.; //street width
			//h = 20.; //avg. building height


			//ANGLE_tilt = 6.;
			////////////////////////////////////////////////
			//MAX_v_angle = atan2((bs_height - ms_height_out), min_distance) * 180 / pi;   // FOR Pre-Calculate v angle
		}
		else if (Configuration_Type == 1) ///// Dense Urban config B
		{
			///// IMT-2020 EVal, Page 20
			num_LINK = num_MS;

			user_speed = 3;  //[3km/h] Indoor user // outdoor user speed -> decided later
			Total_BS_Tx_power = 37; //[dBm for 40MHz]
			//carrier_freq = 30000000000.; //[30GHz]
			Wavelength = light_speed / carrier_freq;

			bs_maxpower    = 37; //[dBm for 40MHz] =Total BS TX power
			micro_bs_power = bs_maxpower - 7;

			ms_maxpower = 23; //[dBm UE power class]
			BS_noisefig = 7;
			MS_noisefig = 10;
			ue_antenna_element_gain = 5;
			max_antgain = 8; //[dBi]  BS antenna element gain

			// bandwidth = 40000000.; //[40MHz]  40 * 10 ^ 6
			// noise = dBm2linear(thermal_noise) * 15000; //[linear]  // 15000 = UE bandwidth (15kHz)
			// noise = dBm2linear(thermal_noise + MS_noisefig) * bandwidth; //[linear]  
			noise = dBm2linear(thermal_noise + (10. * log10(bandwidth)) + MS_noisefig); //[linear]  

			bs_height           = 25.; //[ m]
			macro_bs_height     = 25.; 
            micro_bs_height     = 10.;

			ms_height_in        = 0;   // in LOS prob part
			ms_height_out       = 1.5;
			inter_site_distance = 500; // 200; //32;//57.7;
			min_distance        = 35.; // not decided in M.2412
			ANGLE_tilt          = 0.;  // change by calibration config

			num_floor = 8;
			num_propagation_condition = 3;
		}
		else if (Configuration_Type == 2) ///// Dense Urban config C Macro - 30GHz
		{
			/*------------- IMT-2020 EVal, Page 20 ------------- */ 
			/*---------------   Macro cell case  --------------- */ 
			num_LINK = num_MS;

			//// [3km/h] Indoor user // outdoor user speed -> decided later
			user_speed = 3;  
			Total_BS_Tx_power = 37; //[dBm for 40MHz]
			//// carrier_freq = 30000000000.; //[30GHz]
			Wavelength = light_speed / carrier_freq;

			bs_maxpower = 37; //[dBm for 40MHz]=Total BS TX power
			micro_bs_power = bs_maxpower - 7;

			ms_maxpower = 23; //[dBm UE power class]
			BS_noisefig = 7;
			MS_noisefig = 10;
			ue_antenna_element_gain = 5;
			max_antgain = 8; //[dBi]  BS antenna element gain

			//bandwidth = 40000000.; //[40MHz]  40 * 10 ^ 6
			//noise = dBm2linear(thermal_noise) * 15000; //[linear]  // 15000 = UE bandwidth (15kHz)
			//noise = dBm2linear(thermal_noise + MS_noisefig) * bandwidth; //[linear]  
			noise = dBm2linear(thermal_noise + (10. * log10(bandwidth)) + MS_noisefig); //[linear]  

			bs_height = 25.; //[m]
			macro_bs_height     = 25.;
            micro_bs_height     = 10.;

			ms_height_in = 0;   /// in LOS prob part
			ms_height_out = 1.5;
			inter_site_distance = 200;   // Dense Urban-eMBB macro layer ISD=200m (TR 38.913 Table 6.1.4-1; agency 6G req 2026-06). Was 500m (generic UMa).
			min_distance = 35.; // not decided in M.2412
			ANGLE_tilt = 0.; // change by calibration config

			num_floor = 8;
			num_propagation_condition = 3;

		}
		else if (Configuration_Type == 3) ///// Dense Urban config C Micro - 4GHz
		{
			///// IMT-2020 EVal, Page 20
			///// Micro  cell case
			num_LINK = num_MS;

			user_speed = 3;  //[3km/h] Indoor user // outdoor user speed -> decided later
			Total_BS_Tx_power = 44; //[dBm for 10MHz]
			//carrier_freq = 4000000000.; //[4GHz]
			Wavelength = light_speed / carrier_freq;

			bs_maxpower = 44; //[dBm for 10MHz] = Total BS TX power
			micro_bs_power = bs_maxpower - 7;

			ms_maxpower = 23; //[dBm UE power class]
			BS_noisefig = 7;
			MS_noisefig = 10;
			ue_antenna_element_gain = 0;  // UMi: isotropic UE antenna (TR 38.901 Table 7.8-2A, calibration)
			max_antgain = 8; //[dBi]  BS antenna element gain

			//bandwidth = 10000000.; //[10MHz]  10 * 10 ^ 6
			//noise = dBm2linear(thermal_noise) * 15000; //[linear]  // 15000 = UE bandwidth (15kHz)
			//noise = dBm2linear(thermal_noise + MS_noisefig) * bandwidth; //[linear]  
			noise = dBm2linear(thermal_noise + (10. * log10(bandwidth)) + MS_noisefig); //[linear]  

			bs_height = 10.; //[m]
			macro_bs_height     = 10.;
            micro_bs_height     = 10.;

			ms_height_in = 0;   /// in LOS prob part
			ms_height_out = 1.5;
			inter_site_distance = 200;   // UMi Street Canyon: ISD = 200m (TR 38.901 Table 7.8-2A)

			min_distance = 10.; // UMi: minimum BS-UT distance = 10m (TR 38.901 Table 7.8-2A)
			ANGLE_tilt = 0.; // change by calibration config

			num_floor = 8;
			num_propagation_condition = 3;     /// LOS, NLOS, OUT2IN

		}
		else if (Configuration_Type == 4)
		{
			///// IMT-2020 EVal, Page 20
			///// Macro  cell case
			num_LINK = num_MS;

			user_speed = 3;  //[3km/h] Indoor user // outdoor user speed -> decided later
			Total_BS_Tx_power = 37 +19; //[dBm for 40MHz]
			//carrier_freq = 30000000000.; //[30GHz]
			Wavelength = light_speed / carrier_freq;

			bs_maxpower = 37+19; //[dBm for 40MHz]=Total BS TX power
			micro_bs_power = bs_maxpower - 7;

			ms_maxpower = 23; //[dBm UE power class]
			BS_noisefig = 7;
			MS_noisefig = 10;
			ue_antenna_element_gain = 5;
			max_antgain = 8; //[dBi]  BS antenna element gain

			//bandwidth = 40000000.; //[40MHz]  40 * 10 ^ 6
			//noise = dBm2linear(thermal_noise) * 15000; //[linear]  // 15000 = UE bandwidth (15kHz)
			//noise = dBm2linear(thermal_noise + MS_noisefig) * bandwidth; //[linear]  
			noise = dBm2linear(thermal_noise + (10. * log10(bandwidth)) + MS_noisefig); //[linear]  

			bs_height = 25.; //[m]
			macro_bs_height     = 25.;
            micro_bs_height     = 10.;

			ms_height_in = 0;   /// in LOS prob part
			ms_height_out = 1.5;
			inter_site_distance = 200;
			min_distance = 10.; // not decided in M.2412
			ANGLE_tilt = 0.; // change by calibration config

			num_floor = 8;
			num_propagation_condition = 3;
		}
		else
		{
			cout << "Configuration Type Error !!" << endl;
			getchar();
		}

	}
	else if (TYPE == 13) // RMa
	{
		if (Configuration_Type == 0) ///// RMa config A
		{
			///// IMT-2020 EVal, Page 21

			num_LINK = num_MS;

			user_speed = 3;  //[3km/h] Indoor user // outdoor user speed -> decided later
			Total_BS_Tx_power = 46; //[dBm for 10MHz]
			//carrier_freq = 700000000.; //[700MHz] 700 * 10^6
			Wavelength = light_speed / carrier_freq;

			bs_maxpower = 46; //[dBm for 10MHz] =Total BS TX power
			micro_bs_power = bs_maxpower - 7;

			ms_maxpower = 23; //[dBm UE power class]
			BS_noisefig = 5;
			MS_noisefig = 7;
			ue_antenna_element_gain = 0;
			max_antgain = 8; //[dBi]  BS antenna element gain

			//bandwidth = 10000000.; //[10MHz]  10 * 10 ^ 6
			//noise = dBm2linear(thermal_noise) * 15000; //[linear]  // 15000 = UE bandwidth (15kHz)
			//noise = dBm2linear(thermal_noise + MS_noisefig) * bandwidth; //[linear] 
			noise = dBm2linear(thermal_noise + (10. * log10(bandwidth)) + MS_noisefig); //[linear]  

			bs_height = 35.; //[m]
			macro_bs_height     = 35.;
            micro_bs_height     = 10.;

			ms_height_out = 1.5;
			ms_height_in = 1.5;

			inter_site_distance = 1732;

			min_distance = 10.; // not decided in M.2412
			ANGLE_tilt = 0.; // change by calibration config

			num_floor = 1;
			num_propagation_condition = 3;

		}
		else if (Configuration_Type == 1) ///// RMa config B
		{
			///// IMT-2020 EVal, Page 21

			num_LINK = num_MS;

			user_speed = 3;  //[3km/h] Indoor user // outdoor user speed -> decided later
			Total_BS_Tx_power = 46; //[dBm for 10MHz]
			//carrier_freq = 4000000000.; //[4GHz]
			Wavelength = light_speed / carrier_freq;

			bs_maxpower = 46 + 30; //[dBm for 10MHz] =Total BS TX power
			micro_bs_power = bs_maxpower - 7;

			ms_maxpower = 23; //[dBm UE power class]
			BS_noisefig = 5;
			MS_noisefig = 7;
			ue_antenna_element_gain = 5;
			max_antgain = 8; //[dBi]  BS antenna element gain

			//bandwidth = 10000000.; //[10MHz]  10 * 10 ^ 6
			//noise = dBm2linear(thermal_noise) * 15000; //[linear]  // 15000 = UE bandwidth (15kHz)
			//noise = dBm2linear(thermal_noise + MS_noisefig) * bandwidth; //[linear] 
			noise = dBm2linear(thermal_noise + (10. * log10(bandwidth)) + MS_noisefig); //[linear]  

			bs_height = 35.; //[m]
			macro_bs_height     = 35.;
            micro_bs_height     = 10.;

			ms_height_out = 1.5;
			ms_height_in = 1.5;

			inter_site_distance = 1732;

			min_distance = 10.; // not decided in M.2412
			ANGLE_tilt = 0.; // change by calibration config

			num_floor = 1;
			num_propagation_condition = 3;

		}
		else if (Configuration_Type == 2) ///// RMa config C
		{
			///// IMT-2020 EVal, Page 21

			num_LINK = num_MS;


			user_speed = 3;  //[3km/h] Indoor user // outdoor user speed -> decided later
			Total_BS_Tx_power = 46; //[dBm for 10MHz]
			//carrier_freq = 700000000.; //[700MHz]
			Wavelength = light_speed / carrier_freq;

			bs_maxpower = 46; //[dBm for 10MHz] =Total BS TX power
			micro_bs_power = bs_maxpower - 7;
			
			ms_maxpower = 23; //[dBm UE power class]
			BS_noisefig = 5;
			MS_noisefig = 7;
			ue_antenna_element_gain = 0;
			max_antgain = 8; //[dBi]  BS antenna element gain

			//bandwidth = 10000000.; //[10MHz]   10 * 10 ^ 6
			//noise = dBm2linear(thermal_noise) * 15000; //[linear]  // 15000 = UE bandwidth (15kHz)
			//noise = dBm2linear(thermal_noise + MS_noisefig) * bandwidth; //[linear] 
			noise = dBm2linear(thermal_noise + (10. * log10(bandwidth)) + MS_noisefig); //[linear]  

			bs_height = 35.; //[m]
			macro_bs_height     = 35.;
			micro_bs_height     = 10.;

			ms_height_out = 1.5;
			ms_height_in = 1.5;

			inter_site_distance = 6000;

			min_distance = 10.; // not decided in M.2412
			ANGLE_tilt = 0.; // change by calibration config

			num_floor = 1;
			num_propagation_condition = 3;

		}
		else
		{
			cout << "Configuration Type Error !!" << endl;
			getchar();
		}

	}

	BS_dH = BS_dH * Wavelength;
	BS_dV = BS_dV * Wavelength;
	BS_dgH = BS_dgH * Wavelength;
	BS_dgV = BS_dgV * Wavelength;

	MS_dH = MS_dH * Wavelength;
	MS_dV = MS_dV * Wavelength;
	MS_dgH = MS_dgH * Wavelength;
	MS_dgV = MS_dgV * Wavelength;

	// Override with cfg values if provided
	bool bs_tx_overridden = false;
	bool ut_nf_overridden = false;
	if (cfg_BS_Tx_Power > -9998) {
		Total_BS_Tx_power = cfg_BS_Tx_Power;
		bs_maxpower = cfg_BS_Tx_Power;
		bs_tx_overridden = true;
	}
	if (cfg_UT_Noise_Figure > -9998) {
		MS_noisefig = cfg_UT_Noise_Figure;
		noise = dBm2linear(thermal_noise + (10. * log10(bandwidth)) + MS_noisefig);
		ut_nf_overridden = true;
	}

	// cfg override for UE antenna element gain (e.g., 0=isotropic, 5=directional)
	bool ue_gain_overridden = false;
	if (cfg_UE_antenna_element_gain > -9998) {
		ue_antenna_element_gain = cfg_UE_antenna_element_gain;
		ue_gain_overridden = true;
	}

	// cfg override for inter-site distance (e.g., UMa ISD 200 vs 500 m). -1 = use scenario default.
	if (cfg_inter_site_distance > 0) {
		inter_site_distance = cfg_inter_site_distance;
	}
	// Handheld mode default: force directional path unless cfg explicitly overrides
	if (handheld_mode && !ue_gain_overridden) {
		ue_antenna_element_gain = 5;  // non-zero flag to trigger BF code paths
	}

	cout << "carrier_freq            : " << carrier_freq << endl;
	cout << "inter_site_distance     : " << inter_site_distance << endl;
	cout << "sector_to_ue_min_dist   : " << min_distance << endl;
	cout << "indoor_ue_speed         : " << user_speed << endl;
	cout << "h_BS                    : " << bs_height << endl;
	cout << "bs_tx_power             : " << bs_maxpower << " dBm" << (bs_tx_overridden ? " (overridden by cfg)" : "") << endl;
	cout << "bs_noise_figure         : " << BS_noisefig << " dB" << endl;
	cout << "ue_noise_figure         : " << MS_noisefig << " dB" << (ut_nf_overridden ? " (overridden by cfg)" : "") << endl;
	cout << "bs_antenna_element_gain : " << max_antgain << " dBi" << endl;
	cout << "ue_antenna_element_gain : " << ue_antenna_element_gain << " dBi";
	if (ue_gain_overridden) cout << " (overridden by cfg)";
	else if (handheld_mode) cout << " (handheld default)";
	cout << endl;
	cout << "thermal_noise_level     : " << thermal_noise << " dBm/Hz" << endl;
	cout << "simulation_bandwidth    : " << bandwidth << " Hz" << endl;

	//getchar();
	num_mTRP = 0;

	// 230616 mTRP
	bs = new BS[num_BS + num_mTRP];
	mTRP = &bs[num_BS];

	ms = new MS[num_MS];
	for (int ms_idx = 0; ms_idx < num_MS; ms_idx++)
		ms[ms_idx].self_idx = ms_idx;

	links       = new LINK  [num_LINK];
	sector      = new Sector[num_SECTORS + num_mTRP_SECTORS];

	for (int sec_idx = 0; sec_idx < num_SECTORS + num_mTRP_SECTORS; sec_idx++)
		sector[sec_idx].self_idx = sec_idx;

	// 230615 jhnoh
	int num_channel = num_BS;

	channel = new CHANNEL*[num_channel];
	for (int bs_idx = 0; bs_idx < num_channel; bs_idx++)
	{
		channel[bs_idx] = new CHANNEL[num_MS];
	}

	Set_BS_Location();
}


void Set_BS_Location()
{
	/*
	INPUT: scenario index 
	OUTPUT: BS location
	*/


	if (TYPE == 11) //// 5G InH
	{
		bs[0].loc.x = -50;
		bs[0].loc.y = 35;

		bs[1].loc.x = -30;   
		bs[1].loc.y = 35;

		bs[2].loc.x = -10;
		bs[2].loc.y = 35;

		bs[3].loc.x = 10;
		bs[3].loc.y = 35;

		bs[4].loc.x = 30;
		bs[4].loc.y = 35;

		bs[5].loc.x = 50;
		bs[5].loc.y = 35;

		bs[6].loc.x = -50;
		bs[6].loc.y = 15;

		bs[7].loc.x = -30;
		bs[7].loc.y = 15;

		bs[8].loc.x = -10;
		bs[8].loc.y = 15;

		bs[9].loc.x = 10;
		bs[9].loc.y = 15;

		bs[10].loc.x = 30;
		bs[10].loc.y = 15;

		bs[11].loc.x = 50;
		bs[11].loc.y = 15;

	}
	else  //// hexagonal scenario
	{
		Real cell_radius = inter_site_distance * (1 / sqrt(3));

		// ====================================================================
		// 1-SITE CONFIGURATION: Single BS with 3 interfering replica BSs
		// ====================================================================
		// To properly model inter-sector interference in 1-site mode, we place
		// 4 base stations total:
		//   - BS 0: Central site at (0, 0) with 3 sectors (0°, 120°, 240°)
		//   - BS 1: Interferer for Sector 0 (0°) - placed at distance in 0° direction
		//   - BS 2: Interferer for Sector 1 (120°) - placed at distance in 120° direction
		//   - BS 3: Interferer for Sector 2 (240°) - placed at distance in 240° direction
		//
		// Each interfering BS has the same 3-sector configuration as the central BS,
		// but positioned to create realistic interference for the target sector.
		// ====================================================================
		if (simple_num_BS == 1) {
			// Central BS at origin
			bs[0].loc.x = 0;
			bs[0].loc.y = 0;

			// The actual num_BS will be set to 4 to include interferers
			// This is handled in the configuration logic
		}
		else {
			// 7-site and 19-site configurations
			bs[0].loc.x = 0;
			bs[1].loc.x = 0;
			bs[2].loc.x = -1.5 * cell_radius;
			bs[3].loc.x = -1.5 * cell_radius;
			bs[4].loc.x = 0;
			bs[5].loc.x = 1.5 * cell_radius;
			bs[6].loc.x = 1.5 * cell_radius;
		}
		
		if (simple_num_BS == 19) {
			bs[7].loc.x = 0;
			bs[8].loc.x = -1.5 * cell_radius;
			bs[9].loc.x = -3. * cell_radius;
			bs[10].loc.x = -3. * cell_radius;
			bs[11].loc.x = -3. * cell_radius;
			bs[12].loc.x = -1.5 * cell_radius;
			bs[13].loc.x = 0;
			bs[14].loc.x = 1.5 * cell_radius;
			bs[15].loc.x = 3. * cell_radius;
			bs[16].loc.x = 3. * cell_radius;
			bs[17].loc.x = 3. * cell_radius;
			bs[18].loc.x = 1.5 * cell_radius;
		}

		//Location y
		bs[0].loc.y = 0;
		bs[1].loc.y = inter_site_distance;
		bs[2].loc.y = inter_site_distance * (0.5);
		bs[3].loc.y = inter_site_distance * (-0.5);
		bs[4].loc.y = inter_site_distance * (-1);
		bs[5].loc.y = inter_site_distance * (-0.5);
		bs[6].loc.y = inter_site_distance * (0.5);
		
		if (simple_num_BS == 19) {
			bs[7].loc.y = inter_site_distance * (2.);
			bs[8].loc.y = inter_site_distance * (1.5);
			bs[9].loc.y = inter_site_distance * (1.);
			bs[10].loc.y = 0;
			bs[11].loc.y = inter_site_distance * (-1);
			bs[12].loc.y = inter_site_distance * (-1.5);
			bs[13].loc.y = inter_site_distance * (-2);
			bs[14].loc.y = inter_site_distance * (-1.5);
			bs[15].loc.y = inter_site_distance * (-1);
			bs[16].loc.y = 0;
			bs[17].loc.y = inter_site_distance * (1.);
			bs[18].loc.y = inter_site_distance * (1.5);
		}

	}
}

void ClearDrop()
{
	for(int ms_idx = 0; ms_idx < num_MS; ms_idx++)
	{
		ms[ms_idx].Reset2Default();
		links[ms_idx].Reset2Default();
	}

	for(int sec_idx = 0; sec_idx < num_SECTORS + num_mTRP_SECTORS; sec_idx++)
	{
		sector[sec_idx].Reset2Default();
	}

	for(int bs_idx = 0; bs_idx < num_BS; bs_idx++)
		for(int ms_idx = 0; ms_idx < num_MS; ms_idx++)
			channel[bs_idx][ms_idx].Reset2Default();



}

void Initialdrop()
{

	///Real cell_radius = inter_site_distance / sqrt(3.);


	bool a = 0; //true ==> 원래 크기 faluse==> 크기 축소
	Real n = 20/5.; // 크기 축소 factor

	///////////////////////////////////////////////////////////////// 5G InH
	if (TYPE == 11)          
	{
		// In single_cell_mode, only create users for center BS (BS 0)
		int num_bs_with_users = (single_cell_mode == 1) ? 1 : num_BS;

		for (int bs_idx = 0; bs_idx < num_bs_with_users; bs_idx++) /// 12
		{
			for (int ms_idx = 0; ms_idx < num_Indoor_TRxP * num_MS_persector; ms_idx++)  /// 10 or 30
			{
				ms[ms_idx + bs_idx * num_Indoor_TRxP * num_MS_persector].nearest_bs_idx = bs_idx;   // ms 0~9 -> bs 0, ms 10~19 -> bs 1 ///// or ms 0~29 -> bs 0, ms 30~59 -> bs 1
				do
				{
					ms[ms_idx + bs_idx * num_Indoor_TRxP * num_MS_persector].loc.x = (bs[bs_idx].loc.x - 10) + 20. * randnum.u();

					if (bs_idx < 6)   /// 6 BS  on Top line  (0 ~ 5)
					{
						ms[ms_idx + bs_idx * num_Indoor_TRxP * num_MS_persector].loc.y = (bs[bs_idx].loc.y - 10) + 25. * randnum.u();
					}
					else  /// 6 BS on Bottom line (6 ~ 11)
					{
						ms[ms_idx + bs_idx * num_Indoor_TRxP * num_MS_persector].loc.y = (bs[bs_idx].loc.y - 15) + 25. * randnum.u();
					}

					ms[ms_idx + bs_idx * num_Indoor_TRxP * num_MS_persector].moving_direction = 2 * 180 * randnum.u() - 180;
					ms[ms_idx + bs_idx * num_Indoor_TRxP * num_MS_persector].speed = user_speed * 1000 / 3600;
					ms[ms_idx + bs_idx * num_Indoor_TRxP * num_MS_persector].Indoor = true;

				} while (Get_distance(ms[ms_idx + bs_idx * num_Indoor_TRxP * num_MS_persector].loc, bs[bs_idx].loc) < min_distance);

			}
		}

		if(a==0)
		{
			for (int bs_idx = 0; bs_idx < num_BS; bs_idx++)
			{
				bs[bs_idx].loc.x = bs[bs_idx].loc.x/n;
				bs[bs_idx].loc.y = bs[bs_idx].loc.y/n;
				for (int ms_idx = 0; ms_idx < num_Indoor_TRxP * num_MS_persector; ms_idx++)
				{
					ms[ms_idx + bs_idx * num_Indoor_TRxP * num_MS_persector].loc.x  = ms[ms_idx + bs_idx * num_Indoor_TRxP * num_MS_persector].loc.x/n;
					ms[ms_idx + bs_idx * num_Indoor_TRxP * num_MS_persector].loc.y  = ms[ms_idx + bs_idx * num_Indoor_TRxP * num_MS_persector].loc.y/n;
				}

			}	
		}
	
		
	}
	////////////////////////////////////////////////////////////////// hexagonal scenario
	else
	{
		// In single_cell_mode, only create users for center BS (BS 0)
		int num_bs_with_users = (single_cell_mode == 1) ? 1 : num_BS;

		for (int bs_idx = 0; bs_idx < num_bs_with_users; bs_idx++)
		{
			for (int sector_idx = 0; sector_idx < 3; sector_idx++)
			{
				for (int ms_idx = bs_idx*num_MS_persector * 3 + num_MS_persector*(sector_idx); ms_idx < bs_idx*num_MS_persector * 3 + num_MS_persector*(sector_idx+1); ms_idx++) // ms_idx = 0~9 -> bs[0],sector1, 10~19 -> bs[0],sector2, 20~29 -> bs[0],sector3, 30~39 -> bs[1],sector1 ......
				{
					ms[ms_idx].Configuration(ms_idx, sector_idx);
					ms[ms_idx].nearest_bs_idx = bs_idx;
					ms[ms_idx].loc.x = bs[bs_idx].loc.x + ms[ms_idx].loc.x;
					ms[ms_idx].loc.y = bs[bs_idx].loc.y + ms[ms_idx].loc.y;
					ms[ms_idx].moving_direction = 2 * 180 * randnum.u() - 180; //travel angle (elevation)
					ms[ms_idx].moving_direction_azimuth = 2 * 180 * randnum.u() - 180; // travel azimuth angle
					ms[ms_idx].speed = user_speed * 1000 / 3600;

					if (TYPE == 12)   /// imt 2020 Dense urban
					{
						if (ms[ms_idx].Indoor == false)  //// outdoor speed
						{
							ms[ms_idx].speed = 30. * 1000. / 3600.;
						}
					}
					else if (TYPE == 13)  //// imt 2020 rural
					{
						if (ms[ms_idx].Indoor == false)  /// outdoor speed
						{
							if (Configuration_Type == 2) /// configuration C
							{
								if (randnum.u() < 0.5) 
								{
									ms[ms_idx].speed = user_speed * 1000. / 3600.;  /// outdoor pedestrian 3km/h
								}
								else
								{
									ms[ms_idx].speed = 30. * 1000. / 3600.;  /// outdoor in car
								}
							}
							else
							{
								ms[ms_idx].speed = 120. * 1000. / 3600.;
							}
						}
					}
				}
			}
		}
	}

	cout << "initialdrop.......................DONE" << endl;
	//cout << "  " << endl;
}

void Initialize_CHANNEL()
{
	//검증용 인석 추가
	sector_selected_ue        = 0;
	scheduled_ue_mcs          = 0;
	scheduled_ue_cqi          = 0;
	scheduled_ue_widebandSINR = 0;

	ue_effective_sinr        = 0;
	ue_num_traffic           = 0;
	ue_info_bits             = 0;
	ue_bler                  = 0;

	sector_metric           = 0;
	sector_cqi_read         = 0;
	sector_cqi_avr          = 0;

	Real cell_radius = inter_site_distance / sqrt(3.);
	Real distance;

	LOCATION temp_bs;
	LOCATION nearest_bs;
	
	//////////////////////////////////////////////////////////////////////////////// 5G InH
	if (TYPE == 11) 
	{
		for (int ms_idx = 0; ms_idx < num_MS; ms_idx++)
		{
			for (int bs_idx = 0; bs_idx < num_BS + num_mTRP; bs_idx++)
			{
				nearest_bs.x                      = bs[bs_idx].loc.x;
				nearest_bs.y                      = bs[bs_idx].loc.y;
				
				distance                          = Get_distance(ms[ms_idx].loc, nearest_bs);

				channel[bs_idx][ms_idx].self_bs.x = nearest_bs.x;
				channel[bs_idx][ms_idx].self_bs.y = nearest_bs.y;
				channel[bs_idx][ms_idx].self_ms.x = ms[ms_idx].loc.x;
				channel[bs_idx][ms_idx].self_ms.y = ms[ms_idx].loc.y;

				channel[bs_idx][ms_idx].Set_channel(distance, ms[ms_idx].Indoor, bs_idx, ms_idx, channel[bs_idx][ms_idx].self_bs, channel[bs_idx][ms_idx].self_ms);
			}
		}
	}
	////////////////////////////////////////////////////////////////////////////////// 2 tier wrap around Hexagonal grid , 3 Sector per 1 BS
	else 
	{
		//if (simple_num_BS == 7)
		//{
			for (int ms_idx = 0; ms_idx < num_MS; ms_idx++)
			{
				int ms_to_bs_wrap_idx = (ms_idx / (3*num_MS_persector)); // -> ms_idx�� bs_idx��

				LOCATION center_cell_view;
				center_cell_view.x = ms[ms_idx].loc.x - bs[ms_to_bs_wrap_idx].loc.x;
				center_cell_view.y = ms[ms_idx].loc.y - bs[ms_to_bs_wrap_idx].loc.y;

				for (int wrap_bs_idx = 0; wrap_bs_idx < simple_num_BS; wrap_bs_idx++)
				{
					distance = Get_distance(center_cell_view, bs[wrap_bs_idx].loc);
					int bs_idx = simple_wrap_mat[ms_to_bs_wrap_idx][wrap_bs_idx];

					channel[bs_idx][ms_idx].self_bs.x = bs[wrap_bs_idx].loc.x + bs[ms_to_bs_wrap_idx].loc.x;
					channel[bs_idx][ms_idx].self_bs.y = bs[wrap_bs_idx].loc.y + bs[ms_to_bs_wrap_idx].loc.y;
					channel[bs_idx][ms_idx].self_ms.x = ms[ms_idx].loc.x;
					channel[bs_idx][ms_idx].self_ms.y = ms[ms_idx].loc.y;

					channel[bs_idx][ms_idx].Set_channel(distance, ms[ms_idx].Indoor, bs_idx, ms_idx, channel[bs_idx][ms_idx].self_bs, channel[bs_idx][ms_idx].self_ms);
				}


			}	
		//}
	}
}

// ====================================================================
// Row/Column analog beam allocation (population vote), once per drop.
//
// Replaces the per-link "every UE gets its own best analog beam" optimism
// (physically impossible on one subarray-partitioned panel in MU-MIMO) with
// one shared zenith beam per vertical port row, allocated from the attached
// UEs' RSRP votes:
//   1) each UE votes for up to K zenith beams within X dB of its best beam
//      (measured on its own azimuth / UE-beam slice of the cached RSRP table);
//   2) votes are apportioned over the BS_Mp rows (largest remainder);
//   3) beams are placed sorted by zenith angle, steepest downtilt on the top
//      row (mi = BS_Mp-1) downwards, so equal beams sit on adjacent rows;
//   4) az_mode==2 repeats 1-3 with azimuth beams over the BS_Np columns.
// Serving-link bookkeeping (link_RSRP / link_antgain / str_signal /
// static_gain) is then refreshed with the row-mixed RSRP so downstream
// consumers (Get_CouplingLoss, geometry, RI) stay consistent. Cell
// association itself keeps the best-beam result (broadcast-beam approx).
//
// Deterministic: no RNG draws; UEs are scanned by index (NOT via
// sector[].ue_in_control, whose order is nondeterministic — it is built by
// unsynchronized push_back inside the OMP "Get RSRP" loop).
// ====================================================================

// One UE's vote: +1 to each of at most row_beam_max_cand beams whose RSRP is
// within row_beam_x_db of the best, strongest first (RSRP tie → lower index).
static void RowBeam_Vote(const Real* rsrp_db, int n, Real best_db, Real* count)
{
	bool used[8] = {false,};
	for (int votes = 0; votes < row_beam_max_cand; votes++)
	{
		int pick = -1;
		for (int i = 0; i < n; i++)
		{
			if (used[i]) continue;
			if (rsrp_db[i] < best_db - row_beam_x_db) continue;
			if (pick < 0 || rsrp_db[i] > rsrp_db[pick]) pick = i;
		}
		if (pick < 0) break;
		used[pick] = true;
		count[pick] += REAL(1.0);
	}
}

// Largest-remainder apportionment of `total` slots proportional to count[].
// Remainder ties: larger count first, then lower index. Caller guarantees
// sum(count) > 0.
static void RowBeam_Apportion(const Real* count, int n, int total, int* alloc)
{
	Real sum = 0;
	for (int i = 0; i < n; i++) sum += count[i];

	Real frac[8];
	int given = 0;
	for (int i = 0; i < n; i++)
	{
		Real quota = count[i] * (Real)total / sum;
		alloc[i] = (int)floor(quota + REAL(1e-6));  // epsilon: exact-integer quotas
		frac[i]  = quota - (Real)alloc[i];
		given += alloc[i];
	}
	for (; given < total; given++)
	{
		int pick = 0;
		for (int i = 1; i < n; i++)
		{
			if (frac[i] > frac[pick] + REAL(1e-9)) { pick = i; continue; }
			if (frac[i] > frac[pick] - REAL(1e-9) && count[i] > count[pick]) pick = i;
		}
		alloc[pick] += 1;
		frac[pick] = -REAL(1.0);
	}
}

void Assign_Row_Beams()
{
	if (!row_beam_enable) return;

	const int V  = BS_Mp;
	const int H  = BS_Np;
	const int nZ = tilt_zenith_angle_LCS_size;
	const int nA = tilt_azimuth_angle_LCS_size;

	if (V > 16 || H > 16)
	{
		cout << "ERROR: Assign_Row_Beams supports at most 16 port rows/columns" << endl;
		exit(1);
	}

	// Fallback for sectors with no attached UE: beam nearest boresight
	// (zenith 90 deg; azimuth handled by row_beam_boresight_a).
	int fb_z = 0;
	for (int zi = 1; zi < nZ; zi++)
		if (fabs(bs_tilt_zenith_angle_LCS[zi] - pi / REAL(2.0)) <
		    fabs(bs_tilt_zenith_angle_LCS[fb_z] - pi / REAL(2.0)))
			fb_z = zi;

	for (int flat = 0; flat < num_SECTORS; flat++)
	{
		const int bs_idx  = flat / 3;
		const int sec_idx = flat % 3;

		Real count_z[8] = {0,};
		Real count_a[8] = {0,};
		int  num_attached = 0;

		for (int u = 0; u < num_MS; u++)
		{
			if (links[u]._sector_in_control != flat) continue;
			num_attached++;

			const beam_selection& sel = links[u].analog_beam_selection[flat];

			Real rsrp_z[8];
			Real best_z = -REAL(1e30);
			for (int zi = 0; zi < nZ; zi++)
			{
				rsrp_z[zi] = linear2dB(channel[bs_idx][u].signal_RSRP_gain
					[sec_idx][zi][sel.sector_a][sel.z][sel.a][sel.p]);
				if (rsrp_z[zi] > best_z) best_z = rsrp_z[zi];
			}
			RowBeam_Vote(rsrp_z, nZ, best_z, count_z);

			if (row_beam_az_mode == 2)
			{
				Real rsrp_a[8];
				Real best_a = -REAL(1e30);
				for (int ai = 0; ai < nA; ai++)
				{
					rsrp_a[ai] = linear2dB(channel[bs_idx][u].signal_RSRP_gain
						[sec_idx][sel.sector_z][ai][sel.z][sel.a][sel.p]);
					if (rsrp_a[ai] > best_a) best_a = rsrp_a[ai];
				}
				RowBeam_Vote(rsrp_a, nA, best_a, count_a);
			}
		}

		// --- apportion V rows over zenith beams ---
		int alloc_z[8] = {0,};
		if (num_attached == 0)
		{
			alloc_z[fb_z] = V;
		}
		else if (row_beam_force_uniform)
		{
			int am = 0;
			for (int zi = 1; zi < nZ; zi++) if (count_z[zi] > count_z[am]) am = zi;
			alloc_z[am] = V;
		}
		else
		{
			RowBeam_Apportion(count_z, nZ, V, alloc_z);
		}

		int check_z = 0;
		for (int zi = 0; zi < nZ; zi++) check_z += alloc_z[zi];
		if (check_z != V)
		{
			cout << "ERROR: Assign_Row_Beams zenith allocation " << check_z
			     << " != BS_Mp " << V << " (sector " << flat << ")" << endl;
			exit(1);
		}

		// --- place: steepest downtilt (largest zenith) on the top row, downwards.
		// Element z-coordinate grows with mi (d_tx.z = m*dV), so the top row is
		// mi = V-1; fill from there. Angle tie → lower beam index first.
		{
			bool placed[8] = {false,};
			int mi = V - 1;
			for (int rank = 0; rank < nZ; rank++)
			{
				int pick = -1;
				for (int zi = 0; zi < nZ; zi++)
				{
					if (placed[zi]) continue;
					if (pick < 0 || bs_tilt_zenith_angle_LCS[zi] > bs_tilt_zenith_angle_LCS[pick]) pick = zi;
				}
				placed[pick] = true;
				for (int r = 0; r < alloc_z[pick]; r++) sector[flat].row_beam_z[mi--] = pick;
			}
		}

		// --- apportion H columns over azimuth beams (az_mode==2 only) ---
		if (row_beam_az_mode == 2)
		{
			int alloc_a[8] = {0,};
			if (num_attached == 0)
			{
				alloc_a[row_beam_boresight_a] = H;
			}
			else if (row_beam_force_uniform)
			{
				int am = 0;
				for (int ai = 1; ai < nA; ai++) if (count_a[ai] > count_a[am]) am = ai;
				alloc_a[am] = H;
			}
			else
			{
				RowBeam_Apportion(count_a, nA, H, alloc_a);
			}

			int check_a = 0;
			for (int ai = 0; ai < nA; ai++) check_a += alloc_a[ai];
			if (check_a != H)
			{
				cout << "ERROR: Assign_Row_Beams azimuth allocation " << check_a
				     << " != BS_Np " << H << " (sector " << flat << ")" << endl;
				exit(1);
			}

			// place: azimuth ascending from column ni = 0 upwards
			bool placed[8] = {false,};
			int ni = 0;
			for (int rank = 0; rank < nA; rank++)
			{
				int pick = -1;
				for (int ai = 0; ai < nA; ai++)
				{
					if (placed[ai]) continue;
					if (pick < 0 || bs_tilt_azimuth_angle_LCS[ai] < bs_tilt_azimuth_angle_LCS[pick]) pick = ai;
				}
				placed[pick] = true;
				for (int r = 0; r < alloc_a[pick]; r++) sector[flat].col_beam_a[ni++] = pick;
			}
		}

		cout << "  [RowBeam] sector " << flat << " ue=" << num_attached << " alloc_z=[";
		for (int zi = 0; zi < nZ; zi++) cout << (zi ? "," : "") << alloc_z[zi];
		cout << "] rows(top->bottom)=[";
		for (int mi = V - 1; mi >= 0; mi--) cout << (mi < V - 1 ? "," : "") << sector[flat].row_beam_z[mi];
		cout << "]";
		if (row_beam_az_mode == 2)
		{
			cout << " cols=[";
			for (int ni = 0; ni < H; ni++) cout << (ni ? "," : "") << sector[flat].col_beam_a[ni];
			cout << "]";
		}
		cout << endl;
	}

	// --- refresh serving-link bookkeeping with the row-mixed RSRP ---
	// Per-port RSRP gain depends only on that port's beam, so the array-wide
	// value is exactly the arithmetic mean of the assigned beams' cached gains.
	for (int u = 0; u < num_MS; u++)
	{
		const int flat    = links[u]._sector_in_control;
		const int bs_idx  = flat / 3;
		const int sec_idx = flat % 3;
		const beam_selection& sel = links[u].analog_beam_selection[flat];

		Real g = 0;
		for (int mi = 0; mi < V; mi++)
		{
			if (row_beam_az_mode == 2)
			{
				for (int ni = 0; ni < H; ni++)
					g += channel[bs_idx][u].signal_RSRP_gain
						[sec_idx][sector[flat].row_beam_z[mi]][sector[flat].col_beam_a[ni]][sel.z][sel.a][sel.p];
			}
			else
			{
				const int a_used = (row_beam_az_mode == 1) ? row_beam_boresight_a : sel.sector_a;
				g += channel[bs_idx][u].signal_RSRP_gain
					[sec_idx][sector[flat].row_beam_z[mi]][a_used][sel.z][sel.a][sel.p];
			}
		}
		g /= (Real)((row_beam_az_mode == 2) ? V * H : V);

		const Real rsrp_mix_dB = linear2dB(g);
		links[u].link_RSRP    = rsrp_mix_dB;
		links[u].link_antgain = rsrp_mix_dB;
		links[u].str_signal   = bs_maxpower + rsrp_mix_dB - links[u].link_pathloss;
		links[u].static_gain[flat].first = links[u].str_signal;
	}
}

void Link_configuration()
{
	////////////////////////////////////////////////////////////////////////////////////////////////////
	for (int idx = 0; idx < num_SECTORS; idx++)
	{
		sector[idx].ue_in_control.clear();
		sector[idx].ue_in_comp.clear();
		sector[idx].W = NULL;
		sector[idx].PMI_read = NULL;
		sector[idx].PMI_vector_read = NULL;
	}

	for (int ms_idx = 0; ms_idx < num_MS; ms_idx++)
	{
		// 220808 jhnoh 
		if (drop_idx == 0) {

			links[ms_idx].LS_gain               = new Real          [num_SECTORS];
			links[ms_idx].adj_sector            = new int             [num_SECTORS];
			links[ms_idx].static_gain           = new Real_int_pair [num_SECTORS];
			links[ms_idx].analog_beam_selection = new beam_selection  [num_SECTORS];
			links[ms_idx].intf_w_rnd_RSRP	    = new Real          [num_SECTORS];

			links[ms_idx].rand_sec_a            = new int             [num_SECTORS];
			links[ms_idx].rand_sec_z            = new int             [num_SECTORS];
		}
		else
		{
			for (int sec_idx = 0; sec_idx < num_SECTORS; sec_idx++ )
			{
				links[ms_idx].LS_gain              [sec_idx] = 0;
				links[ms_idx].adj_sector           [sec_idx] = 0;

				links[ms_idx].static_gain          [sec_idx].first  = 0;
				links[ms_idx].static_gain          [sec_idx].second = 0;

				links[ms_idx].analog_beam_selection[sec_idx].a = 0;
				links[ms_idx].analog_beam_selection[sec_idx].z = 0;
				links[ms_idx].analog_beam_selection[sec_idx].p = 0;
				links[ms_idx].analog_beam_selection[sec_idx].sector_a = 0;
				links[ms_idx].analog_beam_selection[sec_idx].sector_z = 0;

				links[ms_idx].intf_w_rnd_RSRP	   [sec_idx] = 0;
				links[ms_idx].rand_sec_a           [sec_idx] = 0;
				links[ms_idx].rand_sec_z           [sec_idx] = 0;
			}
		}
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////
	Set_antenna_location_vector();
	Generate_bs_2D_DFT_beam_precoder();

	///////////////////////////////////////////////////////////////////////////////////////////////////
	cout << "num_threads = " << num_of_threads << endl;
	cout << "Get RSRP.........................."<<std::flush;
	#if ENABLE_MULTITHREADING
	#pragma omp parallel num_threads(num_of_threads)
	{
		#pragma omp for
	#endif
		for (int ms_idx = 0; ms_idx < num_MS; ms_idx++)
		{
			//cout << '\r'<< endl;
			//cout << "Get RSRP : " << ms_idx + 1 << "/" << num_MS << endl;
			links[ms_idx].Configuration(ms_idx);
			sector[links[ms_idx]._sector_in_control].ue_in_control.push_back(ms_idx);
			//sector[links[ms_idx].comp_sector_idx].ue_in_comp.push_back(ms_idx);
		}
	#if ENABLE_MULTITHREADING
	}
	#endif
	cout << "DONE"<<endl;

	Assign_Row_Beams();  // no-op unless row_beam_enable

	cout << "delete link memory................"<<std::flush;
	for (int ms_idx = 0; ms_idx < num_MS; ms_idx++)
	{
		links[ms_idx].Delete_link_memory();
	}
	cout << "DONE" << endl;
}

void InitializeSystem()
{
	CQI_Map_creation_Map_connection();
	Schedule_Map_creation_Connection();
	Type1_Codebook_Gen();
	Type2_Codebook_Gen();
	//Codebook_Generation();
}




/*===================================================================
FUNCTION: CQI_Map_creation_Map_connection()
===================================================================*/
void CQI_Map_creation_Map_connection()
{
	// CQI Map creation
	//ppppPMI_map        = new int        ***[num_MS];
	ppppPMI_map        = new PMI_FEEDBACK ***[num_MS];
	ppppCQI_Map        = new Real       ***[num_MS];
	ppppCQI_comp_Map   = new Real       ***[num_MS];
	ppppPMI_vector_map = new MatrixXcReal    ***[num_MS];
	ppppCSI_matrix_map = new MatrixXcReal    ***[num_MS];  // TDD: Full channel matrices

	for (int ue_idx = 0; ue_idx < num_MS; ue_idx++)
	{
		//ppppPMI_map[ue_idx]        = new int              ** [num_compute_coef];
		ppppPMI_map[ue_idx]        = new PMI_FEEDBACK     ** [num_compute_coef];
		ppppCQI_Map[ue_idx]        = new Real           ** [num_compute_coef];
		ppppCQI_comp_Map[ue_idx]   = new Real           ** [num_compute_coef];
		ppppPMI_vector_map[ue_idx] = new MatrixXcReal        ** [num_compute_coef];
		ppppCSI_matrix_map[ue_idx] = new MatrixXcReal        ** [num_compute_coef];

		for (int coeff_idx = 0; coeff_idx < num_compute_coef; coeff_idx++ )
		{
			//ppppPMI_map       [ue_idx][coeff_idx]        = new int              * [num_rb];
			ppppPMI_map       [ue_idx][coeff_idx]        = new PMI_FEEDBACK     * [num_rb];
			ppppCQI_Map       [ue_idx][coeff_idx]        = new Real           * [num_rb];
			ppppCQI_comp_Map  [ue_idx][coeff_idx]        = new Real           * [num_rb];
			ppppPMI_vector_map[ue_idx][coeff_idx]        = new MatrixXcReal        * [num_rb];
			ppppCSI_matrix_map[ue_idx][coeff_idx]        = new MatrixXcReal        * [num_rb];
			for (int freq_idx = 0; freq_idx < num_rb; freq_idx++)
			{
				//ppppPMI_map       [ue_idx][coeff_idx][freq_idx]        = new       int [cqi_history_length];
				ppppPMI_map       [ue_idx][coeff_idx][freq_idx]        = new PMI_FEEDBACK [cqi_history_length];
				ppppCQI_Map       [ue_idx][coeff_idx][freq_idx]        = new       Real [cqi_history_length];
				ppppCQI_comp_Map  [ue_idx][coeff_idx][freq_idx]        = new       Real [cqi_history_length];
				ppppPMI_vector_map[ue_idx][coeff_idx][freq_idx]        = new    MatrixXcReal [cqi_history_length];
				ppppCSI_matrix_map[ue_idx][coeff_idx][freq_idx]        = new    MatrixXcReal [cqi_history_length];

				// Initialize CSI matrices with proper size (NUM_RX_Port × NUM_TX_Port)
				for (int hist_idx = 0; hist_idx < cqi_history_length; hist_idx++)
				{
					ppppCSI_matrix_map[ue_idx][coeff_idx][freq_idx][hist_idx] = MatrixXcReal::Zero(NUM_RX_Port, NUM_TX_Port);
				}
			}
		}
	}
	// CQI Map to MS connection

	for (int ue_idx = 0; ue_idx < num_MS; ue_idx++)
	{
		ms[ue_idx].PMI = ppppPMI_map[ue_idx];
		ms[ue_idx].CQI = ppppCQI_Map[ue_idx];
		ms[ue_idx].CQI_comp = ppppCQI_comp_Map[ue_idx];
		ms[ue_idx].PMI_vector = ppppPMI_vector_map[ue_idx];
		ms[ue_idx].CSI_matrix = ppppCSI_matrix_map[ue_idx];  // TDD: Connect CSI matrix
	}
}


/*===================================================================
FUNCTION: Schedule_Map_creation_Connection()
===================================================================*/
void Schedule_Map_creation_Connection()
{

	// Schedule Map Creation

	pppSchedule_Map = new SCHEDULE_DECISION **[num_SECTORS + num_mTRP_SECTORS];
	for (int sector_idx = 0; sector_idx < num_SECTORS + num_mTRP_SECTORS; sector_idx++)
	{
		pppSchedule_Map[sector_idx] = new SCHEDULE_DECISION *[num_rb];

		for (int rb_idx = 0; rb_idx < num_rb; rb_idx++)
		{
			pppSchedule_Map[sector_idx][rb_idx] = new SCHEDULE_DECISION[mx_ue_mumimo];
		}
	}

	// Schedule Map to Sector Connection
	for (int sector_idx = 0; sector_idx < num_SECTORS + num_mTRP_SECTORS; sector_idx++)
	{
		sector[sector_idx].ppschedulewrite = pppSchedule_Map[sector_idx];
	}
}


void Type1_Codebook_Gen()
{
	int N_1 = BS_Np;
	int N_2 = BS_Mp;
	
	int P_CSI_RS = 2 * N_1 * N_2;
	
	int num_of_layer = 1;
	const ComplexReal j(0.0, 1.0);

	int O_1 = 4;
	int O_2 = 4;

	RowVectorXcReal u_m((int)N_2);
	MatrixXcReal    v_lm((int)N_1*N_2,1);

	if ( P_CSI_RS == 2 )
	{
		MatrixXcReal zero_mat = MatrixXcReal::Zero(2, 1);

		cb_W_csirs_2_layer_1[0] = zero_mat;
		cb_W_csirs_2_layer_1[1] = zero_mat;
		cb_W_csirs_2_layer_1[2] = zero_mat;
		cb_W_csirs_2_layer_1[3] = zero_mat;

		cb_W_csirs_2_layer_1[0] << 1, 1;
		cb_W_csirs_2_layer_1[1] << 1, j;
		cb_W_csirs_2_layer_1[2] << 1,-1;
		cb_W_csirs_2_layer_1[3] << 1,-j;

		cb_W_csirs_2_layer_1[0] = 1/sqrt(2) * cb_W_csirs_2_layer_1[0]; 
		cb_W_csirs_2_layer_1[1] = 1/sqrt(2) * cb_W_csirs_2_layer_1[1];
		cb_W_csirs_2_layer_1[2] = 1/sqrt(2) * cb_W_csirs_2_layer_1[2];
		cb_W_csirs_2_layer_1[3] = 1/sqrt(2) * cb_W_csirs_2_layer_1[3];
	}
	else
	{
		//cout << "[DEBUG-CB] P_CSI_RS != 2, allocating codebook_W..." << endl;
		//cout << "[DEBUG-CB] Codebook dimensions: [" << N_1*O_1 << "][" << N_2*O_2 << "][4]" << endl;
		//cout << "[DEBUG-CB] Each codebook entry size: " << NUM_TX_Port << " x 1" << endl;
		
		if (codebook_W == NULL) 
		{
			//cout << "[DEBUG-CB] Allocating new codebook_W array..." << endl;
			codebook_W = new MatrixXcReal ** [N_1*O_1];
			for (int l = 0; l < N_1*O_1; l++)
			{
				if (l % 32 == 0) {
					//cout << "[DEBUG-CB] Allocating codebook_W[" << l << "]..." << endl;
				}
				codebook_W[l] = new MatrixXcReal * [N_2*O_2];
				for (int m = 0; m < N_2*O_2; m++)
				{
					codebook_W[l][m] = new MatrixXcReal [4];
					for (int n = 0; n < 4; n++)
					{
						MatrixXcReal temp = MatrixXcReal::Zero(NUM_TX_Port, 1);
						codebook_W[l][m][n] = temp;
					}
				}
			}
		}

		//cout << "[DEBUG-CB] Starting codebook value calculation..." << endl;
		for (int i_11 = 0; i_11 < N_1*O_1; i_11++)
		{
			//if (i_11 % 32 == 0) {
				//cout << "[DEBUG-CB] Processing i_11=" << i_11 << "/" << (N_1*O_1) << endl;
			//}
			int l = i_11;
			for (int i_12 = 0; i_12 < N_2*O_2; i_12++)
			{	
				int m = i_12;
				for (int i = 0; i < N_2; i++)
				{
					u_m(i) = exp(j*(Real)2*pi*(Real)m*(Real)i / (Real)(O_2*N_2) );
				}

				for(int x = 0; x < N_1; x++)
				{
					for(int y = 0; y < N_2; y++)
					{
						v_lm(N_2*x + y) = u_m(y) * exp( (j*(Real)2*pi*(Real)l*(Real)x / (Real)(O_1*N_1)) );
					}
				}

				for ( int i_2 = 0; i_2 < 4; i_2++)
				{
					int n = i_2;
					ComplexReal phase = exp(j*pi*(Real)n / (Real)2);
					if (BS_P == 2) {
						// Dual-pol: W = [v_lm; phase * v_lm]
						codebook_W[l][m][n] << v_lm , phase * v_lm;
					} else {
						// Single-pol (BS_P=1): W = v_lm (no co-phasing)
						codebook_W[l][m][n] = v_lm;
					}
					codebook_W[l][m][n] = (1 / sqrt(NUM_TX_Port)) * codebook_W[l][m][n];
				}
			}
		}
		//cout << "[DEBUG-CB] Codebook value calculation completed" << endl;
	}
	//cout << "[DEBUG-CB] Type1_Codebook_Gen completed successfully" << endl;
}

/*===================================================================
FUNCTION: Type2_Codebook_Gen()

DESCRIPTION:
  Generates and caches the oversampled 2D DFT beam basis used by the
  Type II codebook (TS 38.214 §5.2.2.2.3). The full Type II precoder
  is assembled on-the-fly per UE / per subband during PMI selection
  (see MS::Quantization_of_Ch_Type2), so this routine only builds the
  shared beam dictionary:

        v_(l,m) = u_m ⊗ u_l,
        u_l[x] = exp(j 2π l x / (O1 N1)),   x = 0..N1-1
        u_m[y] = exp(j 2π m y / (O2 N2)),   y = 0..N2-1

  Beam layout expected by channel_update.cpp:
    For a spatial-domain tx antenna element (x, y),
    spatial index = N2 * x + y   (matches Type I convention).

  Skipped when Type II is disabled (g_codebook_type != 2) or when
  P_CSI_RS == 2 (uses the 2-port hardcoded book instead).
===================================================================*/
void Type2_Codebook_Gen()
{
	// BUG FIX 2026-05-25: eType II (g_codebook_type == 3) also uses the same
	// oversampled DFT beam basis. Previously only Rel-15 Type II (=2) generated
	// type2_beam_v, causing eType II runs to silently fall back to Quantization_of_Ch
	// (Type I, rank 1). This wiped out the entire eType II quantization path.
	if (g_codebook_type != 2 && g_codebook_type != 3) return;

	const int N_1 = BS_Np;
	const int N_2 = BS_Mp;
	const int P_CSI_RS = 2 * N_1 * N_2;
	if (P_CSI_RS == 2) {
		cout << "[Type II] P_CSI_RS=2, skipping (2-port book used instead)" << endl;
		return;
	}

	const int O_1 = 4;
	const int O_2 = (N_2 == 1) ? 1 : 4;

	type2_N1 = N_1;
	type2_N2 = N_2;
	type2_O1 = O_1;
	type2_O2 = O_2;

	const ComplexReal j(0.0, 1.0);

	// Free previous allocation (drop-boundary safety)
	if (type2_beam_v != NULL) {
		for (int l = 0; l < O_1 * N_1; l++) {
			delete [] type2_beam_v[l];
		}
		delete [] type2_beam_v;
		type2_beam_v = NULL;
	}

	type2_beam_v = new VectorXcReal * [O_1 * N_1];
	for (int l = 0; l < O_1 * N_1; l++) {
		type2_beam_v[l] = new VectorXcReal[O_2 * N_2];
		for (int m = 0; m < O_2 * N_2; m++) {
			VectorXcReal v_lm = VectorXcReal::Zero(N_1 * N_2);
			for (int x = 0; x < N_1; x++) {
				ComplexReal ul_x = exp(j * (Real)2.0 * pi * (Real)l * (Real)x / (Real)(O_1 * N_1));
				for (int y = 0; y < N_2; y++) {
					ComplexReal um_y = exp(j * (Real)2.0 * pi * (Real)m * (Real)y / (Real)(O_2 * N_2));
					v_lm(N_2 * x + y) = ul_x * um_y;
				}
			}
			type2_beam_v[l][m] = v_lm;  // length N1*N2, not yet normalized
		}
	}

	cout << "[Type II] Codebook beams cached: (N1,N2)=(" << N_1 << "," << N_2
	     << "), (O1,O2)=(" << O_1 << "," << O_2
	     << "), total DFT beams = " << (O_1 * N_1 * O_2 * N_2)
	     << ", L=" << g_type2_L
	     << ", phase=" << g_type2_phase_alphabet
	     << ", SA=" << g_type2_subband_amplitude << endl;
}

void ClearMAPs2StartADrop()
{
	for (int ue_idx = 0; ue_idx < num_MS; ue_idx++)
	{
		for (int coeff_idx = 0; coeff_idx < num_compute_coef; coeff_idx++ )
		{
			for (int freq_idx = 0; freq_idx < num_rb; freq_idx++)
			{
				for (int time_idx = 0; time_idx < cqi_history_length; time_idx++)
				{
					ppppCQI_Map[ue_idx][coeff_idx][freq_idx][time_idx]        = 0.;
					ppppCQI_comp_Map[ue_idx][coeff_idx][freq_idx][time_idx]   = 0.;
					ppppPMI_map[ue_idx][coeff_idx][freq_idx][time_idx]        = {0,0,0};
					ppppPMI_vector_map[ue_idx][coeff_idx][freq_idx][time_idx] = MatrixXcReal::Zero(NUM_TX_Port, g_type2_rank);
				}
			}
		}
	}

	for (int ue_idx = 0; ue_idx < num_MS; ue_idx++)
	{
		ms[ue_idx].CQI        = ppppCQI_Map[ue_idx];
		ms[ue_idx].CQI_comp   = ppppCQI_comp_Map[ue_idx];
		ms[ue_idx].PMI        = ppppPMI_map[ue_idx];
		ms[ue_idx].PMI_vector = ppppPMI_vector_map[ue_idx];
	}	

	SCHEDULE_DECISION no_decision;

	no_decision.ue_selected = NO_UE;
	no_decision.mcs_selected = NO_CQI;
	no_decision.cqi_selected = NO_CQI;
	no_decision.capacity = -1;

	for (int sector_idx = 0; sector_idx < num_SECTORS; sector_idx++)
	{
		for (int freq_idx = 0; freq_idx < num_rb; freq_idx++)
		{
			for (int stream = 0; stream < mx_ue_mumimo; stream++)
			{
				pppSchedule_Map[sector_idx][freq_idx][stream] = no_decision;
			}
		}
	}

	for (int ue_idx = 0; ue_idx < num_MS; ue_idx++)
	{
		ms[ue_idx].ppSchedulerRead = pppSchedule_Map[links[ue_idx]._sector_in_control];
	}

	for (int sector_idx = 0; sector_idx < num_SECTORS; sector_idx++)
	{
		sector[sector_idx].ppschedulewrite = pppSchedule_Map[sector_idx];
	}
}

void Get_CouplingLoss()
{
	for (int ue_idx = 0; ue_idx < num_MS; ue_idx++)
	{
		if((tilt_azimuth_angle_LCS_size == 1) && (tilt_zenith_angle_LCS_size == 1))   /// Indoor Config A, Rural All scenario
		{
			Real linear_signal       = dBm2linear(links[ue_idx].str_signal); //pow(10., (links[ue_idx].str_signal) / 10.);
			Real linear_interference = dBm2linear(links[ue_idx].interference); //pow(10., links[ue_idx].interference / 10.);

			Real sir = linear2dB(linear_signal / linear_interference);  // Wideband SIR (without noise)
			Get_Wideband_SIR(sir, ue_idx);

			linear_interference        = (linear_interference + noise);
			Real geometry            = linear2dB(linear_signal / linear_interference);  // Geometry (SINR with noise)

			Get_Geometry(geometry, ue_idx);
		}
		else if (TYPE == 11)
		{
			//if ((Configuration_Type == 2) && (num_Indoor_TRxP == 1))  // indoor, Config C, 12TRxP
			if (true)  
			{
				Real linear_signal = pow(10., (links[ue_idx].str_signal - bs_maxpower) / 10.);
				Real new_interference = 0;

				int sector_parameter;

				if (TYPE == 11) // indoor
				{
					sector_parameter = num_Indoor_TRxP; // 1 or 3
				}
				else
				{
					sector_parameter = 3;
				}

				Real new_bspower = dBm2linear(bs_maxpower);

				for (int bs_idx = 0; bs_idx < num_BS; bs_idx++)
				{
					for (int sec_idx = 0; sec_idx < sector_parameter; sec_idx++)  ///
					{
						if ((bs_idx == links[ue_idx].self_bs_idx) && (sec_idx == (links[ue_idx].self_sector_idx)))
						{
							// no interference
						}
						else
						{
							int aa = links[ue_idx].rand_sec_a[3*bs_idx+sec_idx]; //(int)((Real)tilt_azimuth_angle_LCS_size * randnum.u());
							int zz = links[ue_idx].rand_sec_z[3*bs_idx+sec_idx]; //(int)((Real)tilt_zenith_angle_LCS_size * randnum.u());
							Real RSRP;

							if (ue_antenna_element_gain == 0)
							{
								RSRP = linear2dB(channel[bs_idx][ue_idx].signal_RSRP_gain[sec_idx][zz][aa][0][0][0]);
							}
							else
							{
								RSRP = linear2dB(channel[bs_idx][ue_idx].signal_RSRP_gain[sec_idx][zz][aa]
								[links[ue_idx].zenith_angle_idx_selected_for_interference]
								[links[ue_idx].azimuth_angle_idx_selected_for_interference]
								[links[ue_idx].panel_idx_selected_for_interference]);
							}

							Real interf = RSRP - channel[bs_idx][ue_idx].pathloss_final;								
							new_interference += pow(10., (interf) / 10.);
						}
					}
				}
				Real sir = linear2dB(linear_signal / new_interference);  // Wideband SIR (without noise)
				Get_Wideband_SIR(sir, ue_idx);
				Real geometry = linear2dB(linear_signal*new_bspower / (new_interference*new_bspower + noise));
				Get_Geometry(geometry, ue_idx);
			}
			else
			{
				Real linear_signal = pow(10., (links[ue_idx].str_signal) / 10.);
				Real linear_interference = pow(10., links[ue_idx].interference / 10.);

				Real sir = linear2dB(linear_signal / linear_interference);  // Wideband SIR (without noise)
				Get_Wideband_SIR(sir, ue_idx);

				linear_interference = (linear_interference + noise);
				Real geometry = linear2dB(linear_signal / linear_interference);

				Get_Geometry(geometry, ue_idx);
			}
		}
		else  // dense urban
		{
			Real linear_signal = pow(10., (links[ue_idx].str_signal - bs_maxpower) / 10.);
			Real new_interference = 0;
			int sector_parameter;

			if (TYPE == 11) // indoor
			{
				sector_parameter = num_Indoor_TRxP; // 1 or 3
			}
			else
			{
				sector_parameter = 3;
			}
			int ms_to_bs_wrap_idx = (ue_idx / (3*num_MS_persector)); // -> ms_idx�� bs_idx��
			Real new_bspower = dBm2linear(bs_maxpower);

			//for (int bs_idx = 0; bs_idx < num_BS; bs_idx++)
			for (int wrap_bs_idx = 0; wrap_bs_idx < simple_num_BS; wrap_bs_idx++)
			{
				int bs_idx = simple_wrap_mat[ms_to_bs_wrap_idx][wrap_bs_idx];
				for (int sec_idx = 0; sec_idx < sector_parameter; sec_idx++)  ///
				{
					if ((bs_idx == links[ue_idx].self_bs_idx) && (sec_idx == (links[ue_idx].self_sector_idx)))
					{
						// no interference
					}
					else
					{
						int aa = links[ue_idx].rand_sec_a[3*bs_idx+sec_idx]; //int aa = (int)((Real)tilt_azimuth_angle_LCS_size * randnum.u());
						int zz = links[ue_idx].rand_sec_z[3*bs_idx+sec_idx]; //int zz = (int)((Real)tilt_zenith_angle_LCS_size * randnum.u());
						Real RSRP;

						if ( !(-1 < aa && aa < tilt_azimuth_angle_LCS_size) )
							cout << "Something wrong with aziumth" << endl;

						if ( !(-1 < zz && zz < tilt_zenith_angle_LCS_size) )
							cout << "Something wrong with aziumth" << endl;

						if (row_beam_enable)
						{
							// Row-beam mode: interference from this sector's ACTUAL
							// assigned beams — average the cached per-beam gains over its
							// port rows (and columns for az_mode==2) instead of one random
							// zenith beam. zz/aa keep being drawn in compute_interference
							// so the RNG stream is identical whether the flag is on or off.
							const int flat_i = 3 * bs_idx + sec_idx;
							int uz = 0, ua = 0, up = 0;
							if (ue_antenna_element_gain != 0)
							{
								uz = links[ue_idx].zenith_angle_idx_selected_for_interference;
								ua = links[ue_idx].azimuth_angle_idx_selected_for_interference;
								up = links[ue_idx].panel_idx_selected_for_interference;
							}
							Real gain = 0;
							for (int mi = 0; mi < BS_Mp; mi++)
							{
								if (row_beam_az_mode == 2)
								{
									for (int ni = 0; ni < BS_Np; ni++)
										gain += channel[bs_idx][ue_idx].signal_RSRP_gain[sec_idx]
											[sector[flat_i].row_beam_z[mi]][sector[flat_i].col_beam_a[ni]][uz][ua][up];
								}
								else
								{
									// az_mode 0: the interferer serves its own UEs on their
									// azimuths — a random azimuth draw models that (keep aa).
									// az_mode 1: every sector transmits at boresight.
									const int a_used = (row_beam_az_mode == 1) ? row_beam_boresight_a : aa;
									gain += channel[bs_idx][ue_idx].signal_RSRP_gain[sec_idx]
										[sector[flat_i].row_beam_z[mi]][a_used][uz][ua][up];
								}
							}
							gain /= (Real)((row_beam_az_mode == 2) ? BS_Mp * BS_Np : BS_Mp);
							RSRP = linear2dB(gain);
						}
						else if (ue_antenna_element_gain == 0)
						{
							RSRP = linear2dB(channel[bs_idx][ue_idx].signal_RSRP_gain[sec_idx][zz][aa][0][0][0]);
						}
						else
						{
							RSRP = linear2dB(channel[bs_idx][ue_idx].signal_RSRP_gain[sec_idx][zz][aa]
							[links[ue_idx].zenith_angle_idx_selected_for_interference]
							[links[ue_idx].azimuth_angle_idx_selected_for_interference]
							[links[ue_idx].panel_idx_selected_for_interference]);
						}
						Real interf = RSRP - channel[bs_idx][ue_idx].pathloss_final;
						Real _intf_w_rnd_RSRP = bs_maxpower + RSRP - channel[bs_idx][ue_idx].pathloss_final;

						//if ( _intf_w_rnd_RSRP != links[ue_idx].intf_w_rnd_RSRP[3*bs_idx+sec_idx])
						//	cout << "Someting wrong with rnd RSRP" << endl;

						new_interference += pow(10., (interf) / 10.);
					}
				}
			}

			Real new_interference_2 = linear2dBm(new_interference * new_bspower);

			//if( new_interference_2 != links[ue_idx].interference)
			//	cout << "Interference computation error" << endl;

			links[ue_idx].interference = new_interference_2;
			//linear2dBm(new_interference * new_bspower);

			Real sir = linear2dB(linear_signal / new_interference);  // Wideband SIR (without noise)
			Get_Wideband_SIR(sir, ue_idx);
			Real geometry = linear2dB(linear_signal*new_bspower / (new_interference*new_bspower + noise));
			Get_Geometry(geometry, ue_idx);
		}
		Real cp = links[ue_idx].str_signal - bs_maxpower;  // str_signal = bs_maxpower + max_antgain + combined antgain - pathloss
		Coupling_Loss_Get(cp, ue_idx);

		//Print_Calib_Debug_Info( ue_idx );
	}
}

void Print_Calib_Debug_Info( void )
{
	if (drop_idx == 0)
	{
		char calib_debug_file_name[100];
		sprintf(calib_debug_file_name,"./%s/Calib_Debug_Info_drop.csv",folder_name);
		Calibration_Debug_info.open(calib_debug_file_name, ios::out);

		Calibration_Debug_info
		<< "drop_idx" <<","
		<< "ue_idx" <<","
		<< "sector_idx" <<","
		<< "bs_x" << ","
		<< "bs_y" << ","
		<< "x-coordinate"<<","
		<< "y-coordinate"<<","
		<< "height"<<","
		<< "distance"<<","
		<< "pathloss"<< ","
		<< "RSRP"<<","
		<< "coupling_loss"<<","
		<< "geometry"<< ","
		<< "indoor"<< ","
		<< "los"<< ","
		<< "high_loss_flag" << ","
		<< "LOS_AOD" << ","
		<< "LOS_AOA" << ","
		<< "LOS_ZOD" << ","
		<< "LOS_ZOA" << ","
		<< "ASD" << ","
		<< "ASA" << ","
		<< "ZSD" << ","
		<< "ZSA" << ","
		<< "SV1_max_dB" << ","
		<< "SV1_min_dB" << ","
		<< "SV1_avg_dB" << ","
		<< "SV2_max_dB" << ","
		<< "SV2_min_dB" << ","
		<< "SV2_avg_dB" << endl;
	}

	for(int ue_idx = 0; ue_idx < num_MS; ue_idx++)
	{
		int serving_bs = links[ue_idx].self_bs_idx;
		CHANNEL* ch_serving = &channel[serving_bs][ue_idx];

		Calibration_Debug_info
		<< drop_idx << ","
		<< ue_idx << ","
		<< links[ue_idx]._sector_in_control << ","
		<< bs[serving_bs].loc.x << ","
		<< bs[serving_bs].loc.y << ","
		<< ms[ue_idx].loc.x << ","
		<< ms[ue_idx].loc.y << ","
		<< ms[ue_idx].MS_HEIGHT_FINAL << ","
		<< links[ue_idx].link_distance<<","
		<< links[ue_idx].link_pathloss << ","
		<< links[ue_idx].link_RSRP << ","
		<< ms[ue_idx].coupling_loss << ","
		<< ms[ue_idx].geometry << ","
		<< links[ue_idx].link_indoor << ","
		<< links[ue_idx].link_los << ","
		<< links[ue_idx].high_loss_flag << ","
		<< ch_serving->LOS_AOD_GCS << ","
		<< ch_serving->LOS_AOA_GCS << ","
		<< ch_serving->LOS_ZOD_GCS << ","
		<< ch_serving->LOS_ZOA_GCS << ","
		<< ch_serving->ASD << ","
		<< ch_serving->ASA << ","
		<< ch_serving->ZSD << ","
		<< ch_serving->ZSA << ",";

		// Per-UE singular value statistics (max, min, avg across RBs)
		if (ms[ue_idx].H_m != NULL && num_rb > 0) {
			Real sv1_max = -1e30, sv1_min = 1e30, sv1_sum = 0;
			Real sv2_max = -1e30, sv2_min = 1e30, sv2_sum = 0;
			int sv_count = 0;

			for (int rb = 0; rb < num_rb; rb++) {
				MatrixXcReal H = ms[ue_idx].H_m[0][rb];
				Eigen::JacobiSVD<MatrixXcReal> svd(H, Eigen::ComputeThinU | Eigen::ComputeThinV);
				auto sv = svd.singularValues();

				if (sv.size() >= 1 && sv(0) > 1e-30) {
					Real sv1_dB = 20.0 * log10(sv(0));
					if (sv1_dB > sv1_max) sv1_max = sv1_dB;
					if (sv1_dB < sv1_min) sv1_min = sv1_dB;
					sv1_sum += sv1_dB;
				}
				if (sv.size() >= 2 && sv(1) > 1e-30) {
					Real sv2_dB = 20.0 * log10(sv(1));
					if (sv2_dB > sv2_max) sv2_max = sv2_dB;
					if (sv2_dB < sv2_min) sv2_min = sv2_dB;
					sv2_sum += sv2_dB;
				}
				sv_count++;
			}

			if (sv_count > 0) {
				Calibration_Debug_info
				<< sv1_max << "," << sv1_min << "," << sv1_sum / sv_count << ","
				<< sv2_max << "," << sv2_min << "," << sv2_sum / sv_count << endl;
			} else {
				Calibration_Debug_info << ",,,,,," << endl;
			}
		} else {
			Calibration_Debug_info << ",,,,,," << endl;
		}
	}

	if (drop_idx == num_drop - 1)
	{
		Calibration_Debug_info.close();
	}

	// Cluster parameter CSV dump
	if (g_dump_cluster_params)
	{
		static ofstream cluster_csv;

		if (drop_idx == 0)
		{
			char cluster_file_name[100];
			sprintf(cluster_file_name, "./%s/Cluster_Params_drop.csv", folder_name);
			cluster_csv.open(cluster_file_name, ios::out);
			cluster_csv << "drop_idx,ue_idx,sector_idx,indoor,los,num_clusters,cluster_idx,power_dB,delay_ns,AOD,AOA,ZOD,ZOA" << endl;
		}

		for (int ue_idx = 0; ue_idx < num_MS; ue_idx++)
		{
			int serving_bs = links[ue_idx].self_bs_idx;
			CHANNEL* ch = &channel[serving_bs][ue_idx];
			int idx1 = ch->strongest_power_idx;
			int idx2 = ch->strongest_power_idx2;

			for (int n = 0; n < ch->num_path; n++)
			{
				// Use powerForAngles which includes LOS K-factor for cluster 0
				// Restore pre-subcluster power: strongest 2 clusters were split 10/6/4
				// so current power is 10/20 of original → multiply by 20/10 = 2
				Real pwr = ch->powerForAngles[n];
				if (n == idx1 || n == idx2)
					pwr *= 2.0;

				cluster_csv
					<< drop_idx << ","
					<< ue_idx << ","
					<< links[ue_idx]._sector_in_control << ","
					<< links[ue_idx].link_indoor << ","
					<< links[ue_idx].link_los << ","
					<< ch->num_path << ","
					<< n << ","
					<< 10.0 * log10(pwr) << ","
					<< ch->delay[n] * 1e9 << ","
					<< ch->AOD[n] << ","
					<< ch->AOA[n] << ","
					<< ch->ZOD[n] << ","
					<< ch->ZOA[n] << endl;
			}
		}

		if (drop_idx == num_drop - 1)
		{
			cluster_csv.close();
		}
	}
}




Real Get_distance(LOCATION a, LOCATION b)
{
	return sqrt((a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y));
}




/*===================================================================
FUNCTION: loading()

Print loading percent
===================================================================*/
void loading(void)
{
	
	//system("clear");
	//system("cls");
	
	if (t < SCHEDULE_DELAY)
	{
		cout << " Initial Setting is going on " << abs(t - SCHEDULE_DELAY) << endl;
	}
	else
	{
		cout << " Pleas Wait Until it get to the 100% " << endl;
		cout << "------------------------" << "\t"
			<< 100.*(Real)(t - SCHEDULE_DELAY + 1 + drop_idx*(N_pf + run_times)) / (Real)((N_pf + run_times)*num_drop) << "%" << "\t"
			<< "------------------------" << endl;
		cout << " Current Drop is " << drop_idx + 1 << "/" << num_drop << endl;
		cout << "------------------------" << "\t"
			<< 100.*(Real)(t - SCHEDULE_DELAY + 1) / (Real)(N_pf + run_times) << "%" << "\t"
			<< "------------------------" << endl;
		cout << " Output File name is `" << file_name << "`" << endl;
	}
	

}

void Parameter_initialization()
{
	//// initialization per drop
	cout << "Parameter initialization..........";
	for (int idx = 0; idx < num_MS; idx++)
	{
		ms[idx].H_m = NULL;
		ms[idx].unscheduled_stack = 0;
	}
	sector_indices_schedule.clear();
	cout << "DONE" << endl;
}


void Set_antenna_location_vector()
{
	/////////////////////  array antenna rotate angle //////////////////////////////////////////////////
	for (int ms_idx = 0; ms_idx < num_MS; ms_idx++)
	{
		if (handheld_mode)
		{
			ms[ms_idx].alpha = 360 * randnum.u() * (pi / 180.);
			ms[ms_idx].beta  = handheld_beta_deg * (pi / 180.);
			ms[ms_idx].gamma = 0.;
		}
		else if (ue_antenna_element_gain == 0)
		{
			ms[ms_idx].alpha = 0.;
			ms[ms_idx].beta  = 0.;
			ms[ms_idx].gamma = 0.;
		}
		else
		{
			ms[ms_idx].alpha = 360 * randnum.u() * (pi / 180.);
			ms[ms_idx].beta  = 0.;
			ms[ms_idx].gamma = 0.;
		}
	}

	for (int bs_idx = 0; bs_idx < num_BS + num_mTRP; bs_idx++)
	{
		Real alpha;
		Real beta  = 0; //(Mechanic_tilt - 90.) * (pi / 180.);
		Real gamma = 0;

		if (TYPE == 11)  //InH
		{
			if (num_Indoor_TRxP == 1)
			{
				alpha = 0;

				for (int sector_idx = 0; sector_idx < 3; sector_idx++)
				{
					bs[bs_idx].ant[sector_idx][0] = alpha;/// [sector num][0=alpha or 1=beta or 2=gamma]
					bs[bs_idx].ant[sector_idx][1] = beta;
					
					bs[bs_idx].ant[sector_idx][2] = gamma;
				}
			}
			else if (num_Indoor_TRxP == 3)
			{
				for (int sector_idx = 0; sector_idx < 3; sector_idx++)
				{
					if (sector_idx == 0)
					{
						alpha = 30. * (pi / 180.);
					}
					else if (sector_idx == 1)
					{
						alpha = 150. * (pi / 180.);
					}
					else
					{
						alpha = 270. * (pi / 180.);
					}

					bs[bs_idx].ant[sector_idx][0] = alpha;/// [sector num][0=alpha or 1=beta or 2=gamma]
					bs[bs_idx].ant[sector_idx][1] = beta;
					bs[bs_idx].ant[sector_idx][2] = gamma;
				}
			}
		}
		else  // Outdoor scenario
		{
			for (int sector_idx = 0; sector_idx < 3; sector_idx++)
			{
				if (sector_idx == 0)
				{
					alpha = 30. * (pi / 180.);
				}
				else if (sector_idx == 1)
				{
					alpha = 150. * (pi / 180.);
				}
				else
				{
					alpha = 270. * (pi / 180.);
				}

				bs[bs_idx].ant[sector_idx][0] = alpha;/// [sector num][0=alpha or 1=beta or 2=gamma]
				bs[bs_idx].ant[sector_idx][1] = beta;
				bs[bs_idx].ant[sector_idx][2] = gamma;
			}
		}
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////
	// Handheld UT: 8 candidate antenna positions on 15cm×7cm rectangular device (TR 38.901 Table 7.8-2A)
	const Real HH_L = 0.15;  // device length (m)
	const Real HH_W = 0.07;  // device width (m)
	const LOCATION3D handheld_positions[8] = {
		{-HH_L/2, -HH_W/2, 0},  // position 1: bottom-left corner
		{      0, -HH_W/2, 0},  // position 2: bottom-center
		{ HH_L/2, -HH_W/2, 0},  // position 3: bottom-right corner
		{ HH_L/2,       0, 0},  // position 4: right-center
		{ HH_L/2,  HH_W/2, 0},  // position 5: top-right corner
		{      0,  HH_W/2, 0},  // position 6: top-center
		{-HH_L/2,  HH_W/2, 0},  // position 7: top-left corner
		{-HH_L/2,       0, 0},  // position 8: left-center
	};

	for (int ue_idx = 0; ue_idx < num_MS; ue_idx++)
	{
		if (handheld_mode)
		{
			// Handheld: place antenna elements at physical device positions
			// MS_M = handheld_num_ports, MS_N=1, MS_P=1
			for (int port = 0; port < handheld_num_ports; port++)
			{
				int ant_idx = handheld_port_indices[port] - 1;  // 1-based to 0-based
				LOCATION3D pos = handheld_positions[ant_idx];
				// Transform from device LCS to GCS using UE orientation
				ms[ue_idx].d_rx[port][0][0][0][0] = Transpose_LCS_to_GCS_location(
					ms[ue_idx].alpha, ms[ue_idx].beta, ms[ue_idx].gamma, pos);
			}
		}
		else
		{
		for (int mg = 0; mg < MS_Mg; mg++)
		{
			for (int ng = 0; ng < MS_Ng; ng++)
			{
				for (int m = 0; m<MS_M; m++)
				{
					for (int n = 0; n<MS_N; n++)
					{
						for (int p = 0; p < MS_P; p++)
						{
							ms[ue_idx].d_rx[m][n][p][mg][ng].x = 0;
							ms[ue_idx].d_rx[m][n][p][mg][ng].y = n * MS_dH + ng * MS_dgH;
							ms[ue_idx].d_rx[m][n][p][mg][ng].z = m * MS_dV + mg * MS_dgV;


							if (MS_Ng == 1)
							{
								ms[ue_idx].d_rx[m][n][p][mg][ng] = Transpose_LCS_to_GCS_location(ms[ue_idx].alpha, ms[ue_idx].beta, ms[ue_idx].gamma, ms[ue_idx].d_rx[m][n][p][mg][ng]);
							}
							else
							{
								if (ng == 0)
								{
									ms[ue_idx].d_rx[m][n][p][mg][ng] = Transpose_LCS_to_GCS_location(ms[ue_idx].alpha, ms[ue_idx].beta, ms[ue_idx].gamma, ms[ue_idx].d_rx[m][n][p][mg][ng]);
								}
								else
								{
									ms[ue_idx].d_rx[m][n][p][mg][ng] = Transpose_LCS_to_GCS_location(ms[ue_idx].alpha, ms[ue_idx].beta, ms[ue_idx].gamma, ms[ue_idx].d_rx[m][n][p][mg][ng]);
								}
							}
						}
					}
				}
			}
		}
		} // else (non-handheld)
	}

	for (int bs_idx = 0; bs_idx < num_BS + num_mTRP; bs_idx++)
	{
		for (int sec_idx = 0; sec_idx < 3; sec_idx++)
		{
			for (int mg = 0; mg < BS_Mg; mg++)
			{
				for (int ng = 0; ng < BS_Ng; ng++)
				{
					for (int m = 0; m<BS_M; m++)
					{
						for (int n = 0; n<BS_N; n++)
						{
							for (int p = 0; p<BS_P; p++)
							{
								bs[bs_idx].d_tx[sec_idx][m][n][p][mg][ng].x = 0;
								bs[bs_idx].d_tx[sec_idx][m][n][p][mg][ng].y = n * BS_dH + ng * BS_dgH;
								bs[bs_idx].d_tx[sec_idx][m][n][p][mg][ng].z = m * BS_dV + mg * BS_dgV;

								
								bs[bs_idx].d_tx[sec_idx][m][n][p][mg][ng] = Transpose_LCS_to_GCS_location(bs[bs_idx].ant[sec_idx][0], bs[bs_idx].ant[sec_idx][1], bs[bs_idx].ant[sec_idx][2], bs[bs_idx].d_tx[sec_idx][m][n][p][mg][ng]);
								//cout << bs[bs_idx].d_tx[sec_idx][m][n][p][mg][ng].x << " " << bs[bs_idx].d_tx[sec_idx][m][n][p][mg][ng].y << " " << bs[bs_idx].d_tx[sec_idx][m][n][p][mg][ng].z << endl;
							}
						}
					}
				}
			}
		}		
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}

LOCATION3D Transpose_LCS_to_GCS_location(Real alpha, Real beta, Real gamma, LOCATION3D X)
{
	LOCATION3D Y;

	Y.x = cos(alpha)*cos(beta) * X.x + (cos(alpha)*sin(beta)*sin(gamma) - sin(alpha)*cos(gamma)) * X.y + (cos(alpha)*sin(beta)*cos(gamma) + sin(alpha)*sin(gamma)) * X.z;
	Y.y = sin(alpha)*cos(beta) * X.x + (sin(alpha)*sin(beta)*sin(gamma) + cos(alpha)*cos(gamma)) * X.y + (sin(alpha)*sin(beta)*cos(gamma) - cos(alpha)*sin(gamma)) * X.z;
	Y.z = -sin(beta)          * X.x + cos(beta)*sin(gamma)                                      * X.y + cos(beta)*cos(gamma)                                      * X.z;

	return Y;
}

void Generate_bs_2D_DFT_beam_precoder()
{
	if (TYPE == 11)  // Indoor
	{
		if (Configuration_Type == 0)
		{
			bs_tilt_azimuth_angle_LCS[0] = 0.;
			bs_tilt_zenith_angle_LCS[0] = 90. * (pi / 180.);

			tilt_azimuth_angle_LCS_size = 1;
			tilt_zenith_angle_LCS_size = 1;
		}
		else if (Configuration_Type == 1)
		{
			/// for calibration
			//bs
			bs_tilt_azimuth_angle_LCS[0] = (-3 * pi / 8);
			bs_tilt_azimuth_angle_LCS[1] = (-1 * pi / 8);
			bs_tilt_azimuth_angle_LCS[2] = (pi / 8);
			bs_tilt_azimuth_angle_LCS[3] = (3 * pi / 8);


			bs_tilt_zenith_angle_LCS[0] = (pi / 4);
			bs_tilt_zenith_angle_LCS[1] = (3 * pi / 4);

			tilt_azimuth_angle_LCS_size = 4;
			tilt_zenith_angle_LCS_size = 2;

			//ue
			ue_tilt_azimuth_angle_LCS[0] = (-3 * pi / 8);
			ue_tilt_azimuth_angle_LCS[1] = (-pi / 8);
			ue_tilt_azimuth_angle_LCS[2] = (pi / 8);
			ue_tilt_azimuth_angle_LCS[3] = (3 * pi / 8);

			ue_tilt_zenith_angle_LCS[0] = (pi / 4);
			ue_tilt_zenith_angle_LCS[1] = (3 * pi / 4);

			ue_tilt_azimuth_angle_LCS_size = 4;
			ue_tilt_zenith_angle_LCS_size = 2;
		}
		else if (Configuration_Type == 2)
		{
			// for calibration
			
			//bs
			bs_tilt_azimuth_angle_LCS[0] = (-7 * pi / 16);
			bs_tilt_azimuth_angle_LCS[1] = (-5 * pi / 16);
			bs_tilt_azimuth_angle_LCS[2] = (-3 * pi / 16);
			bs_tilt_azimuth_angle_LCS[3] = (-1 * pi / 16);
			bs_tilt_azimuth_angle_LCS[4] = (1 * pi / 16);
			bs_tilt_azimuth_angle_LCS[5] = (3 * pi / 16);
			bs_tilt_azimuth_angle_LCS[6] = (5 * pi / 16);
			bs_tilt_azimuth_angle_LCS[7] = (7 * pi / 16);

			bs_tilt_zenith_angle_LCS[0] = (1 * pi / 8);
			bs_tilt_zenith_angle_LCS[1] = (3 * pi / 8);
			bs_tilt_zenith_angle_LCS[2] = (5 * pi / 8);
			bs_tilt_zenith_angle_LCS[3] = (7 * pi / 8);

			tilt_azimuth_angle_LCS_size = 8;
			tilt_zenith_angle_LCS_size = 4;

			//ue
			ue_tilt_azimuth_angle_LCS[0] = (-3 * pi / 8);
			ue_tilt_azimuth_angle_LCS[1] = (-pi / 8);
			ue_tilt_azimuth_angle_LCS[2] = (pi / 8);
			ue_tilt_azimuth_angle_LCS[3] = (3 * pi / 8);

			ue_tilt_zenith_angle_LCS[0] = (pi / 4);
			ue_tilt_zenith_angle_LCS[1] = (3 * pi / 4);

			ue_tilt_azimuth_angle_LCS_size = 4;
			ue_tilt_zenith_angle_LCS_size = 2;
			
		}
	}
	else if (TYPE == 12)  // Dense Urban
	{
		if (Configuration_Type == 0) // config A
		{
			if (Calibration_mode == 0)
			{
				#if 0
				bs_tilt_azimuth_angle_LCS[0] = (-5 * pi / 16);
				bs_tilt_azimuth_angle_LCS[1] = (-3 * pi / 16);
				bs_tilt_azimuth_angle_LCS[2] = (-pi / 16);
				bs_tilt_azimuth_angle_LCS[3] = (pi / 16);
				bs_tilt_azimuth_angle_LCS[4] = (3 * pi / 16);
				bs_tilt_azimuth_angle_LCS[5] = (5 * pi / 16);

				bs_tilt_zenith_angle_LCS[0] = (5 * pi / 8);
				bs_tilt_zenith_angle_LCS[1] = (7 * pi / 8);

				tilt_azimuth_angle_LCS_size = 6;
				tilt_zenith_angle_LCS_size = 2;
				#endif
				bs_tilt_azimuth_angle_LCS[0] = 0.;
				bs_tilt_zenith_angle_LCS[0] = 102. * (pi / 180.);

				tilt_azimuth_angle_LCS_size = 1;
				tilt_zenith_angle_LCS_size = 1;			
			}
			else
			{
				//bs_tilt_azimuth_angle_LCS[0] = 0.;
				//bs_tilt_zenith_angle_LCS[0] = 105. * (pi / 180.);

				bs_tilt_azimuth_angle_LCS[0] = 0.;
				bs_tilt_zenith_angle_LCS[0] = 102. * (pi / 180.);

				tilt_azimuth_angle_LCS_size = 1;
				tilt_zenith_angle_LCS_size = 1;	

			}

		}
		else
		{
			if (Calibration_mode == 1)
			{
				// Full Calibration: single fixed beam at mechanical bore-sight
				// TR 38.901: "panning (0,0) deg" = no electronic steering,
				// total beam direction determined by mechanical tilt only.
				// Same convention as UMa (Configuration_Type==0) calibration.
				//bs
				bs_tilt_azimuth_angle_LCS[0] = 0.;
				bs_tilt_zenith_angle_LCS[0] = 102. * (pi / 180.);

				tilt_azimuth_angle_LCS_size = 1;
				tilt_zenith_angle_LCS_size = 1;

				// [DISABLED] Previous 32-beam DFT grid (8 azimuth x 4 zenith)
				// — caused ~10 dB CL / SIR mismatch because none of the DFT
				//   beams aligned with the calibration bore-sight (theta_LCS=102°),
				//   and interfering BSs selected random beams from this grid.
				//// bs
				//bs_tilt_azimuth_angle_LCS[0] = (-7 * pi / 16);
				//bs_tilt_azimuth_angle_LCS[1] = (-5 * pi / 16);
				//bs_tilt_azimuth_angle_LCS[2] = (-3 * pi / 16);
				//bs_tilt_azimuth_angle_LCS[3] = (-1 * pi / 16);
				//bs_tilt_azimuth_angle_LCS[4] = (1 * pi / 16);
				//bs_tilt_azimuth_angle_LCS[5] = (3 * pi / 16);
				//bs_tilt_azimuth_angle_LCS[6] = (5 * pi / 16);
				//bs_tilt_azimuth_angle_LCS[7] = (7 * pi / 16);
				//bs_tilt_zenith_angle_LCS[0] = (1 * pi / 8);
				//bs_tilt_zenith_angle_LCS[1] = (3 * pi / 8);
				//bs_tilt_zenith_angle_LCS[2] = (5 * pi / 8);
				//bs_tilt_zenith_angle_LCS[3] = (7 * pi / 8);
				//tilt_azimuth_angle_LCS_size = 8;
				//tilt_zenith_angle_LCS_size = 4;
				//// ue
				//ue_tilt_azimuth_angle_LCS[0] = (-3 * pi / 8);
				//ue_tilt_azimuth_angle_LCS[1] = (-pi / 8);
				//ue_tilt_azimuth_angle_LCS[2] = (pi / 8);
				//ue_tilt_azimuth_angle_LCS[3] = (3 * pi / 8);
				//ue_tilt_zenith_angle_LCS[0] = (pi / 4);
				//ue_tilt_zenith_angle_LCS[1] = (3 * pi / 4);
				//ue_tilt_azimuth_angle_LCS_size = 4;
				//ue_tilt_zenith_angle_LCS_size = 2;
			}
			else
			{
				//bs
				bs_tilt_azimuth_angle_LCS[0] = (-7 * pi / 16);
				bs_tilt_azimuth_angle_LCS[1] = (-5 * pi / 16);
				bs_tilt_azimuth_angle_LCS[2] = (-3 * pi / 16);
				bs_tilt_azimuth_angle_LCS[3] = (-1 * pi / 16);
				bs_tilt_azimuth_angle_LCS[4] = (1 * pi / 16);
				bs_tilt_azimuth_angle_LCS[5] = (3 * pi / 16);
				bs_tilt_azimuth_angle_LCS[6] = (5 * pi / 16);
				bs_tilt_azimuth_angle_LCS[7] = (7 * pi / 16);

				bs_tilt_zenith_angle_LCS[0] = (1 * pi / 8);
				bs_tilt_zenith_angle_LCS[1] = (3 * pi / 8);
				bs_tilt_zenith_angle_LCS[2] = (5 * pi / 8);
				bs_tilt_zenith_angle_LCS[3] = (7 * pi / 8);

				tilt_azimuth_angle_LCS_size = 8;
				tilt_zenith_angle_LCS_size = 4;

				//ue
				ue_tilt_azimuth_angle_LCS[0] = (-3 * pi / 8);
				ue_tilt_azimuth_angle_LCS[1] = (-pi / 8);
				ue_tilt_azimuth_angle_LCS[2] = (pi / 8);
				ue_tilt_azimuth_angle_LCS[3] = (3 * pi / 8);

				ue_tilt_zenith_angle_LCS[0] = (pi / 4);
				ue_tilt_zenith_angle_LCS[1] = (3 * pi / 4);

				ue_tilt_azimuth_angle_LCS_size = 4;
				ue_tilt_zenith_angle_LCS_size = 2;
			}
		}
		
	}
	else if (TYPE == 13) // Rural
	{
		if (Configuration_Type == 2) // config C
		{
			bs_tilt_azimuth_angle_LCS[0] = 0.;
			bs_tilt_zenith_angle_LCS[0] = 96. * (pi / 180.);

			tilt_azimuth_angle_LCS_size = 1;
			tilt_zenith_angle_LCS_size = 1;
		}
		else
		{
			//bs
			bs_tilt_azimuth_angle_LCS[0]   = (-7 * pi / 16);
			bs_tilt_azimuth_angle_LCS[1]   = (-5 * pi / 16);
			bs_tilt_azimuth_angle_LCS[2]   = (-3 * pi / 16);
			bs_tilt_azimuth_angle_LCS[3]   = (-1 * pi / 16);
			bs_tilt_azimuth_angle_LCS[4]   = (1 * pi / 16);
			bs_tilt_azimuth_angle_LCS[5]   = (3 * pi / 16);
			bs_tilt_azimuth_angle_LCS[6]   = (5 * pi / 16);
			bs_tilt_azimuth_angle_LCS[7]   = (7 * pi / 16);

			bs_tilt_zenith_angle_LCS[0]    = (1 * pi / 8);
			bs_tilt_zenith_angle_LCS[1]    = (3 * pi / 8);
			bs_tilt_zenith_angle_LCS[2]    = (5 * pi / 8);
			bs_tilt_zenith_angle_LCS[3]    = (7 * pi / 8);

			tilt_azimuth_angle_LCS_size    = 8;
			tilt_zenith_angle_LCS_size     = 4;

			//ue
			ue_tilt_azimuth_angle_LCS[0]   = (-3 * pi / 8);
			ue_tilt_azimuth_angle_LCS[1]   = (-pi / 8);
			ue_tilt_azimuth_angle_LCS[2]   = (pi / 8);
			ue_tilt_azimuth_angle_LCS[3]   = (3 * pi / 8);

			ue_tilt_zenith_angle_LCS[0]    = (pi / 4);
			ue_tilt_zenith_angle_LCS[1]    = (3 * pi / 4);

			ue_tilt_azimuth_angle_LCS_size = 4;
			ue_tilt_zenith_angle_LCS_size  = 2;
		}
	}

	// Row-beam zenith grid override: replace the scenario's hardcoded zenith beam
	// set with B uniformly spaced beams in [min, max] deg. Azimuth grid unchanged.
	// Only active with row_beam_enable, so the legacy grids above stay bit-exact.
	// (InH is rejected at cfg parse time, so this only reaches the DU/Rural grids.)
	if (row_beam_enable && row_beam_num_zenith > 0)
	{
		int B = row_beam_num_zenith;
		for (int zi = 0; zi < B; zi++)
		{
			Real deg = (B == 1)
			         ? (row_beam_zenith_min_deg + row_beam_zenith_max_deg) / REAL(2.0)
			         : row_beam_zenith_min_deg
			           + (row_beam_zenith_max_deg - row_beam_zenith_min_deg) * zi / REAL(B - 1);
			bs_tilt_zenith_angle_LCS[zi] = deg * (pi / REAL(180.0));
		}
		tilt_zenith_angle_LCS_size = B;
	}

	// Azimuth grid index nearest boresight (0 rad in LCS), used by az_mode==1.
	// Tie (symmetric grids have no exact 0) resolves to the lower index.
	if (row_beam_enable)
	{
		row_beam_boresight_a = 0;
		for (int ai = 1; ai < tilt_azimuth_angle_LCS_size; ai++)
			if (fabs(bs_tilt_azimuth_angle_LCS[ai]) < fabs(bs_tilt_azimuth_angle_LCS[row_beam_boresight_a]))
				row_beam_boresight_a = ai;
	}

	// Handheld mode: single UE beam direction (no UE beam search needed)
	// Handheld per-port pattern is computed in Get_UE_antenna_pattern, not through BF weights
	if (handheld_mode) {
		ue_tilt_azimuth_angle_LCS[0] = 0.;
		ue_tilt_zenith_angle_LCS[0]  = pi / 2.;  // 90° broadside (placeholder, not used by handheld)
		ue_tilt_azimuth_angle_LCS_size = 1;
		ue_tilt_zenith_angle_LCS_size  = 1;
	}

	// ====================================================================
	// Validate beam quantization parameters
	// ====================================================================
	if (tilt_azimuth_angle_LCS_size > MAX_RSRP_SEC_AZIMUTH) {
		cout << "ERROR: tilt_azimuth_angle_LCS_size (" << tilt_azimuth_angle_LCS_size
		     << ") exceeds MAX_RSRP_SEC_AZIMUTH (" << MAX_RSRP_SEC_AZIMUTH << ")" << endl;
		exit(1);
	}
	if (tilt_zenith_angle_LCS_size > MAX_RSRP_SEC_ZENITH) {
		cout << "ERROR: tilt_zenith_angle_LCS_size (" << tilt_zenith_angle_LCS_size
		     << ") exceeds MAX_RSRP_SEC_ZENITH (" << MAX_RSRP_SEC_ZENITH << ")" << endl;
		exit(1);
	}
	if (ue_tilt_azimuth_angle_LCS_size > MAX_RSRP_UE_AZIMUTH) {
		cout << "ERROR: ue_tilt_azimuth_angle_LCS_size (" << ue_tilt_azimuth_angle_LCS_size
		     << ") exceeds MAX_RSRP_UE_AZIMUTH (" << MAX_RSRP_UE_AZIMUTH << ")" << endl;
		exit(1);
	}
	if (ue_tilt_zenith_angle_LCS_size > MAX_RSRP_UE_ZENITH) {
		cout << "ERROR: ue_tilt_zenith_angle_LCS_size (" << ue_tilt_zenith_angle_LCS_size
		     << ") exceeds MAX_RSRP_UE_ZENITH (" << MAX_RSRP_UE_ZENITH << ")" << endl;
		exit(1);
	}

	int K = BS_M / BS_Mp;    /////   4 / 4 = 1
	int L = BS_N / BS_Np;   //////   4/4 = 1

	Real dH;
	Real dV;

	dH = BS_dH;
	dV = BS_dV;


	for (int zenith_idx = 0; zenith_idx < tilt_zenith_angle_LCS_size; zenith_idx++)
	{
		for (int k = 0; k < K; k++) // vertical element
		{
			w[zenith_idx][k] = Generate_vertical_virtualization_weight(bs_tilt_zenith_angle_LCS[zenith_idx], dV, Real(k), Real(K));
		}
	}

	for (int zenith_idx = 0; zenith_idx < tilt_zenith_angle_LCS_size; zenith_idx++)
	{
		for (int azimuth_idx = 0; azimuth_idx < tilt_azimuth_angle_LCS_size; azimuth_idx++)
		{
			for (int l = 0; l < L; l++) // horizontal element
			{
				v[zenith_idx][azimuth_idx][l] = Generate_horizontal_virtualization_weight(bs_tilt_zenith_angle_LCS[zenith_idx], bs_tilt_azimuth_angle_LCS[azimuth_idx], dH, Real(l), Real(L));
			}
		}
	}
	for (int zenith_idx = 0; zenith_idx<tilt_zenith_angle_LCS_size; zenith_idx++)
	{
		for (int azimuth_idx = 0; azimuth_idx<tilt_azimuth_angle_LCS_size; azimuth_idx++)
		{
			for (int k = 0; k < K; k++) // vertical element
			{
				for (int l = 0; l < L; l++) // horizontal element
				{
					virtualization_weight_wv[zenith_idx][azimuth_idx][k][l] = w[zenith_idx][k] * v[zenith_idx][azimuth_idx][l];
				}
			}
		}
	}


	if (ue_antenna_element_gain != 0)  // UE antenna
	{
		int ue_K = MS_M / MS_Mp;    /////   4 / 4 = 1
		int ue_L = MS_N / MS_Np;   //////   4/4 = 1

		Real ue_dH;
		Real ue_dV;

		ue_dH = MS_dH;
		ue_dV = MS_dV;


		for (int zenith_idx = 0; zenith_idx < ue_tilt_zenith_angle_LCS_size; zenith_idx++)
		{
			for (int k = 0; k < ue_K; k++) // vertical element
			{
				ue_w[zenith_idx][k] = Generate_vertical_virtualization_weight(ue_tilt_zenith_angle_LCS[zenith_idx], ue_dV, Real(k), Real(ue_K));
			}
		}

		for (int zenith_idx = 0; zenith_idx < ue_tilt_zenith_angle_LCS_size; zenith_idx++)
		{
			for (int azimuth_idx = 0; azimuth_idx < ue_tilt_azimuth_angle_LCS_size; azimuth_idx++)
			{
				for (int l = 0; l < ue_L; l++) // horizontal element
				{
					ue_v[zenith_idx][azimuth_idx][l] = Generate_horizontal_virtualization_weight(ue_tilt_zenith_angle_LCS[zenith_idx], ue_tilt_azimuth_angle_LCS[azimuth_idx], ue_dH, Real(l), Real(ue_L));
				}
			}
		}
		for (int zenith_idx = 0; zenith_idx<ue_tilt_zenith_angle_LCS_size; zenith_idx++)
		{
			for (int azimuth_idx = 0; azimuth_idx<ue_tilt_azimuth_angle_LCS_size; azimuth_idx++)
			{
				for (int k = 0; k < ue_K; k++) // vertical element
				{
					for (int l = 0; l < ue_L; l++) // horizontal element
					{
						ue_virtualization_weight_wv[zenith_idx][azimuth_idx][k][l] = ue_w[zenith_idx][k] * ue_v[zenith_idx][azimuth_idx][l];
					}
				}
			}
		}
	}

}

ComplexReal  Generate_vertical_virtualization_weight(Real theta, Real dv, Real k, Real K)
{
	ComplexReal  j(0.0, 1.0);
	ComplexReal  w;

	w = REAL(1.0) / sqrt(Real(K))*exp(-j * REAL(2.0) * pi / Wavelength * k * dv*cos(theta));

	return w;
}


ComplexReal  Generate_horizontal_virtualization_weight(Real theta, Real phi, Real dH, Real l, Real L)
{
	ComplexReal  j(0.0, 1.0);
	ComplexReal  v;

	v = REAL(1.0) / sqrt(Real(L))*exp(-j * REAL(2.0) * pi / Wavelength * l * dH * sin(theta) * sin(phi));

	return v;
}

