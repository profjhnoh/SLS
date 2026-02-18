#include <iostream>
#include <fstream>
#include <algorithm>

using namespace std;

#include "const.h"
#include "common.h"
#include "Logging_point.h"

/*===================================================================
FUNCTION: Averaging_Sinr()


===================================================================*/
void Averaging_Sinr(int ue_idx, Real sinr)
{
	if (avr_sinr == NULL)
	{
		avr_sinr = new Real[num_drop*num_MS];
		for (int ue_idx = 0; ue_idx < num_drop*num_MS; ue_idx++)
		{
			avr_sinr[ue_idx] = 0;
		}

		num_of_added = new int[num_drop*num_MS];
		for (int ue_idx = 0; ue_idx < num_drop*num_MS; ue_idx++)
		{
			num_of_added[ue_idx] = 0;
		}
	}

	if (!isinf(sinr)) // sinr not infinite
	{
		avr_sinr[drop_idx*num_MS + ue_idx] += sinr;
		num_of_added[drop_idx*num_MS + ue_idx]++;

	}
}

/*===================================================================
FUNCTION: Coupling_Loss_Get()


===================================================================*/
void Coupling_Loss_Get(Real cp_ue, int ue_idx)
{
	if (drop_idx == 0 && ue_idx == 0)
	{
		coupling_loss = new Real[num_drop*num_MS];
	}

	coupling_loss[ue_idx + drop_idx * num_MS] = cp_ue;
	//cout  << coupling_loss[ue_idx + drop_idx * num_MS] << endl;
	ms[ue_idx].coupling_loss = cp_ue;
}


/*===================================================================
FUNCTION: Coupling_Loss_CDF()


===================================================================*/
void Coupling_Loss_CDF(void)
{
	STASTICS cp_stastics;
	if (scenario == 11)  ///////InH
	{
		if (num_Indoor_TRxP == 1)
		{
			if (Configuration_Type == 0)
			{
				if (Channel_Model_Type == 0)
				{
					cp_stastics.setup(coupling_loss, num_drop*num_MS, (char*)"DL_CP_loss_InH_config_A_Model_A");
				}
				else
				{
					cp_stastics.setup(coupling_loss, num_drop*num_MS, (char*)"DL_CP_loss_InH_config_A_Model_B");
				}
			}
			else if (Configuration_Type == 1)
			{
				if (Channel_Model_Type == 0)
				{
					cp_stastics.setup(coupling_loss, num_drop*num_MS, (char*)"DL_CP_loss_InH_config_B_Model_A");
				}
				else
				{
					cp_stastics.setup(coupling_loss, num_drop*num_MS, (char*)"DL_CP_loss_InH_config_B_Model_B");
				}

			}
			else if (Configuration_Type == 2)
			{
				if (Channel_Model_Type == 0)
				{
					cp_stastics.setup(coupling_loss, num_drop*num_MS, (char*)"DL_CP_loss_InH_config_C_Model_A");
				}
				else
				{
					cp_stastics.setup(coupling_loss, num_drop*num_MS, (char*)"DL_CP_loss_InH_config_C_Model_B");
				}

			}
		}
		else if (num_Indoor_TRxP == 3)
		{
			if (Configuration_Type == 0)
			{
				if (Channel_Model_Type == 0)
				{
					cp_stastics.setup(coupling_loss, num_drop*num_MS, (char*)"DL_CP_loss_InH_36TRxP_config_A_Model_A");
				}
				else if (Channel_Model_Type == 2)
				{
					cp_stastics.setup(coupling_loss, num_drop*num_MS, (char*)"DL_CP_loss_InH_36TRxP_config_A_Model_ETRI");
				}
				else
				{
					cp_stastics.setup(coupling_loss, num_drop*num_MS, (char*)"DL_CP_loss_InH_36TRxP_config_A_Model_B");
				}
			}
			else if (Configuration_Type == 1)
			{
				if (Channel_Model_Type == 0)
				{
					cp_stastics.setup(coupling_loss, num_drop*num_MS, (char*)"DL_CP_loss_InH_36TRxP_config_B_Model_A");
				}
				else if (Channel_Model_Type == 2)
				{
					cp_stastics.setup(coupling_loss, num_drop*num_MS, (char*)"DL_CP_loss_InH_36TRxP_config_B_Model_ETRI");
				}
				else
				{
					cp_stastics.setup(coupling_loss, num_drop*num_MS, (char*)"DL_CP_loss_InH_36TRxP_config_B_Model_B");
				}
			}
			else if (Configuration_Type == 2)
			{
				if (Channel_Model_Type == 0)
				{
					cp_stastics.setup(coupling_loss, num_drop*num_MS, (char*)"DL_CP_loss_InH_36TRxP_config_C_Model_A");
				}
				else if (Channel_Model_Type == 2)
				{
					cp_stastics.setup(coupling_loss, num_drop*num_MS, (char*)"DL_CP_loss_InH_36TRxP_config_C_Model_ETRI");
				}				
				else
				{
					cp_stastics.setup(coupling_loss, num_drop*num_MS, (char*)"DL_CP_loss_InH_36TRxP_config_C_Model_B");
				}

			}
		}
		else
		{
			cout << " = num_Indoor_TRxP = ERROR! --> Coupling_Loss_CDF " << endl;
			getchar();
		}



	}
	else if (scenario == 12)  ///////// DenseUrban
	{
		if (Configuration_Type == 0)
		{
			if (Channel_Model_Type == 0)
			{
				cp_stastics.setup(coupling_loss, num_drop*num_MS, (char*)"DL_CP_loss_DenseUrban_config_A_Model_A");
			}
			else if (Channel_Model_Type == 2)
			{
				cp_stastics.setup(coupling_loss, num_drop*num_MS, (char*)"DL_CP_loss_DenseUrban_config_A_Model_ETRI");
			}			
			else
			{
				cp_stastics.setup(coupling_loss, num_drop*num_MS, (char*)"DL_CP_loss_DenseUrban_config_A_Model_B");
			}
		}
		else if (Configuration_Type == 1)
		{
			if (Channel_Model_Type == 0)
			{
				cp_stastics.setup(coupling_loss, num_drop*num_MS, (char*)"DL_CP_loss_DenseUrban_config_B_Model_A");
			}
			else if (Channel_Model_Type == 2)
			{
				cp_stastics.setup(coupling_loss, num_drop*num_MS, (char*)"DL_CP_loss_DenseUrban_config_B_Model_ETRI");
			}				
			else
			{
				cp_stastics.setup(coupling_loss, num_drop*num_MS, (char*)"DL_CP_loss_DenseUrban_config_B_Model_B");
			}

		}
		else if (Configuration_Type == 2 || Configuration_Type == 3)
		{
			if (Channel_Model_Type == 0)
			{
				cp_stastics.setup(coupling_loss, num_drop*num_MS, (char*)"DL_CP_loss_DenseUrban_config_C_Model_A");
			}
			else if (Channel_Model_Type == 2)
			{
				cp_stastics.setup(coupling_loss, num_drop*num_MS, (char*)"DL_CP_loss_DenseUrban_config_C_Model_ETRI");
			}				
			else
			{
				cp_stastics.setup(coupling_loss, num_drop*num_MS, (char*)"DL_CP_loss_DenseUrban_config_C_Model_B");
			}

		}
	}
	else if (scenario == 13)   ////// Rural
	{
		if (Configuration_Type == 0)
		{
			if (Channel_Model_Type == 0)
			{
				cp_stastics.setup(coupling_loss, num_drop*num_MS, (char*)"DL_CP_loss_Rural_config_A_Model_A");
			}
			else
			{
				cp_stastics.setup(coupling_loss, num_drop*num_MS, (char*)"DL_CP_loss_Rural_config_A_Model_B");
			}
		}
		else if (Configuration_Type == 1)
		{
			if (Channel_Model_Type == 0)
			{
				cp_stastics.setup(coupling_loss, num_drop*num_MS, (char*)"DL_CP_loss_Rural_config_B_Model_A");
			}
			else
			{
				cp_stastics.setup(coupling_loss, num_drop*num_MS, (char*)"DL_CP_loss_Rural_config_B_Model_B");
			}

		}
		else if (Configuration_Type == 2)
		{
			if (Channel_Model_Type == 0)
			{
				cp_stastics.setup(coupling_loss, num_drop*num_MS, (char*)"DL_CP_loss_Rural_config_C_Model_A");
			}
			else
			{
				cp_stastics.setup(coupling_loss, num_drop*num_MS, (char*)"DL_CP_loss_Rural_config_C_Model_B");
			}

		}
	}
	////////////////////////cp_stastics.setup(coupling_loss, num_drop*num_MS, (char*)"Coupling_loss_DL");
	cp_stastics.get_cdf();
	cp_stastics.print_cdf();
	cout << " check Coupling Loss CDF " << endl;
}


/*===================================================================
FUNCTION: Geometry_CDF()



===================================================================*/
void Get_Geometry(Real geometry_ue, int ue_idx)
{
	if (drop_idx == 0 && ue_idx == 0)
	{
		geometry = new Real[num_drop*num_MS];
	}

	geometry[ue_idx + drop_idx * num_MS] = geometry_ue;
	ms[ue_idx].geometry = geometry_ue;
	//cout << "get geometry[" << ue_idx + drop_idx * num_MS <<  "] = " <<  geometry[ue_idx + drop_idx * num_MS] << endl;
}

/*===================================================================
FUNCTION: Get_Wideband_SIR() / Wideband_SIR_CDF()
  3GPP TR 38.901 calibration: Wideband SIR = S/I (without noise)
===================================================================*/
void Get_Wideband_SIR(Real sir_ue, int ue_idx)
{
	if (drop_idx == 0 && ue_idx == 0)
	{
		wideband_sir = new Real[num_drop*num_MS];
	}

	wideband_sir[ue_idx + drop_idx * num_MS] = sir_ue;
}

void Wideband_SIR_CDF(void)
{
	STASTICS sir_stastics;
	sir_stastics.setup(wideband_sir, num_drop*num_MS, (char*)"DL_Wideband_SIR");
	sir_stastics.get_cdf();
	sir_stastics.print_cdf();
	cout << " check Wideband SIR CDF " << endl;
}

/*===================================================================
FUNCTION: Geometry_CDF()



===================================================================*/
void Geometry_CDF(void)
{
	STASTICS geo_stastics;
	
	//geo_stastics.setup(geometry, num_drop*num_MS, (char*)"Geometry_DL");

	if (TYPE == 11)  ///////InH
	{

		if (num_Indoor_TRxP == 1)
		{
			if (Configuration_Type == 0)
			{
				if (Channel_Model_Type == 0)
				{
					geo_stastics.setup(geometry, num_drop*num_MS, (char*)"DL_Geometry_InH_config_A_Model_A");
				}
				else
				{
					geo_stastics.setup(geometry, num_drop*num_MS, (char*)"DL_Geometry_InH_config_A_Model_B");
				}
			}
			else if (Configuration_Type == 1)
			{
				if (Channel_Model_Type == 0)
				{
					geo_stastics.setup(geometry, num_drop*num_MS, (char*)"DL_Geometry_InH_config_B_Model_A");
				}
				else
				{
					geo_stastics.setup(geometry, num_drop*num_MS, (char*)"DL_Geometry_InH_config_B_Model_B");
				}

			}
			else if (Configuration_Type == 2)
			{
				if (Channel_Model_Type == 0)
				{
					geo_stastics.setup(geometry, num_drop*num_MS, (char*)"DL_Geometry_InH_config_C_Model_A");
				}
				else
				{
					geo_stastics.setup(geometry, num_drop*num_MS, (char*)"DL_Geometry_InH_config_C_Model_B");
				}

			}
		}
		else if (num_Indoor_TRxP == 3)
		{
			if (Configuration_Type == 0)
			{
				if (Channel_Model_Type == 0)
				{
					geo_stastics.setup(geometry, num_drop*num_MS, (char*)"DL_Geometry_InH_36TRxP_config_A_Model_A");
				}
				else
				{
					geo_stastics.setup(geometry, num_drop*num_MS, (char*)"DL_Geometry_InH_36TRxP_config_A_Model_B");
				}
			}
			else if (Configuration_Type == 1)
			{
				if (Channel_Model_Type == 0)
				{
					geo_stastics.setup(geometry, num_drop*num_MS, (char*)"DL_Geometry_InH_36TRxP_config_B_Model_A");
				}
				else
				{
					geo_stastics.setup(geometry, num_drop*num_MS, (char*)"DL_Geometry_InH_36TRxP_config_B_Model_B");
				}

			}
			else if (Configuration_Type == 2)
			{
				if (Channel_Model_Type == 0)
				{
					geo_stastics.setup(geometry, num_drop*num_MS, (char*)"DL_Geometry_InH_36TRxP_config_C_Model_A");
				}
				else
				{
					geo_stastics.setup(geometry, num_drop*num_MS, (char*)"DL_Geometry_InH_36TRxP_config_C_Model_B");
				}

			}
		}
		else
		{
			cout << " = num_Indoor_TRxP = ERROR! --> Coupling_Loss_CDF " << endl;
			getchar();
		}


		
	}
	else if (TYPE == 12)  ///////// DenseUrban
	{
		if (Configuration_Type == 0)
		{
			if (Channel_Model_Type == 0)
			{
				geo_stastics.setup(geometry, num_drop*num_MS, (char*)"DL_Geometry_DenseUrban_config_A_Model_A");
			}
			else
			{
				geo_stastics.setup(geometry, num_drop*num_MS, (char*)"DL_Geometry_DenseUrban_config_A_Model_B");
			}
		}
		else if (Configuration_Type == 1)
		{
			if (Channel_Model_Type == 0)
			{
				geo_stastics.setup(geometry, num_drop*num_MS, (char*)"DL_Geometry_DenseUrban_config_B_Model_A");
			}
			else
			{
				geo_stastics.setup(geometry, num_drop*num_MS, (char*)"DL_Geometry_DenseUrban_config_B_Model_B");
			}

		}
		else if (Configuration_Type == 2)
		{
			if (Channel_Model_Type == 0)
			{
				geo_stastics.setup(geometry, num_drop*num_MS, (char*)"DL_Geometry_DenseUrban_config_C_Model_A");
			}
			else
			{
				geo_stastics.setup(geometry, num_drop*num_MS, (char*)"DL_Geometry_DenseUrban_config_C_Model_B");
			}
		}
		else if (Configuration_Type == 3)
		{
			if (Channel_Model_Type == 0)
			{
				geo_stastics.setup(geometry, num_drop*num_MS, (char*)"DL_Geometry_DenseUrban_Micro_config_C_Model_A");
			}
			else
			{
				geo_stastics.setup(geometry, num_drop*num_MS, (char*)"DL_Geometry_DenseUrban_Micro_config_C_Model_B");
			}

		}		
	}
	else if (TYPE == 13)   ////// Rural
	{
		if (Configuration_Type == 0)
		{
			if (Channel_Model_Type == 0)
			{
				geo_stastics.setup(geometry, num_drop*num_MS, (char*)"DL_Geometry_Rural_config_A_Model_A");
			}
			else
			{
				geo_stastics.setup(geometry, num_drop*num_MS, (char*)"DL_Geometry_Rural_config_A_Model_B");
			}
		}
		else if (Configuration_Type == 1)
		{
			if (Channel_Model_Type == 0)
			{
				geo_stastics.setup(geometry, num_drop*num_MS, (char*)"DL_Geometry_Rural_config_B_Model_A");
			}
			else
			{
				geo_stastics.setup(geometry, num_drop*num_MS, (char*)"DL_Geometry_Rural_config_B_Model_B");
			}

		}
		else if (Configuration_Type == 2)
		{
			if (Channel_Model_Type == 0)
			{
				geo_stastics.setup(geometry, num_drop*num_MS, (char*)"DL_Geometry_Rural_config_C_Model_A");
			}
			else
			{
				geo_stastics.setup(geometry, num_drop*num_MS, (char*)"DL_Geometry_Rural_config_C_Model_B");
			}

		}
	}

	geo_stastics.get_cdf();
	geo_stastics.print_cdf();
	
	cout << " check Geometry CDF " << endl;
}


/*===================================================================
FUNCTION: Geometry_CDF()



===================================================================*/
void Get_RMS_delay_spread(Real RMS_delay_spread_ue, int ue_idx)
{
	if (drop_idx == 0 && ue_idx == 0)
	{
		RMS_delay_spread = new Real[num_drop*num_MS];
	}

	RMS_delay_spread[ue_idx + drop_idx * num_MS] = RMS_delay_spread_ue;
	//cout << "get geometry[" << ue_idx + drop_idx * num_MS <<  "] = " <<  geometry[ue_idx + drop_idx * num_MS] << endl;
}

/*===================================================================
FUNCTION: Geometry_CDF()



===================================================================*/
void RMS_delay_spread_CDF(void)
{
	STASTICS RMS_stastics;

	RMS_stastics.setup(RMS_delay_spread, num_drop*num_MS, (char*)"RMS_delay_spread_DL");
	RMS_stastics.get_cdf();
	RMS_stastics.print_cdf();

}



/*===================================================================
FUNCTION: Geometry_CDF()


===================================================================*/
void Get_circular_angle_spread_AOA(Real angle_spread_AOA_ue, int ue_idx)
{
	if (drop_idx == 0 && ue_idx == 0)
	{
		circular_angle_spread_AOA = new Real[num_drop*num_MS];
	}

	circular_angle_spread_AOA[ue_idx + drop_idx * num_MS] = angle_spread_AOA_ue;
	//cout << "get geometry[" << ue_idx + drop_idx * num_MS <<  "] = " <<  geometry[ue_idx + drop_idx * num_MS] << endl;
}


void circular_angle_spread_AOA_CDF(void)
{
	STASTICS AOA_stastics;

	AOA_stastics.setup(circular_angle_spread_AOA, num_drop*num_MS, (char*)"angle_spread_AOA_DL");
	AOA_stastics.get_cdf();
	AOA_stastics.print_cdf();

}

void Get_circular_angle_spread_AOD(Real angle_spread_AOD_ue, int ue_idx)
{
	if (drop_idx == 0 && ue_idx == 0)
	{
		circular_angle_spread_AOD = new Real[num_drop*num_MS];
	}

	circular_angle_spread_AOD[ue_idx + drop_idx * num_MS] = angle_spread_AOD_ue;
	//cout << "get geometry[" << ue_idx + drop_idx * num_MS <<  "] = " <<  geometry[ue_idx + drop_idx * num_MS] << endl;
}

void circular_angle_spread_AOD_CDF(void)
{
	STASTICS AOD_stastics;

	AOD_stastics.setup(circular_angle_spread_AOD, num_drop*num_MS, (char*)"angle_spread_AOD_DL");
	AOD_stastics.get_cdf();
	AOD_stastics.print_cdf();

}


/*===================================================================
FUNCTION: Collect_LSP_from_ServingCell()
Description: Collect LSP values (DS, ASD, ASA, ZSD, ZSA) from each
             UE's serving cell based on propagation condition.
===================================================================*/
void Collect_LSP_from_ServingCell()
{
	for (int ue_idx = 0; ue_idx < num_MS; ue_idx++)
	{
		int serv_bs = links[ue_idx].self_bs_idx;
		int prop = channel[serv_bs][ue_idx].Propagation;

		// DS: raw LSP (seconds -> nanoseconds)
		Real ds;
		if (prop == LOS_propagation)       ds = ms[ue_idx].LSPs[serv_bs](LOS_DS);
		else if (prop == NLOS_propagation)  ds = ms[ue_idx].LSPs[serv_bs](NLOS_DS);
		else                                ds = ms[ue_idx].LSPs[serv_bs](OUT2IN_DS);
		Get_LSP_DS(ds * 1e9, ue_idx);

		// ASD/ASA/ZSD/ZSA: circular angular spread from SSP clusters (degrees)
		Get_LSP_ASD(channel[serv_bs][ue_idx].circular_angle_spread_AOD, ue_idx);
		Get_LSP_ASA(channel[serv_bs][ue_idx].circular_angle_spread_AOA, ue_idx);
		Get_LSP_ZSD(channel[serv_bs][ue_idx].circular_angle_spread_ZOD, ue_idx);
		Get_LSP_ZSA(channel[serv_bs][ue_idx].circular_angle_spread_ZOA, ue_idx);
	}
}

/*===================================================================
FUNCTION: Get_LSP_DS / LSP_DS_CDF
===================================================================*/
void Get_LSP_DS(Real val, int ue_idx)
{
	if (drop_idx == 0 && ue_idx == 0)
	{
		LSP_DS = new Real[num_drop*num_MS];
	}
	LSP_DS[ue_idx + drop_idx * num_MS] = val;
}

void LSP_DS_CDF(void)
{
	STASTICS stat;
	stat.setup(LSP_DS, num_drop*num_MS, (char*)"LSP_DS_ns_DL");
	stat.get_cdf();
	stat.print_cdf();
	cout << " check LSP DS CDF " << endl;
}

/*===================================================================
FUNCTION: Get_LSP_ASD / LSP_ASD_CDF
===================================================================*/
void Get_LSP_ASD(Real val, int ue_idx)
{
	if (drop_idx == 0 && ue_idx == 0)
	{
		LSP_ASD = new Real[num_drop*num_MS];
	}
	LSP_ASD[ue_idx + drop_idx * num_MS] = val;
}

void LSP_ASD_CDF(void)
{
	STASTICS stat;
	stat.setup(LSP_ASD, num_drop*num_MS, (char*)"LSP_ASD_DL");
	stat.get_cdf();
	stat.print_cdf();
	cout << " check LSP ASD CDF " << endl;
}

/*===================================================================
FUNCTION: Get_LSP_ASA / LSP_ASA_CDF
===================================================================*/
void Get_LSP_ASA(Real val, int ue_idx)
{
	if (drop_idx == 0 && ue_idx == 0)
	{
		LSP_ASA = new Real[num_drop*num_MS];
	}
	LSP_ASA[ue_idx + drop_idx * num_MS] = val;
}

void LSP_ASA_CDF(void)
{
	STASTICS stat;
	stat.setup(LSP_ASA, num_drop*num_MS, (char*)"LSP_ASA_DL");
	stat.get_cdf();
	stat.print_cdf();
	cout << " check LSP ASA CDF " << endl;
}

/*===================================================================
FUNCTION: Get_LSP_ZSD / LSP_ZSD_CDF
===================================================================*/
void Get_LSP_ZSD(Real val, int ue_idx)
{
	if (drop_idx == 0 && ue_idx == 0)
	{
		LSP_ZSD = new Real[num_drop*num_MS];
	}
	LSP_ZSD[ue_idx + drop_idx * num_MS] = val;
}

void LSP_ZSD_CDF(void)
{
	STASTICS stat;
	stat.setup(LSP_ZSD, num_drop*num_MS, (char*)"LSP_ZSD_DL");
	stat.get_cdf();
	stat.print_cdf();
	cout << " check LSP ZSD CDF " << endl;
}

/*===================================================================
FUNCTION: Get_LSP_ZSA / LSP_ZSA_CDF
===================================================================*/
void Get_LSP_ZSA(Real val, int ue_idx)
{
	if (drop_idx == 0 && ue_idx == 0)
	{
		LSP_ZSA = new Real[num_drop*num_MS];
	}
	LSP_ZSA[ue_idx + drop_idx * num_MS] = val;
}

void LSP_ZSA_CDF(void)
{
	STASTICS stat;
	stat.setup(LSP_ZSA, num_drop*num_MS, (char*)"LSP_ZSA_DL");
	stat.get_cdf();
	stat.print_cdf();
	cout << " check LSP ZSA CDF " << endl;
}


/*===================================================================
FUNCTION: Plos()



===================================================================*/
void Plos(void)
{
	int los_count = 0;
	for (int ue_idx = 0; ue_idx < num_MS; ue_idx++)
	{
		/*
		if (los[ue_idx][ms[ue_idx]._sector_in_control] == 1)
		{
			los_count++;
		}
		*/
	}
}

/*===================================================================
FUNCTION: SINR_CDF()


===================================================================*/
void SINR_CDF(void)
{



	for (int ue_idx = 0; ue_idx < num_drop*num_MS; ue_idx++)
	{

		avr_sinr[ue_idx] = linear2dB(avr_sinr[ue_idx] / num_of_added[ue_idx]);
	}

	STASTICS sinr_stastics;
	sinr_stastics.setup(avr_sinr, num_drop*num_MS, (char*)"SINR_DL");
	sinr_stastics.get_cdf();
	sinr_stastics.print_cdf();

}

#ifdef ENABLE_PROGRESSBAR
/*===================================================================
FUNCTION: Progressbar()
Author : jhnoh
Create Date : 22-03-29
===================================================================*/
void ResetProgressbar( int drop_idx )
{
	//bar.set_progress(0);
	char msg[50];
	sprintf(msg,"  Drop %d/%d ", drop_idx+1, num_drop);
	//bar.set_option(option::PrefixText{msg});

	int total = N_pf + run_times + SCHEDULE_DELAY;
    bar = new ProgressBar{
	option::BarWidth{50},
	option::Start{"["},
	option::Fill{"■"},
	option::Lead{"■"},
	option::Remainder{"-"},
	option::End{" ]"},
	option::PrefixText{msg},
	option::ShowElapsedTime{true},
	option::ShowRemainingTime{true},
	option::MaxProgress{total},
	option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}
	};
	old_time_percent = -1;


}

/*===================================================================
FUNCTION: Progressbar()
Author : jhnoh
Create Date : 22-03-29
Modified: Display t value and percentage with 2 decimal places every tick
===================================================================*/
void Progressbar( int t )
{
	int total = N_pf + run_times + SCHEDULE_DELAY;
	Real percent = ((Real)(t + 1) / (Real)total) * 100.0;

	// Update postfix text to show t and percentage with 2 decimal places
	char postfix_msg[100];
	sprintf(postfix_msg, " t=%d (%.2f%%)", t, percent);
	bar->set_option(option::PostfixText{postfix_msg});

	// Set progress to current position (t+1 because we're setting progress after completing t)
	bar->set_progress(t + 1);
}
#endif

/*===================================================================
FUNCTION: Singular Value CDF Collection Functions
Author : Claude
Description: 3GPP compliant PRB Singular Value Ratio CDF collection

             3GPP Definition:
             - PRB singular values = eigenvalues of mean covariance matrix in PRB
             - Output: ratio of (2nd, 3rd, ..., xth) to 1st (largest) SV in dB
             - x = min(BS antenna ports, UE antenna ports)
             - Collected at t=0 (initial time)

             Formula: SV_Ratio_n = 10*log10(SV_n / SV_1) [dB]
===================================================================*/
vector<vector<Real>> singular_value_samples;      // Absolute SV values (for backward compatibility)
vector<vector<Real>> singular_value_ratio_samples; // SV ratios in dB (3GPP definition)
int num_singular_values = 0;

void Init_Singular_Value_Collection()
{
	if (!g_collect_singular_values) return;

	num_singular_values = min(NUM_RX_Port, NUM_TX_Port);
	singular_value_samples.clear();
	singular_value_samples.resize(num_singular_values);
	singular_value_ratio_samples.clear();
	singular_value_ratio_samples.resize(num_singular_values - 1);  // Ratios: SV2/SV1, SV3/SV1, ...

	// Reserve memory to avoid frequent reallocations
	// Expected samples: num_drops * num_MS * num_rb (only at t=0)
	int expected_samples = num_drop * num_MS * num_rb;
	for (int i = 0; i < num_singular_values; i++) {
		singular_value_samples[i].reserve(expected_samples);
	}
	for (int i = 0; i < num_singular_values - 1; i++) {
		singular_value_ratio_samples[i].reserve(expected_samples);
	}

	cout << "Singular Value Collection: ENABLED (3GPP compliant, " << num_singular_values << " SVs, "
	     << (num_singular_values - 1) << " ratio CDFs)" << endl;
}

void Collect_Singular_Values_All()
{
	if (!g_collect_singular_values) return;

	// 3GPP: Collect only at t=0 (initial time)
	if (t != 0) return;

	// DEBUG: Print first UE's H_m matrix and singular values for diagnosis
	static bool debug_printed = false;
	if (!debug_printed && num_MS > 0 && ms[30].H_m != NULL) {
		debug_printed = true;
		cout << "\n=== DEBUG: Singular Value Diagnosis (UE 0, RB 0) ===" << endl;
		MatrixXcReal H = ms[30].H_m[0][0];
		cout << "H_m matrix size: " << H.rows() << " x " << H.cols() << endl;
		cout << "H_m Frobenius norm: " << H.norm() << endl;
		cout << "H_m matrix (first 4x4 block):" << endl;
		int max_r = min((int)H.rows(), 4);
		int max_c = min((int)H.cols(), 4);
		for (int r = 0; r < max_r; r++) {
			cout << "  ";
			for (int c = 0; c < max_c; c++) {
				cout << "(" << H(r,c).real() << "," << H(r,c).imag() << ") ";
			}
			cout << endl;
		}

		// Compute and print singular values
		Eigen::JacobiSVD<MatrixXcReal> svd(H, Eigen::ComputeThinU | Eigen::ComputeThinV);
		auto sv = svd.singularValues();
		cout << "Singular values: ";
		for (int i = 0; i < sv.size(); i++) {
			cout << sv(i) << " ";
		}
		cout << endl;
		cout << "Condition number (SV1/SVn): " << sv(0) / sv(sv.size()-1) << endl;

		// Check multiple RBs
		cout << "\nSV ratios across first 5 RBs:" << endl;
		for (int rb = 0; rb < min(num_rb, num_rb); rb++) {
			MatrixXcReal H_rb = ms[30].H_m[0][rb];
			Eigen::JacobiSVD<MatrixXcReal> svd_rb(H_rb, Eigen::ComputeThinU | Eigen::ComputeThinV);
			auto sv_rb = svd_rb.singularValues();
			cout << "  RB " << rb << ": SV=[";
			for (int i = 0; i < sv_rb.size(); i++) {
				cout << sv_rb(i);
				if (i < sv_rb.size()-1) cout << ", ";
			}
			cout << "], ratio_dB=[";
			for (int i = 1; i < sv_rb.size(); i++) {
				cout << 10.0*log10(sv_rb(i)/sv_rb(0));
				if (i < sv_rb.size()-1) cout << ", ";
			}
			cout << "]" << endl;
		}
		cout << "================================================\n" << endl;
	}

	// Thread-local storage for each thread's collected singular values
	// Each thread collects into its own vector to avoid contention
	int num_threads = num_of_threads;

	// Pre-allocate thread-local storage for absolute SVs and ratios
	vector<vector<vector<Real>>> thread_local_samples(num_threads);
	vector<vector<vector<Real>>> thread_local_ratios(num_threads);
	for (int tid = 0; tid < num_threads; tid++) {
		thread_local_samples[tid].resize(num_singular_values);
		thread_local_ratios[tid].resize(num_singular_values - 1);
		// Estimate samples per thread
		int samples_per_thread = (num_MS * num_rb) / num_threads + 1;
		for (int sv = 0; sv < num_singular_values; sv++) {
			thread_local_samples[tid][sv].reserve(samples_per_thread);
		}
		for (int r = 0; r < num_singular_values - 1; r++) {
			thread_local_ratios[tid][r].reserve(samples_per_thread);
		}
	}

	// Parallel collection across all UEs
	#if ENABLE_MULTITHREADING
	#pragma omp parallel num_threads(num_threads)
	{
		int tid = omp_get_thread_num();
		#pragma omp for schedule(dynamic)
	#else
		int tid = 0;
	#endif
		for (int ue_idx = 0; ue_idx < num_MS; ue_idx++) {
			// Check if channel matrix is available
			if (ms[ue_idx].H_m == NULL) continue;

			// Collect singular values from each RB
			for (int rb_idx = 0; rb_idx < num_rb; rb_idx++) {
				MatrixXcReal H = ms[ue_idx].H_m[0][rb_idx];

				// 3GPP: PRB singular values = eigenvalues of mean covariance matrix
				// For single PRB: R = H^H * H (channel correlation matrix at RX side)
				// Singular values of H = sqrt(eigenvalues of H^H * H)
				// Using SVD directly: singular values are already sqrt(eigenvalues)
				Eigen::JacobiSVD<MatrixXcReal> svd(H, Eigen::ComputeThinU | Eigen::ComputeThinV);
				auto sv = svd.singularValues();

				if(std::isnan(sv(0)) || sv(0) > 1000 || sv(0) < -1000 )
				{
					for (int r = 0; r < 2; r++) {
						cout << "  ";
						for (int c = 0; c < 2; c++) {
							cout << "(" << H(r,c).real() << "," << H(r,c).imag() << ") ";
						}
						cout << endl;
					}
					cout << "---------------------------" << endl;
				}

				// Store each singular value in thread-local storage
				int num_sv = min((int)sv.size(), num_singular_values);
				for (int i = 0; i < num_sv; i++) {
					thread_local_samples[tid][i].push_back(sv(i));
				}

				// 3GPP: Compute SV ratios (SV_n / SV_1) in dB scale
				// SV_1 is the largest (sv(0)), SV_2 is second largest (sv(1)), etc.
				if (num_sv > 1 && sv(0) > 1e-10) {  // Avoid division by zero
					Real sv1 = sv(0);  // Largest singular value
					for (int r = 0; r < num_sv - 1; r++) {
						// Ratio: SV_{r+2} / SV_1 in dB
						Real ratio_dB = 20.0 * log10(sv1/sv(r + 1));
						thread_local_ratios[tid][r].push_back(ratio_dB);
					}
				}
			}
		}
	#if ENABLE_MULTITHREADING
	}
	#endif

	// Merge thread-local results into global storage
	for (int tid = 0; tid < num_threads; tid++) {
		// Absolute singular values
		for (int sv = 0; sv < num_singular_values; sv++) {
			singular_value_samples[sv].insert(
				singular_value_samples[sv].end(),
				thread_local_samples[tid][sv].begin(),
				thread_local_samples[tid][sv].end()
			);
		}
		// SV ratios in dB
		for (int r = 0; r < num_singular_values - 1; r++) {
			singular_value_ratio_samples[r].insert(
				singular_value_ratio_samples[r].end(),
				thread_local_ratios[tid][r].begin(),
				thread_local_ratios[tid][r].end()
			);
		}
	}
}

void Singular_Value_CDF()
{
	if (!g_collect_singular_values) return;

	cout << "Computing Singular Value CDFs (3GPP compliant)..." << endl;

	// Output absolute singular values (linear)
	for (int sv_idx = 0; sv_idx < num_singular_values; sv_idx++) {
		if (singular_value_samples[sv_idx].empty()) {
			cout << "  SV_" << (sv_idx + 1) << ": No samples collected" << endl;
			continue;
		}

		char name[32];
		sprintf(name, "SV_%d", sv_idx + 1);

		STASTICS sv_stat;
		sv_stat.setup(singular_value_samples[sv_idx].data(),
		              singular_value_samples[sv_idx].size(), name);
		sv_stat.get_cdf();
		sv_stat.print_cdf();

		cout << "  " << name << "_cdf.dat: " << singular_value_samples[sv_idx].size()
		     << " samples" << endl;
	}

	// Output 3GPP compliant SV ratio CDFs: SV_n / SV_1 in dB
	cout << "Computing SV Ratio CDFs (3GPP: 10*log10(SV_n/SV_1))..." << endl;

	for (int r = 0; r < num_singular_values - 1; r++) {
		if (singular_value_ratio_samples[r].empty()) {
			cout << "  SV_Ratio_" << (r + 2) << "_1: No samples collected" << endl;
			continue;
		}

		char name[64];
		sprintf(name, "SV_Ratio_%d_1_dB", r + 2);  // e.g., SV_Ratio_2_1_dB = SV2/SV1

		STASTICS ratio_stat;
		ratio_stat.setup(singular_value_ratio_samples[r].data(),
		                 singular_value_ratio_samples[r].size(), name);
		ratio_stat.get_cdf();
		ratio_stat.print_cdf();

		cout << "  " << name << "_cdf.dat: " << singular_value_ratio_samples[r].size()
		     << " samples (ratio SV_" << (r + 2) << "/SV_1 in dB)" << endl;
	}
}

/*===================================================================
FUNCTION: Precoding-based Metrics CDF Collection Functions
Author : Claude
Description: Collect wideband SINR and coupling loss with TX precoding (W)
             and RX MMSE receiver (u) applied (full digital beamforming)

             For each UE, accumulate metrics over all time slots, then
             compute time-averaged values for CDF output.

             Called from Receive_DL_mTRP after Compute_RBs_SINR() to reuse
             already-computed SINR values and avoid redundant computation.
===================================================================*/
vector<Real> precoding_wideband_sinr_samples;    // Final CDF samples (time-averaged per UE)
vector<Real> precoding_coupling_loss_samples;    // Final CDF samples (time-averaged per UE)

// Per-UE accumulator arrays for time averaging
// Index: [drop_idx * num_MS + ue_idx]
Real* precoding_sinr_sum = NULL;      // Sum of wideband SINR (linear) over time
Real* precoding_cl_sum = NULL;        // Sum of coupling loss (dB) over time
int*  precoding_sample_count = NULL;  // Number of samples per UE

void Init_Precoding_Metrics()
{
	precoding_wideband_sinr_samples.clear();
	precoding_coupling_loss_samples.clear();

	// Reserve memory for final CDF: num_drops * num_MS
	int expected_samples = num_drop * num_MS;
	precoding_wideband_sinr_samples.reserve(expected_samples);
	precoding_coupling_loss_samples.reserve(expected_samples);

	// Allocate per-UE accumulator arrays
	int total_ues = num_drop * num_MS;
	precoding_sinr_sum = new Real[total_ues]();      // () for zero-initialization
	precoding_cl_sum = new Real[total_ues]();
	precoding_sample_count = new int[total_ues]();

	cout << "Precoding Metrics Collection: ENABLED (time-averaged per UE)" << endl;
}

// Called from Receive_DL_mTRP after SINR computation
// wideband_sinr_linear: average SINR across RBs (linear scale)
// coupling_loss_dB: precoding coupling loss (dB)
void Accumulate_Precoding_Metrics(int ue_idx, Real wideband_sinr_linear, Real coupling_loss_dB)
{
	if (precoding_sinr_sum == NULL) return;  // Not initialized

	int idx = drop_idx * num_MS + ue_idx;

	// Thread-safe accumulation (called from parallel Measure loop)
	#if ENABLE_MULTITHREADING
	#pragma omp atomic
	#endif
	precoding_sinr_sum[idx] += wideband_sinr_linear;

	#if ENABLE_MULTITHREADING
	#pragma omp atomic
	#endif
	precoding_cl_sum[idx] += coupling_loss_dB;

	#if ENABLE_MULTITHREADING
	#pragma omp atomic
	#endif
	precoding_sample_count[idx]++;
}

// Call this at the end of simulation to compute time-averaged CDF
void Finalize_Precoding_Metrics()
{
	int total_ues = num_drop * num_MS;

	for (int idx = 0; idx < total_ues; idx++) {
		if (precoding_sample_count[idx] > 0) {
			// Time-averaged SINR (average in linear domain, then convert to dB)
			Real avg_sinr_linear = precoding_sinr_sum[idx] / precoding_sample_count[idx];
			Real avg_sinr_dB = 10.0 * log10(avg_sinr_linear);
			precoding_wideband_sinr_samples.push_back(avg_sinr_dB);

			// Time-averaged coupling loss (average in dB domain)
			Real avg_cl_dB = precoding_cl_sum[idx] / precoding_sample_count[idx];
			precoding_coupling_loss_samples.push_back(avg_cl_dB);
		}
	}

	// Clean up accumulator arrays
	delete[] precoding_sinr_sum;
	delete[] precoding_cl_sum;
	delete[] precoding_sample_count;
	precoding_sinr_sum = NULL;
	precoding_cl_sum = NULL;
	precoding_sample_count = NULL;
}

void Precoding_SINR_CDF()
{
	if (precoding_wideband_sinr_samples.empty()) {
		cout << "Precoding SINR CDF: No samples collected" << endl;
		return;
	}

	cout << "Computing Precoding Wideband SINR CDF..." << endl;

	STASTICS sinr_stat;
	sinr_stat.setup(precoding_wideband_sinr_samples.data(),
	                precoding_wideband_sinr_samples.size(),
	                (char*)"Precoding_Wideband_SINR");
	sinr_stat.get_cdf();
	sinr_stat.print_cdf();

	cout << "  Precoding_Wideband_SINR_cdf.dat: "
	     << precoding_wideband_sinr_samples.size() << " samples" << endl;
}

void Precoding_Coupling_Loss_CDF()
{
	if (precoding_coupling_loss_samples.empty()) {
		cout << "Precoding Coupling Loss CDF: No samples collected" << endl;
		return;
	}

	cout << "Computing Precoding Coupling Loss CDF..." << endl;

	STASTICS cl_stat;
	cl_stat.setup(precoding_coupling_loss_samples.data(),
	              precoding_coupling_loss_samples.size(),
	              (char*)"Precoding_Coupling_Loss");
	cl_stat.get_cdf();
	cl_stat.print_cdf();

	cout << "  Precoding_Coupling_Loss_cdf.dat: "
	     << precoding_coupling_loss_samples.size() << " samples" << endl;
}


// ====================================================================
// Element-Level Singular Value Collection (ns-3 style)
// ====================================================================
// Collects SVD of H_m_elem[rb](totalRx, totalTx) — element-level
// channel matrices without beamforming, for comparison with 3GPP
// calibration results.
// ====================================================================

static vector<vector<Real>> elem_sv_samples;       // [sv_idx] -> samples
static vector<vector<Real>> elem_sv_ratio_samples;  // [ratio_idx] -> samples (dB)
static int elem_num_sv = 0;
static bool elem_sv_initialized = false;

void Collect_Singular_Values_ElementLevel()
{
	if (!g_collect_singular_values) return;

	int totalTx = BS_M * BS_N * BS_P;
	int totalRx = MS_M * MS_N * MS_P;
	int min_dim = min(totalRx, totalTx);

	// Initialize on first call
	if (!elem_sv_initialized) {
		elem_num_sv = min_dim;
		elem_sv_samples.clear();
		elem_sv_samples.resize(elem_num_sv);
		elem_sv_ratio_samples.clear();
		if (elem_num_sv > 1)
			elem_sv_ratio_samples.resize(elem_num_sv - 1);

		int expected = num_drop * num_MS * num_rb;
		for (int i = 0; i < elem_num_sv; i++)
			elem_sv_samples[i].reserve(expected);
		for (int i = 0; i < elem_num_sv - 1; i++)
			elem_sv_ratio_samples[i].reserve(expected);

		elem_sv_initialized = true;

		cout << "Element-Level SV Collection: ENABLED (" << elem_num_sv << " SVs, "
		     << "totalRx=" << totalRx << ", totalTx=" << totalTx << ")" << endl;
	}

	// Collect from all UEs
	for (int ue_idx = 0; ue_idx < num_MS; ue_idx++) {
		if (ms[ue_idx].H_m_elem == NULL) continue;

		for (int rb_idx = 0; rb_idx < num_rb; rb_idx++) {
			MatrixXcReal& H = ms[ue_idx].H_m_elem[rb_idx];

			if (H.rows() == 0 || H.cols() == 0) continue;

			Eigen::JacobiSVD<MatrixXcReal> svd(H, Eigen::ComputeThinU | Eigen::ComputeThinV);
			auto sv = svd.singularValues();

			int num_sv_avail = min((int)sv.size(), elem_num_sv);
			for (int i = 0; i < num_sv_avail; i++) {
				// Standard: PRB singular values = eigenvalues of H^H·H = sigma^2
				elem_sv_samples[i].push_back(sv(i) * sv(i));
			}

			// SV ratios in dB: 10*log10(lambda_1 / lambda_n)  (power domain)
			if (num_sv_avail > 1 && sv(0) > 1e-10) {
				Real lambda1 = sv(0) * sv(0);
				for (int r = 0; r < num_sv_avail - 1; r++) {
					Real lambda_n = sv(r + 1) * sv(r + 1);
					Real ratio_dB = 10.0 * log10(lambda1 / lambda_n);
					elem_sv_ratio_samples[r].push_back(ratio_dB);
				}
			}
		}
	}

	// Debug: Print first UE info
	if (ms[0].H_m_elem != NULL) {
		MatrixXcReal& H0 = ms[0].H_m_elem[0];
		cout << "  Element-level H_m size: " << H0.rows() << " x " << H0.cols()
		     << ", Frobenius norm: " << H0.norm() << endl;

		Eigen::JacobiSVD<MatrixXcReal> svd0(H0, Eigen::ComputeThinU | Eigen::ComputeThinV);
		auto sv0 = svd0.singularValues();
		cout << "  Element-level eigenvalues of H^H*H (UE0, RB0): ";
		for (int i = 0; i < sv0.size(); i++) cout << sv0(i) * sv0(i) << " ";
		cout << endl;
	}
}

void Singular_Value_ElementLevel_CDF()
{
	if (!g_collect_singular_values || !elem_sv_initialized) return;

	cout << "Computing Element-Level Singular Value CDFs..." << endl;

	// Absolute SVs
	for (int sv_idx = 0; sv_idx < elem_num_sv; sv_idx++) {
		if (elem_sv_samples[sv_idx].empty()) {
			cout << "  Elem_SV_" << (sv_idx + 1) << ": No samples" << endl;
			continue;
		}

		char name[64];
		sprintf(name, "Elem_SV_%d", sv_idx + 1);

		STASTICS sv_stat;
		sv_stat.setup(elem_sv_samples[sv_idx].data(),
		              elem_sv_samples[sv_idx].size(), name);
		sv_stat.get_cdf();
		sv_stat.print_cdf();

		cout << "  " << name << "_cdf.dat: " << elem_sv_samples[sv_idx].size()
		     << " samples" << endl;
	}

	// SV ratios
	if (elem_num_sv > 1) {
		cout << "Computing Element-Level SV Ratio CDFs..." << endl;
		for (int r = 0; r < elem_num_sv - 1; r++) {
			if (elem_sv_ratio_samples[r].empty()) continue;

			char name[64];
			sprintf(name, "Elem_SV_Ratio_%d_1_dB", r + 2);

			STASTICS ratio_stat;
			ratio_stat.setup(elem_sv_ratio_samples[r].data(),
			                 elem_sv_ratio_samples[r].size(), name);
			ratio_stat.get_cdf();
			ratio_stat.print_cdf();

			cout << "  " << name << "_cdf.dat: " << elem_sv_ratio_samples[r].size()
			     << " samples" << endl;
		}
	}
}

