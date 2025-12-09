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


