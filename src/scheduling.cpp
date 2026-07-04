#include "common.h"
#include "nr_l2sm.h"
#include <numeric>
#include <algorithm>
#include <random>
#include <Eigen/LU>

/*===================================================================
FUNCTION: int getIndex(vector<int> v, int K)
===================================================================*/
int getIndex(vector<int> v, int K)
{
    auto it = find(v.begin(), v.end(), K);
    // If element was found
    if (it != v.end()) 
    {
        // calculating the index
        // of K
        int index = it - v.begin();
        return index;
    }
    else {
        // If the element is not
        // present in the vector
        return -1;
    }
}

/*===================================================================
FUNCTION: Scheduling( int )
===================================================================*/
void Scheduling( void )
{
	// Determine sector range based on single_cell_mode
	// In single_cell_mode, only schedule center BS (BS 0) = sectors 0, 1, 2
	int start_sector = 0;
	int end_sector = (single_cell_mode == 1) ? 3 : (num_SECTORS + num_mTRP_SECTORS);

	if (Scheduling_Type == 0)
	{
		for (int sec_idx = start_sector; sec_idx < end_sector; sec_idx++)
			sector[sec_idx].Scheduling_RoundRobin();
	}
	else if (Scheduling_Type == 1)
	{
		// Parallel processing across sectors
		#if ENABLE_MULTITHREADING
		#pragma omp parallel num_threads(num_of_threads)
		{
			#pragma omp for
		#endif
			for (int sec_idx = start_sector; sec_idx < end_sector; sec_idx++)
			{
				sector[sec_idx].Scheduling_PF_MU_MIMO();
			}
		#if ENABLE_MULTITHREADING
		}
		#endif
	}
	else if (Scheduling_Type == 2 || Scheduling_Type == 3)
	{
		#if ENABLE_MULTITHREADING
		#pragma omp parallel num_threads(num_of_threads)
		{
			#pragma omp for
		#endif
			for (int sec_idx = start_sector; sec_idx < end_sector; sec_idx++)
			{
				// TDD mode: Read full channel matrices, FDD mode: Read PMI vectors
				if (TDD_mode == 1) {
					sector[sec_idx].Read_CSI_Feedback();
				} else {
					sector[sec_idx].Read_Ch_Feedback();
				}
				sector[sec_idx].Initialize_schedule_map();
			}
		#if ENABLE_MULTITHREADING
		}
		#endif

		if ( sector_indices_schedule.size() == 0) {
			vector <int> lc_sector_indices_schedule(end_sector - start_sector);

			iota   (lc_sector_indices_schedule.begin(), lc_sector_indices_schedule.end(), start_sector);
			shuffle(lc_sector_indices_schedule.begin(), lc_sector_indices_schedule.end(), e);

			sector_indices_schedule = lc_sector_indices_schedule;
		}


		//for (int idx = 0; idx < sector_indices_schedule.size(); idx++)
		for (int idx = 0; idx < sector_indices_schedule.size(); idx++)
		{
			int sec_idx = sector_indices_schedule[idx];
			if( sec_idx < num_SECTORS)
			{
				if (sector[sec_idx].schedule_done_flag == 0 && sector[sec_idx].ue_in_control.size() > 0)
				{
					sector[sec_idx].Scheduling_NCJT( UE_NOT_SPECIFIED );
					SCHEDULE_DECISION no_scheduled{UE_NOT_SPECIFIED,UE_NOT_SPECIFIED,UE_NOT_SPECIFIED,UE_NOT_SPECIFIED};
					sector[sec_idx].Set_AVR_Cqi_NCJT( no_scheduled );
				}
			}
			else
			{
				int mTRP_idx = (int)(sec_idx/3);
				if (sector[3*mTRP_idx  ].schedule_done_flag == 0 &&
				    sector[3*mTRP_idx+1].schedule_done_flag == 0 &&
					sector[3*mTRP_idx+2].schedule_done_flag == 0 &&
				    sector[sec_idx].ue_in_control.size() > 0 )
					{
						sector[sec_idx].Scheduling_NCJT( UE_NOT_SPECIFIED );
						SCHEDULE_DECISION no_scheduled{UE_NOT_SPECIFIED,UE_NOT_SPECIFIED,UE_NOT_SPECIFIED,UE_NOT_SPECIFIED};
						sector[sec_idx].Set_AVR_Cqi_NCJT( no_scheduled );
					}				
			}
		}

		for (int sec_idx = start_sector; sec_idx < end_sector; sec_idx++)
			sector[sec_idx].schedule_done_flag = false;
	}
	sector_indices_schedule.clear();
}

/*===================================================================
FUNCTION: void Sector::Scheduling_NCJT( int ms_idx )
===================================================================*/
void Sector::Scheduling_NCJT( int ms_idx )
{	
	if ( ms_idx == UE_NOT_SPECIFIED ) // Initial Scheduling
	{
		Real max_metric = 0;
		int    max_ue     = 0;

		int subband_start_rb_idx = 0; // All rbs are subband
		if ( ppschedulewrite[subband_start_rb_idx][0].ue_selected == -1 ) // if none of resource are allocated yet
		{
			for(int ue_idx = 0; ue_idx < ue_in_control.size(); ue_idx++)
			{
				Real avr_metric = 0;
				for(int rb_idx = 0; rb_idx < num_rb; rb_idx++)
				{
					avr_metric += METRIC[ue_idx][rb_idx]/num_rb; // average all metric for interested resources
				}

				if ( ue_idx == 0 || avr_metric > max_metric)
				{
					max_metric = avr_metric;
					max_ue     = ue_idx;
				}
			}

			for(int rb_idx = 0; rb_idx < num_rb; rb_idx++)
			{
				if ( ppschedulewrite[rb_idx][0].ue_selected == -1 )
				{
					ppschedulewrite[rb_idx][0].ue_selected  = ue_in_control[max_ue];
					ppschedulewrite[rb_idx][0].mcs_selected = MCS_decision[max_ue][rb_idx];
					SU_MIMO_Precoding( max_ue , rb_idx );
				}
			}

		}
	}
	else
	{	
		if (!schedule_done_flag)
		{
			int comp_ue_idx = getIndex( ue_in_comp , ms_idx );
			if ( comp_ue_idx == -1 )
			{
				cout << "Something Wrong with mTRP scheduling : scheduling.cpp 132" << endl;
				getchar();
			}
			else 
			{		
				for(int rb_idx = 0; rb_idx < num_rb; rb_idx++)
				{
					ppschedulewrite[rb_idx][0].ue_selected  = ue_in_comp[comp_ue_idx];
					ppschedulewrite[rb_idx][0].mcs_selected = MCS_decision[ue_in_control.size() + comp_ue_idx][rb_idx];
					SU_MIMO_Precoding( ue_in_control.size() + comp_ue_idx , rb_idx );
				}
			}
		}
	}
}

void Sector::SU_MIMO_Precoding( int scheduled_ue, int rb_idx )
{
	// Use dominant layer for Scheduling_Type=2 (mTRP NCJT). For rank > 1 a
	// multi-layer variant could use all R columns of PMI_vector_read.
	W[rb_idx] = PMI_vector_read[scheduled_ue][rb_idx].col(0);
}


void Sector::Set_AVR_Cqi_NCJT( SCHEDULE_DECISION comp_schedule_info )
{
	if ( comp_schedule_info.ue_selected == UE_NOT_SPECIFIED )
	{
		for (int ue_idx = 0; ue_idx < (int)ue_in_control.size(); ue_idx++)
		{
			int count          = 0;
			Real avr_SINR    = 0.;
			//int comp_sec_idx   = links[ue_in_control[ue_idx]].adj_sector[1];
			int comp_sec_idx   = links[ue_in_control[ue_idx]].comp_sector_idx;
			int comp_mode_flag = 0;

			int comp_ue_idx_in_adjsec = sector[comp_sec_idx].ue_in_control.size() 
			            +getIndex( sector[comp_sec_idx].ue_in_comp , ue_in_control[ue_idx] );

			int cqi_test_result;
			int determined_cqi;
			int determined_mcs;
			Real avr_capacity;

			for (int freq_idx = 0; freq_idx < num_rb; freq_idx++)
			{
				for (int stream = 0; stream < mx_ue_mumimo; stream++)
				{
					if (ppschedulewrite[freq_idx][stream].ue_selected == ue_in_control[ue_idx])
					{
						if (Scheduling_Type == 2 ) 
						{							
							if (comp_mode_flag == 1)
							{
								avr_SINR += CQI_comp_read[ue_idx][freq_idx];
								count++;

								Real ** adj_sector_cqi_read = sector[comp_sec_idx].CQI_read;
								avr_SINR += adj_sector_cqi_read[comp_ue_idx_in_adjsec][freq_idx];
								count++;
							}
							else
							{
								avr_SINR += CQI_read[ue_idx][freq_idx];
								count++;
							}
						}
						else if (Scheduling_Type == 3 )
						{
							if (comp_mode_flag == 1)
							{							
								Real ** adj_sector_cqi_read = sector[comp_sec_idx].CQI_read;
								avr_SINR += adj_sector_cqi_read[comp_ue_idx_in_adjsec][freq_idx];
								count++;
							}
							else
							{
								avr_SINR += CQI_read[ue_idx][freq_idx];
								count++;								
							}
						}
					}
				}
			}

			avr_SINR     = 10. * log10(avr_SINR/count);
			int temp_cqi = determine_CQI(avr_SINR)-1;// + 3;
			count          = 0;
			avr_SINR       = 0.;

			for (int freq_idx = 0; freq_idx < num_rb; freq_idx++)
			{
				for (int stream = 0; stream < mx_ue_mumimo; stream++)
				{
					if (ppschedulewrite[freq_idx][stream].ue_selected == ue_in_control[ue_idx])
					{
						if (Scheduling_Type == 2 )
						{
							if ( comp_mode_flag )
							{
								avr_SINR += exp(-CQI_comp_read[ue_idx][freq_idx] / Beta[temp_cqi]);
								count++;

								Real ** adj_sector_cqi_read = sector[comp_sec_idx].CQI_read;
								avr_SINR += exp(-adj_sector_cqi_read[comp_ue_idx_in_adjsec][freq_idx] / Beta[temp_cqi]);
								count++;
							}
							else
							{
								avr_SINR += exp(-CQI_read[ue_idx][freq_idx] / Beta[temp_cqi]);
								count++;
							}
						} 
						else if (Scheduling_Type == 3)
						{
							if ( comp_mode_flag )
							{
								Real ** adj_sector_cqi_read = sector[comp_sec_idx].CQI_read;
								avr_SINR += exp(-adj_sector_cqi_read[comp_ue_idx_in_adjsec][freq_idx] / Beta[temp_cqi]);
								count++;								
							}
							else
							{
								avr_SINR += exp(-CQI_read[ue_idx][freq_idx] / Beta[temp_cqi]);
								count++;
							}
						}
					}
				}
			}
				
			avr_SINR           = -Beta[temp_cqi] * log(avr_SINR / (Real)(count));	
			Real avr_SINR_dB = 10. * log10(avr_SINR);
			cqi_test_result    = determine_CQI(avr_SINR_dB);// + 3;

			while ( temp_cqi > cqi_test_result-1 ) 
			{
				avr_SINR = 0;
				avr_SINR_dB = 0;
				count = 0;
				temp_cqi = cqi_test_result-1;
				for (int freq_idx = 0; freq_idx < num_rb; freq_idx++)
				{
					for (int stream = 0; stream < mx_ue_mumimo; stream++)
					{
						if (ppschedulewrite[freq_idx][stream].ue_selected == ue_in_control[ue_idx])
						{
							if (Scheduling_Type == 2 )
							{
								if ( comp_mode_flag )
								{
									avr_SINR += exp(-CQI_comp_read[ue_idx][freq_idx] / Beta[temp_cqi]);
									count++;

									Real ** adj_sector_cqi_read = sector[comp_sec_idx].CQI_read;
									avr_SINR += exp(-adj_sector_cqi_read[comp_ue_idx_in_adjsec][freq_idx] / Beta[temp_cqi]);
									count++;
								}
								else
								{
									avr_SINR += exp(-CQI_read[ue_idx][freq_idx] / Beta[temp_cqi]);
									count++;
								}
							} 
							else if (Scheduling_Type == 3)
							{
								if ( comp_mode_flag )
								{
									Real ** adj_sector_cqi_read = sector[comp_sec_idx].CQI_read;
									avr_SINR += exp(-adj_sector_cqi_read[comp_ue_idx_in_adjsec][freq_idx] / Beta[temp_cqi]);
									count++;								
								}
								else
								{
									avr_SINR += exp(-CQI_read[ue_idx][freq_idx] / Beta[temp_cqi]);
									count++;
								}
							}
						}
					}
				}
				avr_SINR        = -Beta[temp_cqi] * log(avr_SINR / (Real)(count));	
				avr_SINR_dB     = 10. * log10(avr_SINR);
				cqi_test_result = determine_CQI(avr_SINR_dB);// + 3;
			}
			determined_cqi = temp_cqi + 1;

			if ( ms[ue_in_control[ue_idx]].re_tx_failed_cnt > 0)
			{
				ms[ue_in_control[ue_idx]].cqi_offset = ms[ue_in_control[ue_idx]].cqi_offset - 1;
				ms[ue_in_control[ue_idx]].re_tx_failed_cnt = 0;
			}
			determined_cqi = determined_cqi + ms[ue_in_control[ue_idx]].cqi_offset;
			
			if ( determined_cqi < 1)
				determined_cqi = 1;

			if      (determined_cqi == 1)   { determined_mcs = 0;  }
			else if (determined_cqi == 2)   { determined_mcs = 1;  } //0
			else if (determined_cqi == 3)   { determined_mcs = 2;  }
			else if (determined_cqi == 4)   { determined_mcs = 4;  }
			else if (determined_cqi == 5)   { determined_mcs = 6;  }
			else if (determined_cqi == 6)   { determined_mcs = 8;  }
			else if (determined_cqi == 7)   { determined_mcs = 11; }
			else if (determined_cqi == 8)   { determined_mcs = 13; }
			else if (determined_cqi == 9)   { determined_mcs = 16; }
			else if (determined_cqi == 10)  { determined_mcs = 18; }
			else if (determined_cqi == 11)  { determined_mcs = 21; }
			else if (determined_cqi == 12)  { determined_mcs = 23; }
			else if (determined_cqi == 13)  { determined_mcs = 25; }
			else if (determined_cqi == 14)  { determined_mcs = 26; } //27
			else if (determined_cqi == 15)  { determined_mcs = 27; }
			

			if (determined_mcs > 27)
			{
				determined_mcs = 27;
			}
			else if (determined_mcs < 0)
			{
				determined_mcs = 0;
			}

			for (int freq_idx = 0; freq_idx < num_rb; freq_idx++)
			{
				for (int stream = 0; stream < mx_ue_mumimo; stream++)
				{
					if (ppschedulewrite[freq_idx][stream].ue_selected == ue_in_control[ue_idx])
					{
						ppschedulewrite[freq_idx][stream].mcs_selected = determined_mcs;
						ppschedulewrite[freq_idx][stream].cqi_selected = determined_cqi;
						ppschedulewrite[freq_idx][stream].capacity     = avr_SINR;
						ppschedulewrite[freq_idx][stream].temp_cqi     = temp_cqi;
					}
				}
			}

			if (comp_mode_flag 
			    && ppschedulewrite[0][0].ue_selected == ue_in_control[ue_idx] )
				//&& sector[comp_sec_idx].ppschedulewrite[0][0].ue_selected == ue_in_control[ue_idx])
			{
				SCHEDULE_DECISION comp_info{ue_in_control[ue_idx],determined_mcs,determined_cqi,avr_capacity};
				sector[comp_sec_idx].Set_AVR_Cqi_NCJT( comp_info );
			}
		}
		schedule_done_flag = true;	
	} 
	else
	{
		int comp_ue_idx = getIndex( ue_in_comp , comp_schedule_info.ue_selected );
		for (int freq_idx = 0; freq_idx < num_rb; freq_idx++)
		{
			for (int stream = 0; stream < mx_ue_mumimo; stream++)
			{
				if (ppschedulewrite[freq_idx][stream].ue_selected == comp_schedule_info.ue_selected )
				{
					ppschedulewrite[freq_idx][stream].mcs_selected = comp_schedule_info.mcs_selected;
					ppschedulewrite[freq_idx][stream].cqi_selected = comp_schedule_info.cqi_selected;
					ppschedulewrite[freq_idx][stream].cqi_selected = comp_schedule_info.capacity;
				}
			}
		}
		schedule_done_flag = true;	
	}
}

/*===================================================================
FUNCTION: Sector::Scheduling_RoundRobin( int )
===================================================================*/
void Sector::Scheduling_RoundRobin(void)
{
	/*===================================================================
	FUNCTION : Read_Ch_Feedback()
	INPUT    : pppPMI_map,pppCQI_map
	OUTPUT   : PMI_read, CQI_read, MCS_decision,CQI_AVR,METRIC
	===================================================================*/
	Read_Ch_Feedback();

	/*===================================================================
	FUNCTION : Initialize_schedule_map()
	INPUT    : ue_selected, mcs_selected, cqi_selected
	OUTPUT   : none
	===================================================================*/
	Initialize_schedule_map();

	/*===================================================================
	FUNCTION : Scheduling_algorithm_module_RoundRobin();
	INPUT    : ue_In_control
	OUTPUT   : ue_selected
	===================================================================*/
	Scheduling_algorithm_module_RoundRobin();

	/*===================================================================
	FUNCTION : Set_AVR_Cqi() or Set_AVR_Cqi_Precoding_Based()
	INPUT    : ue_selected, precoding matrix W (if USE_PRECODING_BASED_SINR=1)
	OUTPUT   : mcs_selected, cqi_selected
	===================================================================*/
	if (USE_PRECODING_BASED_SINR == 1)
	{
		// Calculate SINR from actual precoding matrix and channel
		Set_AVR_Cqi_Precoding_Based();
	}
	else
	{
		// Use feedback CQI (default)
		Set_AVR_Cqi();
	}

}

#if 0  // unused: Scheduling_mTRP — never called
void Sector::Scheduling_mTRP( int ms_idx )
{
	//Read_Ch_Feedback();
	//Initialize_schedule_map();
	Scheduling_NCJT( ms_idx );
	Set_AVR_Cqi();
	AVR_Cqi_update();

	schedule_done_flag = true;
}
#endif  // unused: Scheduling_mTRP



void Sector::Scheduling_PF_MU_MIMO(void)
{
	/*===================================================================
	FUNCTION : Read_Ch_Feedback() or Read_CSI_Feedback()
	INPUT    : pppPMI_map,pppCQI_map (FDD) or ppppCSI_matrix_map (TDD)
	OUTPUT   : PMI_read, CQI_read, MCS_decision,CQI_AVR,METRIC
	===================================================================*/
	// TDD mode: Read full channel matrices, FDD mode: Read PMI vectors
	if (TDD_mode == 1) {
		Read_CSI_Feedback();
	} else {
		Read_Ch_Feedback();
	}

	/*===================================================================
	FUNCTION : Initialize_schedule_map()
	INPUT    : ue_selected, mcs_selected, cqi_selected
	OUTPUT   : none
	===================================================================*/
	Initialize_schedule_map();

	/*===================================================================
	FUNCTION : Scheduling_algorithm_module_MU_MIMO();
	INPUT    : ue_In_control
	OUTPUT   : ue_selected
	===================================================================*/
	Scheduling_algorithm_module_MU_MIMO();
		
	/*===================================================================
	FUNCTION : Set_AVR_Cqi() or Set_AVR_Cqi_Precoding_Based()
	INPUT    : ue_selected, precoding matrix W (if USE_PRECODING_BASED_SINR=1)
	OUTPUT   : mcs_selected, cqi_selected
	===================================================================*/
	if (USE_PRECODING_BASED_SINR == 1)
	{
		// Calculate SINR from actual precoding matrix and channel
		Set_AVR_Cqi_Precoding_Based();
	}
	else
	{
		// Use feedback CQI (default)
		Set_AVR_Cqi();
	}

	/*===================================================================
	FUNCTION : AVR_Cqi_update()
	INPUT    : ue_selected
	OUTPUT   : CQI_AVR
	===================================================================*/
	AVR_Cqi_update(); 
}

void Sector::Read_Ch_Feedback(void)
{

	if (t < N_pf + SCHEDULE_DELAY)
	{
		_N_pf = (t + 1) - SCHEDULE_DELAY;
	}
	else
	{
		_N_pf = N_pf;
	}

	Real alpha = 1;

	const int fixed_mem = 80;
	if (PMI_read == NULL)
	{
		PMI_read         = new PMI_FEEDBACK  * [fixed_mem];  //[ue_in_control.size() + ue_in_comp.size()]; //30
		CQI_read         = new Real       * [fixed_mem];  //[ue_in_control.size() + ue_in_comp.size()];
		CQI_comp_read    = new Real       * [fixed_mem];  //[ue_in_control.size() + ue_in_comp.size()];
		MCS_decision     = new int          * [fixed_mem];  //[ue_in_control.size() + ue_in_comp.size()];
		CQIndex_decision = new int          * [fixed_mem];  //[ue_in_control.size() + ue_in_comp.size()];
		CQI_AVR          = new Real       * [fixed_mem];  //[ue_in_control.size() + ue_in_comp.size()];
		METRIC           = new Real       * [fixed_mem];  //[ue_in_control.size() + ue_in_comp.size()];
		PMI_vector_read  = new MatrixXcReal    * [fixed_mem];  //[ue_in_control.size() + ue_in_comp.size()];

		//for (int ue_idx = 0; ue_idx < (int)ue_in_control.size() + ue_in_comp.size(); ue_idx++)
		for (int ue_idx = 0; ue_idx < fixed_mem; ue_idx++)
		{
			PMI_read[ue_idx]         = new PMI_FEEDBACK [num_rb];
			CQI_read[ue_idx]         = new    Real [num_rb];
			CQI_comp_read[ue_idx]    = new    Real [num_rb];
			MCS_decision[ue_idx]     = new       int [num_rb];
			CQIndex_decision[ue_idx] = new       int [num_rb];
			CQI_AVR[ue_idx]          = new    Real [num_rb];
			METRIC[ue_idx]           = new    Real [num_rb];
			PMI_vector_read[ue_idx]  = new MatrixXcReal [num_rb];
		}

		// Initialize all allocated memory (fixed_mem, not just 40!)
		for (int ue_idx = 0; ue_idx < fixed_mem; ue_idx++)
		{
			for (int rb_idx = 0; rb_idx < num_rb; rb_idx++)
			{
				PMI_read        [ue_idx][rb_idx] = {0, 0, 0};
				CQI_read        [ue_idx][rb_idx] = 0.;
				CQI_comp_read   [ue_idx][rb_idx] = 0.;
				MCS_decision    [ue_idx][rb_idx] = 0;
				CQIndex_decision[ue_idx][rb_idx] = 0;
				CQI_AVR         [ue_idx][rb_idx] = 0;
				METRIC          [ue_idx][rb_idx] = 0;
				PMI_vector_read [ue_idx][rb_idx] = MatrixXcReal::Zero(NUM_TX_Port, g_type2_rank);
			}
		}
	}

	if ( ue_in_control.size() + ue_in_comp.size() > fixed_mem ) {
		cout << " Needs more memory !!!" <<" sec_idx = "<< self_idx<<"  ue_in_control.size() = "<< ue_in_control.size() <<"\t"<< "ue_in_comp.size() = "<< ue_in_comp.size() << endl;
		getchar();
	}

	if (t == 0) {
		//for (int ue_idx = 0; ue_idx < (int)ue_in_control.size() + ue_in_comp.size(); ue_idx++)
		for (int ue_idx = 0; ue_idx < 40 ; ue_idx++)
		{
			for (int rb_idx = 0; rb_idx < num_rb; rb_idx++)
			{
				PMI_read        [ue_idx][rb_idx] = {0, 0, 0};
				CQI_read        [ue_idx][rb_idx] = 0.;
				CQI_comp_read   [ue_idx][rb_idx] = 0.;
				MCS_decision    [ue_idx][rb_idx] = 0;
				CQIndex_decision[ue_idx][rb_idx] = 0;
				CQI_AVR         [ue_idx][rb_idx] = 0;
				METRIC          [ue_idx][rb_idx] = 0;
				PMI_vector_read [ue_idx][rb_idx] = MatrixXcReal::Zero(NUM_TX_Port, g_type2_rank);
			}
		}
	}

	for (int ue_idx = 0; ue_idx < (int)ue_in_control.size() + ue_in_comp.size(); ue_idx++)
	{
		int ms_idx;
		if (ue_idx < ue_in_control.size())
			ms_idx = ue_in_control[ue_idx];
		else
			ms_idx = ue_in_comp[ue_idx-ue_in_control.size()];

		for (int rb_idx = 0; rb_idx < num_rb; rb_idx++)
		{
			int adj_sec_idx = 0;
			if (ue_idx  >= ue_in_control.size())
				adj_sec_idx = 1;

			PMI_read       [ue_idx][rb_idx] = ppppPMI_map       [ms_idx][adj_sec_idx][rb_idx][(cqi_history_length + t - SCHEDULE_DELAY) % cqi_history_length];
			CQI_read       [ue_idx][rb_idx] = ppppCQI_Map       [ms_idx][adj_sec_idx][rb_idx][(cqi_history_length + t - SCHEDULE_DELAY) % cqi_history_length]; 
			CQI_comp_read  [ue_idx][rb_idx] = ppppCQI_comp_Map  [ms_idx][adj_sec_idx][rb_idx][(cqi_history_length + t - SCHEDULE_DELAY) % cqi_history_length]; 
			PMI_vector_read[ue_idx][rb_idx] = ppppPMI_vector_map[ms_idx][adj_sec_idx][rb_idx][(cqi_history_length + t - SCHEDULE_DELAY) % cqi_history_length];

			Real cqi_offset_added = CQI_read[ue_idx][rb_idx];
			int mcs_up = determine_MCS(10*log10(cqi_offset_added));//+ mcs_offset; /// 
			//int cqi_up = determine_CQI((cqi_offset_added));//+ mcs_offset; /// 
			int cqi_up = mcs_up;

			if (true)  // Use MCS index
			{
				if (mcs_up > 27)
				{
					MCS_decision[ue_idx][rb_idx] = 27;
				}
				else if (mcs_up <0)
				{
					MCS_decision[ue_idx][rb_idx] = 0;
				}
				else
				{
					MCS_decision[ue_idx][rb_idx] = mcs_up;
				}

			}
			else  ///// Use CQI index 
			{
				// save CQI index
				MCS_decision[ue_idx][rb_idx] = mcs_up;

			}

			if (cqi_up > 15)
			{
				CQIndex_decision[ue_idx][rb_idx] = 15;
			}
			else if (cqi_up < 1)
			{
				CQIndex_decision[ue_idx][rb_idx] = 1;
			}
			else
			{
				CQIndex_decision[ue_idx][rb_idx] = cqi_up;
			}
		}
	}

	//if ( t > 100)
	//{
	//	cout << "Scheduling Check !!" << endl;
	//}

	for (int ue_idx = 0; ue_idx < (int)ue_in_control.size() + ue_in_comp.size(); ue_idx++) 
	{				
		for (int rb_idx = 0; rb_idx < num_rb; rb_idx++)
		{
			if (CQI_AVR[ue_idx][rb_idx] > 0 )
			{
				if( ue_idx < (int)ue_in_control.size() ) {
					if( ppschedulewrite[rb_idx][0].ue_selected == ue_in_control[ue_idx] ) 
					{
						Real current_capacity_on_rb = log2(1 + CQI_read[ue_idx][rb_idx]);
						int comp_sec_idx = links[ue_in_control[ue_idx]].comp_sector_idx;
						if ( Scheduling_Type == 2 && sector[comp_sec_idx].ppschedulewrite[rb_idx][0].ue_selected == ue_in_control[ue_idx] )
						{	
							int adjsector_comp_ue_idx = sector[comp_sec_idx].ue_in_control.size() + 
														getIndex( sector[comp_sec_idx].ue_in_comp , ue_in_control[ue_idx] );
						
							current_capacity_on_rb  = log2(1 + CQI_comp_read[ue_idx][rb_idx]);
							current_capacity_on_rb += log2(1 + sector[comp_sec_idx].CQI_read[adjsector_comp_ue_idx][rb_idx]);
						}

						CQI_AVR[ue_idx][rb_idx] = (1. / _N_pf) * current_capacity_on_rb 
											+ (1 - 1. / _N_pf) * CQI_AVR[ue_idx][rb_idx];
					}
					else
					{
						CQI_AVR[ue_idx][rb_idx] =  (1 - 1. / _N_pf) * CQI_AVR[ue_idx][rb_idx];
					}						
				}
				else
				{
					if( ppschedulewrite[rb_idx][0].ue_selected == ue_in_comp[ue_idx-ue_in_control.size()] )  
					{
						CQI_AVR[ue_idx][rb_idx] = (1. / _N_pf) * log2(1 + CQI_read[ue_idx][rb_idx]) 
											+ (1 - 1. / _N_pf) * CQI_AVR[ue_idx][rb_idx];
					}
					else
					{
						CQI_AVR[ue_idx][rb_idx] =  (1 - 1. / _N_pf) * CQI_AVR[ue_idx][rb_idx];
					}	
				}					
			}
			else  ///// first loop 
			{
				CQI_AVR[ue_idx][rb_idx] = log2(1 + CQI_read[ue_idx][rb_idx]);
			}
		}
	}

	if (Scheduling_Type == 0)
	{
		for (int ue_idx = 0; ue_idx < (int)ue_in_control.size() + ue_in_comp.size(); ue_idx++)
		{
			for (int rb_idx = 0; rb_idx < num_rb; rb_idx++)
			{
				METRIC[ue_idx][rb_idx] = log2(1 + CQI_read[ue_idx][rb_idx]);//
			}
		}
	}
	else if (Scheduling_Type == 1 || Scheduling_Type == 2 || Scheduling_Type == 3 )  // 
	{
		for (int ue_idx = 0; ue_idx < (int)ue_in_control.size() + ue_in_comp.size(); ue_idx++)
		{
			Real alpha = 1;
			Real beta  = 1;
			//if ( t < 500 )
			//alpha = 0;

			for (int rb_idx = 0; rb_idx < num_rb; rb_idx++)
			{
				Real instantaneous_rate = log2(1 + CQI_read[ue_idx][rb_idx]);
				Real average_rate = CQI_AVR[ue_idx][rb_idx];

				// Avoid division by same value (happens on first scheduling)
				// When instantaneous == average, all UEs get METRIC=1
				if (std::abs(instantaneous_rate - average_rate) < 1e-6) {
					// First scheduling: use instantaneous rate only
					METRIC[ue_idx][rb_idx] = pow(instantaneous_rate, alpha);
				} else {
					// Normal Proportional Fair metric
					METRIC[ue_idx][rb_idx] = pow(instantaneous_rate, alpha)
					                        / pow(average_rate, beta);
				}
			}
		}
	}
}


int Sector::determine_MCS(Real sinr_estimate)
{
	// Throughput-maximizing selection over all 28 MCS (TBS/segmentation aware);
	// see the free determine_MCS in receive_downlink.cpp. Both copies hooked.
	if (g_matlab_bler && g_matlab_tput_mcs && TputMCS_Grid_ready())
		return nr_mcs_maxtput_from_dB(sinr_estimate);

	if (false)  //// TS 38.214 , Table 5.1.3.1-2:, use spectral efficiency
	{
		int _MCS_decision = 0;


		if (sinr_estimate < 0.2344)       { _MCS_decision = 0; }
		else if (sinr_estimate < 0.3770)  { _MCS_decision = 1; }
		else if (sinr_estimate < 0.6016)  { _MCS_decision = 2; }
		else if (sinr_estimate < 0.8770)  { _MCS_decision = 3; }
		else if (sinr_estimate < 1.1758)  { _MCS_decision = 4; }
		else if (sinr_estimate < 1.4766)  { _MCS_decision = 5; }
		else if (sinr_estimate < 1.6953)  { _MCS_decision = 6; }
		else if (sinr_estimate < 1.9141)  { _MCS_decision = 7; }
		else if (sinr_estimate < 2.1602)  { _MCS_decision = 8; }
		else if (sinr_estimate < 2.4063)  { _MCS_decision = 9; }
		else if (sinr_estimate < 2.5703)  { _MCS_decision = 10; }
		else if (sinr_estimate < 2.7305)  { _MCS_decision = 11; }
		else if (sinr_estimate < 3.0293)  { _MCS_decision = 12; }
		else if (sinr_estimate < 3.3223)  { _MCS_decision = 13; }
		else if (sinr_estimate < 3.6094)  { _MCS_decision = 14; }
		else if (sinr_estimate < 3.9023)  { _MCS_decision = 15; }
		else if (sinr_estimate < 4.2129)  { _MCS_decision = 16; }
		else if (sinr_estimate < 4.5234)  { _MCS_decision = 17; }
		else if (sinr_estimate < 4.8164)  { _MCS_decision = 18; }
		else if (sinr_estimate < 5.1152)  { _MCS_decision = 19; }
		else if (sinr_estimate < 5.3320)  { _MCS_decision = 20; }
		else if (sinr_estimate < 5.5547)  { _MCS_decision = 21; }
		else if (sinr_estimate < 5.8906)  { _MCS_decision = 22; }
		else if (sinr_estimate < 6.2266)  { _MCS_decision = 23; }
		else if (sinr_estimate < 6.5703)  { _MCS_decision = 24; }
		else if (sinr_estimate < 6.9141)  { _MCS_decision = 25; }
		else if (sinr_estimate < 7.1602)  { _MCS_decision = 26; }
		else if (sinr_estimate < 7.4063)  { _MCS_decision = 27; }
		else                              { _MCS_decision = 27; }


		return _MCS_decision;
	}

	else if (false)  /// CQI
	{
		int _CQI_decision = 0;


		if (sinr_estimate < -9.478)  { _CQI_decision = 1; }
		else if (sinr_estimate < -6.658)  { _CQI_decision = 2; }
		else if (sinr_estimate < -4.098)  { _CQI_decision = 3; }
		else if (sinr_estimate < -1.798)  { _CQI_decision = 4; }
		else if (sinr_estimate < 0.399)   { _CQI_decision = 5; }
		else if (sinr_estimate < 2.424)   { _CQI_decision = 6; }
		else if (sinr_estimate < 4.489)   { _CQI_decision = 7; }
		else if (sinr_estimate < 6.367)   { _CQI_decision = 8; }
		else if (sinr_estimate < 8.456)   { _CQI_decision = 9; }
		else if (sinr_estimate < 10.266)  { _CQI_decision = 10; }
		else if (sinr_estimate < 12.218)  { _CQI_decision = 11; }
		else if (sinr_estimate < 14.122)  { _CQI_decision = 12; }
		else if (sinr_estimate < 15.849)  { _CQI_decision = 13; }
		else if (sinr_estimate < 17.786)  { _CQI_decision = 14; }
		else                              { _CQI_decision = 15; }


		return _CQI_decision;

	}

	else if (true)  
	{
		int _CQI_decision = 0;

		int _MCS_decision = 0;
			/// VIENNA curve

		if (sinr_estimate > SINR_threshold_dB[14])       { _CQI_decision = 15; }
		else if (sinr_estimate > SINR_threshold_dB[13])  { _CQI_decision = 14; }
		else if (sinr_estimate > SINR_threshold_dB[12])  { _CQI_decision = 13; }
		else if (sinr_estimate > SINR_threshold_dB[11])  { _CQI_decision = 12; }
		else if (sinr_estimate > SINR_threshold_dB[10])  { _CQI_decision = 11; }
		else if (sinr_estimate > SINR_threshold_dB[9])  { _CQI_decision = 10; }
		else if (sinr_estimate > SINR_threshold_dB[8])  { _CQI_decision = 9; }
		else if (sinr_estimate > SINR_threshold_dB[7])  { _CQI_decision = 8; }
		else if (sinr_estimate > SINR_threshold_dB[6])  { _CQI_decision = 7; }
		else if (sinr_estimate > SINR_threshold_dB[5])  { _CQI_decision = 6; }
		else if (sinr_estimate > SINR_threshold_dB[4])  { _CQI_decision = 5; }
		else if (sinr_estimate > SINR_threshold_dB[3])  { _CQI_decision = 4; }
		else if (sinr_estimate > SINR_threshold_dB[2])  { _CQI_decision = 3; }
		else if (sinr_estimate > SINR_threshold_dB[1])  { _CQI_decision = 2; }
		else if (sinr_estimate > SINR_threshold_dB[0])  { _CQI_decision = 1; }
		else                                            { _CQI_decision = 1; }

			//_CQI_decision = 8;
			// TS 36.101  TABLE A.4-7 �� MCS Scheme 1
			if      (_CQI_decision == 1)   { _MCS_decision = 0;  }
			else if (_CQI_decision == 2)   { _MCS_decision = 0;  }
			else if (_CQI_decision == 3)   { _MCS_decision = 2;  }
			else if (_CQI_decision == 4)   { _MCS_decision = 4;  }
			else if (_CQI_decision == 5)   { _MCS_decision = 6;  }
			else if (_CQI_decision == 6)   { _MCS_decision = 8;  }
			else if (_CQI_decision == 7)   { _MCS_decision = 11; }
			else if (_CQI_decision == 8)   { _MCS_decision = 13; }
			else if (_CQI_decision == 9)   { _MCS_decision = 16; }
			else if (_CQI_decision == 10)  { _MCS_decision = 18; }
			else if (_CQI_decision == 11)  { _MCS_decision = 21; }
			else if (_CQI_decision == 12)  { _MCS_decision = 23; }
			else if (_CQI_decision == 13)  { _MCS_decision = 25; }
			else if (_CQI_decision == 14)  { _MCS_decision = 27; }
			else if (_CQI_decision == 15)  { _MCS_decision = 27; }



			return _MCS_decision;

	
	}


	else  //  old code
	{
	}


	
}





int Sector::determine_CQI(Real sinr_estimate)
{
	if (true)  /// CQI
	{
		int _CQI_decision = 0;
		if (sinr_estimate > SINR_threshold_dB[14])       { _CQI_decision = 15; }
		else if (sinr_estimate > SINR_threshold_dB[13])  { _CQI_decision = 14; }
		else if (sinr_estimate > SINR_threshold_dB[12])  { _CQI_decision = 13; }
		else if (sinr_estimate > SINR_threshold_dB[11])  { _CQI_decision = 12; }
		else if (sinr_estimate > SINR_threshold_dB[10])  { _CQI_decision = 11; }
		else if (sinr_estimate > SINR_threshold_dB[9])  { _CQI_decision = 10; }
		else if (sinr_estimate > SINR_threshold_dB[8])  { _CQI_decision = 9; }
		else if (sinr_estimate > SINR_threshold_dB[7])  { _CQI_decision = 8; }
		else if (sinr_estimate > SINR_threshold_dB[6])  { _CQI_decision = 7; }
		else if (sinr_estimate > SINR_threshold_dB[5])  { _CQI_decision = 6; }
		else if (sinr_estimate > SINR_threshold_dB[4])  { _CQI_decision = 5; }
		else if (sinr_estimate > SINR_threshold_dB[3])  { _CQI_decision = 4; }
		else if (sinr_estimate > SINR_threshold_dB[2])  { _CQI_decision = 3; }
		else if (sinr_estimate > SINR_threshold_dB[1])  { _CQI_decision = 2; }
		else if (sinr_estimate > SINR_threshold_dB[0])  { _CQI_decision = 1; }
		else                                            { _CQI_decision = 1; }
		return _CQI_decision;
	}
}


/*===================================================================
FUNCTION: Sector::Read_CSI_Feedback()

DESCRIPTION:
  TDD reciprocity approach: Read full channel matrices from UEs
  instead of quantized PMI vectors. This allows BS to directly
  compute ZF/MMSE precoding using actual channel knowledge.

INPUT:  ms[].CSI_matrix - Stored channel matrices from UEs
OUTPUT: CSI_matrix_read[ue_idx][rb_idx] - Channel matrices for precoding
===================================================================*/
void Sector::Read_CSI_Feedback(void)
{
	const int fixed_mem = 80;

	// Allocate memory on first call
	if (CSI_matrix_read == NULL)
	{
		CSI_matrix_read = new MatrixXcReal * [fixed_mem];
		for (int ue_idx = 0; ue_idx < fixed_mem; ue_idx++)
		{
			CSI_matrix_read[ue_idx] = new MatrixXcReal [num_rb];

			// Initialize each matrix with proper size
			for (int rb_idx = 0; rb_idx < num_rb; rb_idx++)
			{
				CSI_matrix_read[ue_idx][rb_idx] = MatrixXcReal::Zero(NUM_RX_Port, NUM_TX_Port);
			}
		}		
	}

	if (t < N_pf + SCHEDULE_DELAY)
	{
		_N_pf = (t + 1) - SCHEDULE_DELAY;
	}
	else
	{
		_N_pf = N_pf;
	}

	Real alpha = 1;

	if (PMI_read == NULL)
	{
		PMI_read         = new PMI_FEEDBACK  * [fixed_mem];  //[ue_in_control.size() + ue_in_comp.size()]; //30
		CQI_read         = new Real         * [fixed_mem];  //[ue_in_control.size() + ue_in_comp.size()];
		CQI_comp_read    = new Real         * [fixed_mem];  //[ue_in_control.size() + ue_in_comp.size()];
		MCS_decision     = new int          * [fixed_mem];  //[ue_in_control.size() + ue_in_comp.size()];
		CQIndex_decision = new int          * [fixed_mem];  //[ue_in_control.size() + ue_in_comp.size()];
		CQI_AVR          = new Real         * [fixed_mem];  //[ue_in_control.size() + ue_in_comp.size()];
		METRIC           = new Real         * [fixed_mem];  //[ue_in_control.size() + ue_in_comp.size()];
		PMI_vector_read  = new MatrixXcReal    * [fixed_mem];  //[ue_in_control.size() + ue_in_comp.size()];

		//for (int ue_idx = 0; ue_idx < (int)ue_in_control.size() + ue_in_comp.size(); ue_idx++)
		for (int ue_idx = 0; ue_idx < fixed_mem; ue_idx++)
		{
			PMI_read[ue_idx]         = new PMI_FEEDBACK [num_rb];
			CQI_read[ue_idx]         = new      Real [num_rb];
			CQI_comp_read[ue_idx]    = new      Real [num_rb];
			MCS_decision[ue_idx]     = new       int [num_rb];
			CQIndex_decision[ue_idx] = new       int [num_rb];
			CQI_AVR[ue_idx]          = new      Real [num_rb];
			METRIC[ue_idx]           = new      Real [num_rb];
			PMI_vector_read[ue_idx]  = new MatrixXcReal [num_rb];
		}

		// Initialize all allocated memory (fixed_mem, not just 40!)
		for (int ue_idx = 0; ue_idx < fixed_mem; ue_idx++)
		{
			for (int rb_idx = 0; rb_idx < num_rb; rb_idx++)
			{
				PMI_read        [ue_idx][rb_idx] = {0, 0, 0};
				CQI_read        [ue_idx][rb_idx] = 0.;
				CQI_comp_read   [ue_idx][rb_idx] = 0.;
				MCS_decision    [ue_idx][rb_idx] = 0;
				CQIndex_decision[ue_idx][rb_idx] = 0;
				CQI_AVR         [ue_idx][rb_idx] = 0;
				METRIC          [ue_idx][rb_idx] = 0;
				PMI_vector_read [ue_idx][rb_idx] = MatrixXcReal::Zero(NUM_TX_Port, g_type2_rank);
			}
		}
	}

	// Read channel matrices and CQI from UEs (with SCHEDULE_DELAY)
	for (int ue_idx = 0; ue_idx < (int)ue_in_control.size() ; ue_idx++)
	{
		int ms_idx;
		if (ue_idx < (int)(ue_in_control.size()))
			ms_idx = ue_in_control[ue_idx];
		else
			ms_idx = ue_in_comp[ue_idx - ue_in_control.size()];

		for (int rb_idx = 0; rb_idx < num_rb; rb_idx++)
		{
			// Determine which coefficient index to use (serving sector)
			int adj_sec_idx = 0;

			// Read channel matrix with delay (similar to PMI feedback delay)
			// ppppCSI_matrix_map[ms_idx][adj_sec_idx][rb_idx][history_idx]
			CSI_matrix_read[ue_idx][rb_idx] = ppppCSI_matrix_map[ms_idx][adj_sec_idx][rb_idx][(cqi_history_length + t - SCHEDULE_DELAY) % cqi_history_length];

			// Read PMI vector (dominant eigenvector) with delay
			PMI_vector_read[ue_idx][rb_idx] = ppppPMI_vector_map[ms_idx][adj_sec_idx][rb_idx][(cqi_history_length + t - SCHEDULE_DELAY) % cqi_history_length];

			// Read CQI values (same as FDD mode)
			CQI_read[ue_idx][rb_idx] = ppppCQI_Map[ms_idx][adj_sec_idx][rb_idx][(cqi_history_length + t - SCHEDULE_DELAY) % cqi_history_length];

			// Compute MCS and CQI index decision (same logic as Read_Ch_Feedback)
			Real cqi_offset_added = CQI_read[ue_idx][rb_idx];
			int mcs_up = determine_MCS(10*log10(cqi_offset_added));
			int cqi_up = mcs_up;

			if (true)  // Use MCS index
			{
				if (mcs_up > 27)
				{
					MCS_decision[ue_idx][rb_idx] = 27;
				}
				else if (mcs_up < 0)
				{
					MCS_decision[ue_idx][rb_idx] = 0;
				}
				else
				{
					MCS_decision[ue_idx][rb_idx] = mcs_up;
				}
			}
			else  // Use CQI index
			{
				MCS_decision[ue_idx][rb_idx] = mcs_up;
			}

			if (cqi_up > 15)
			{
				CQIndex_decision[ue_idx][rb_idx] = 15;
			}
			else if (cqi_up < 1)
			{
				CQIndex_decision[ue_idx][rb_idx] = 1;
			}
			else
			{
				CQIndex_decision[ue_idx][rb_idx] = cqi_up;
			}
		}
	}

	// Update CQI_AVR for proportional fair scheduling
	for (int ue_idx = 0; ue_idx < (int)ue_in_control.size(); ue_idx++)
	{
		for (int rb_idx = 0; rb_idx < num_rb; rb_idx++)
		{
			if (CQI_AVR[ue_idx][rb_idx] > 0)
			{
				if (ue_idx < (int)ue_in_control.size())
				{
					if (ppschedulewrite[rb_idx][0].ue_selected == ue_in_control[ue_idx])
					{
						Real current_capacity_on_rb = log2(1 + CQI_read[ue_idx][rb_idx]);
						CQI_AVR[ue_idx][rb_idx] = (1. / _N_pf) * current_capacity_on_rb
											+ (1 - 1. / _N_pf) * CQI_AVR[ue_idx][rb_idx];
					}
					else
					{
						CQI_AVR[ue_idx][rb_idx] = (1 - 1. / _N_pf) * CQI_AVR[ue_idx][rb_idx];
					}
				}
			}
			else  // first loop
			{
				CQI_AVR[ue_idx][rb_idx] = log2(1 + CQI_read[ue_idx][rb_idx]);
			}
		}
	}

	// Calculate proportional fair scheduling metric
	if (Scheduling_Type == 0)
	{
		for (int ue_idx = 0; ue_idx < (int)ue_in_control.size(); ue_idx++)
		{
			for (int rb_idx = 0; rb_idx < num_rb; rb_idx++)
			{
				METRIC[ue_idx][rb_idx] = log2(1 + CQI_read[ue_idx][rb_idx]);
			}
		}
	}
	else if (Scheduling_Type == 1 || Scheduling_Type == 2 || Scheduling_Type == 3)
	{
		for (int ue_idx = 0; ue_idx < (int)ue_in_control.size(); ue_idx++)
		{
			Real alpha = 1;
			Real beta  = 1;

			for (int rb_idx = 0; rb_idx < num_rb; rb_idx++)
			{
				Real instantaneous_rate = log2(1 + CQI_read[ue_idx][rb_idx]);
				Real average_rate = CQI_AVR[ue_idx][rb_idx];

				// Avoid division by same value (happens on first scheduling)
				if (std::abs(instantaneous_rate - average_rate) < 1e-6) {
					// First scheduling: use instantaneous rate only
					METRIC[ue_idx][rb_idx] = pow(instantaneous_rate, alpha);
				} else {
					// Normal Proportional Fair metric
					METRIC[ue_idx][rb_idx] = pow(instantaneous_rate, alpha)
					                        / pow(average_rate, beta);
				}
			}
		}
	}
}


/*===================================================================
FUNCTION: Sector::Initialize_schedule_map( int )
===================================================================*/
void Sector::Initialize_schedule_map(void)
{
	for (int rb_idx = 0; rb_idx < num_rb; rb_idx++)
	{
		for (int stream = 0; stream < mx_ue_mumimo; stream++)
		{
			ppschedulewrite[rb_idx][stream].ue_selected = -1;
			ppschedulewrite[rb_idx][stream].mcs_selected = -1; // intitialize
			ppschedulewrite[rb_idx][stream].cqi_selected = -1; // intitialize
			ppschedulewrite[rb_idx][stream].capacity = -1; // intitialize
		}
	}


	if (Scheduling_Type == 0)
	{
		////////////////////////////////////////////////////////////////////
		if (W == NULL)
		{
			W = new MatrixXcReal[num_rb];

			for (int freq_idx = 0; freq_idx < num_rb; freq_idx++)
			{
				MatrixXcReal temp = MatrixXcReal::Zero(NUM_TX_Port, mx_ue_mumimo);
				W[freq_idx] = temp;
			}
		}
		////////////////////////////////////////////////////////////////////
	}
	else
	{
		////////////////////////////////////////////////////////////////////
		if (W == NULL)
		{
			W = new MatrixXcReal[num_rb];
		}
		////////////////////////////////////////////////////////////////////
	}

}

/*===================================================================
FUNCTION: Sector::Set_AVR_Cqi( void )
DESCRIPTION: Check the cqi for every allocated resource block, and Set it with minimum Cqi.
===================================================================*/
void Sector::Set_AVR_Cqi(void)   ////
{
	for (int ue_idx = 0; ue_idx < (int)ue_in_control.size(); ue_idx++)
	{
		int count = 0;
		Real sum = 0;

		int count_i = 0;
		Real sum_i = 0;

		int count_cap = 0;
		Real sum_cap = 0;
		Real e_sinr[15] = { 0. };
		int decision_CQI_index = 0;
		
		int determined_cqi;
		int determined_mcs;
		Real avr_cqi;
		Real avr_cqi_interfer = 0.;
		Real avr_capacity;

		MatrixXcReal identity_RX = MatrixXcReal::Identity(NUM_RX_Port, NUM_RX_Port);
		MatrixXcReal identity_TX = MatrixXcReal::Identity(NUM_TX_Port, NUM_TX_Port);

		MatrixXcReal CAP = MatrixXcReal::Zero(NUM_RX_Port, NUM_RX_Port);
		Real A;

		avr_capacity = 1.;


		Real linear_signal = pow(10., links[ue_in_control[ue_idx]].str_signal / 10.);
		MatrixXcReal B;
		MatrixXcReal eye = MatrixXcReal::Identity(NUM_RX_Port, NUM_RX_Port);
		Real linear_interference = pow(10., links[ue_in_control[ue_idx]].interference / 10.) + noise;
		B = eye;   // 2x2

		///////////////////////////////////////
		//// Average SINR
		//////////////////////////////////////
		for (int freq_idx = 0; freq_idx < num_rb; freq_idx++)
		{
			for (int stream = 0; stream < mx_ue_mumimo; stream++)
			{
				if (ppschedulewrite[freq_idx][stream].ue_selected == ue_in_control[ue_idx])
				{
					sum += CQI_read[ue_idx][freq_idx];
					count++;
				}
			}
		}

		avr_cqi = (sum / (Real)count);   //// <- linear average signal power
		avr_cqi_interfer = 1.;

		Real avr_SINR = 10. * log10(avr_cqi);

		// Apply OLLA offset if enabled
		if (g_olla_enable) {
			avr_SINR += ms[ue_in_control[ue_idx]].olla_offset;
		}
		// Realized-ESINR feedback correction (measurement-driven, see MS.h)
		if (g_matlab_bler && g_matlab_esinr_fb) {
			avr_SINR += ms[ue_in_control[ue_idx]].matlab_sinr_corr;
		}

		determined_mcs = determine_MCS(avr_SINR);// + 3;
		determined_cqi = determine_CQI(avr_SINR);// + 3;

		if (determined_mcs > 27)
		{
			determined_mcs = 27;
		}
		else if (determined_mcs < 0)
		{
			determined_mcs = 0;
		}

		for (int freq_idx = 0; freq_idx < num_rb; freq_idx++)
		{
			for (int stream = 0; stream < mx_ue_mumimo; stream++)
			{
				if (ppschedulewrite[freq_idx][stream].ue_selected == ue_in_control[ue_idx])
				{
					ppschedulewrite[freq_idx][stream].mcs_selected = determined_mcs;
					ppschedulewrite[freq_idx][stream].cqi_selected = determined_cqi;
					ppschedulewrite[freq_idx][stream].capacity = avr_capacity;

				}
			}
		}
	}
}


/*===================================================================
FUNCTION: Sector::Set_AVR_Cqi_Precoding_Based( void )
DESCRIPTION: Calculate SINR/CQI/MCS based on actual precoding matrix W
             OPTIMIZED: Exploits subband structure (same UEs and precoding per subband)
             Assumes inter-cell interference is known (perfect or with margin)
===================================================================*/
void Sector::Set_AVR_Cqi_Precoding_Based(void)
{
	// Subband configuration (must match Scheduling_algorithm_module_MU_MIMO)
	const int subband_size = 16;  // Number of RBs per subband
	int num_subband = (int)(num_rb / subband_size);
	if (num_rb % subband_size != 0) {
		num_subband += 1;
	}

	// Apply inter-cell interference margin (for future use)
	const Real interference_margin_linear = pow(10.0, INTERCELL_INTERFERENCE_MARGIN_DB / 10.0);

	// Pre-allocate Identity matrix (reused across all calculations)
	const MatrixXcReal Identity = MatrixXcReal::Identity(NUM_RX_Port, NUM_RX_Port);

	for (int ue_idx = 0; ue_idx < (int)ue_in_control.size(); ue_idx++)
	{
		const int ms_idx = ue_in_control[ue_idx];

		// Count number of scheduled subbands for this UE
		int count_subbands = 0;
		Real sum_sinr_linear = 0.0;
		Real min_sinr_linear = INFINITY;  // For minimum SINR mode

		// Calculate power levels ONCE per UE (reused across all subbands)
		const Real noise = dBm2linear(thermal_noise + (10. * log10(bandwidth)) + MS_noisefig);
		const Real linear_signal = dBm2linear(bs_maxpower - channel[(int)(self_idx/3)][ms_idx].pathloss_final);
		const Real linear_interference = dBm2linear(links[ms_idx].interference) * interference_margin_linear;

		// ============================================================================
		// PER-LAYER MCS (g_per_layer_mcs==1): give each of the UE's layer-streams its
		// own MCS from its own post-MMSE SINR, instead of the dominant-layer MCS for all.
		// Indexed by SCHEDULE_DECISION.layer_idx (0..R-1), matching the receiver side.
		// ============================================================================
		if (g_per_layer_mcs == 1)
		{
			Real sum_sinr_l[4] = {0,0,0,0};
			Real min_sinr_l[4] = {INFINITY,INFINITY,INFINITY,INFINITY};
			int  cnt_l[4]      = {0,0,0,0};

			for (int subband_idx = 0; subband_idx < num_subband; subband_idx++)
			{
				int rb_start = subband_idx * subband_size;
				int rb_end   = rb_start + subband_size;
				if (rb_end > num_rb) rb_end = num_rb;
				const int rb_idx = rb_start + (rb_end - rb_start) / 2;

				// Is this UE scheduled in this subband?
				bool present = false;
				for (int stream = 0; stream < mx_ue_mumimo; stream++)
					if (ppschedulewrite[rb_idx][stream].ue_selected == ms_idx) { present = true; break; }
				if (!present) continue;

				// H_bar and Ryy^-1 are shared across this UE's layers (depend only on all streams).
				const MatrixXcReal& H     = ms[ms_idx].H_m[0][rb_idx];
				const MatrixXcReal& W_rb  = W[rb_idx];
				const MatrixXcReal  H_bar = sqrt(linear_signal) * H * W_rb;
				const MatrixXcReal  Ryy_inv = (H_bar * H_bar.adjoint() + linear_interference * Identity).inverse();

				for (int stream = 0; stream < mx_ue_mumimo; stream++)
				{
					if (ppschedulewrite[rb_idx][stream].ue_selected != ms_idx) continue;
					int L = ppschedulewrite[rb_idx][stream].layer_idx;
					if (L < 0 || L >= 4) L = 0;

					const MatrixXcReal h_u = H_bar.col(stream);
					const MatrixXcReal u   = h_u.adjoint() * Ryy_inv;
					const Real signal_power = std::norm((u * h_u)(0, 0));
					Real interference_power = 0.0;
					for (int s2 = 0; s2 < mx_ue_mumimo; s2++)
					{
						if (s2 == stream || ppschedulewrite[rb_idx][s2].ue_selected == -1) continue;
						interference_power += std::norm((u * H_bar.col(s2))(0, 0));
					}
					interference_power += (linear_interference + noise) * (u * u.adjoint()).norm();
					const Real sinr = signal_power / interference_power;
					if (isnan(sinr) || isinf(sinr) || sinr <= 0) continue;

					sum_sinr_l[L] += sinr;
					if (sinr < min_sinr_l[L]) min_sinr_l[L] = sinr;
					cnt_l[L]++;
				}
			}

			for (int L = 0; L < 4; L++)
			{
				if (cnt_l[L] == 0) continue;
				Real sel = (g_use_min_sinr_for_mcs == 1) ? min_sinr_l[L] : sum_sinr_l[L] / (Real)cnt_l[L];
				Real dB  = 10.0 * log10(sel);
				if (g_olla_enable) dB += ms[ms_idx].olla_offset;
				// Realized-ESINR feedback correction (measurement-driven, see MS.h)
				if (g_matlab_bler && g_matlab_esinr_fb) dB += ms[ms_idx].matlab_sinr_corr;
				int mcs = determine_MCS(dB);
				if (mcs > 27) mcs = 27; else if (mcs < 0) mcs = 0;
				int cqi = determine_CQI(dB);

				for (int rb2 = 0; rb2 < num_rb; rb2++)
					for (int stream = 0; stream < mx_ue_mumimo; stream++)
						if (ppschedulewrite[rb2][stream].ue_selected == ms_idx &&
						    ppschedulewrite[rb2][stream].layer_idx == L)
						{
							ppschedulewrite[rb2][stream].mcs_selected = mcs;
							ppschedulewrite[rb2][stream].cqi_selected = cqi;
							ppschedulewrite[rb2][stream].capacity     = sel;
						}
			}
			continue;  // per-layer path done for this UE; skip the legacy single-MCS block
		}

		// Process each subband independently
		for (int subband_idx = 0; subband_idx < num_subband; subband_idx++)
		{
			// Determine RB range for this subband
			int rb_start = subband_idx * subband_size;
			int rb_end = rb_start + subband_size;
			if (rb_end > num_rb) {
				rb_end = num_rb;
			}

			// Representative RB for this subband (middle RB, same as used in scheduling)
			const int schedule_RB = rb_start + (rb_end - rb_start) / 2;

			// Check if this UE is scheduled in this subband (check representative RB)
			int stream_num = -1;
			for (int stream = 0; stream < mx_ue_mumimo; stream++)
			{
				if (ppschedulewrite[schedule_RB][stream].ue_selected == ms_idx)
				{

					stream_num = stream;
					break;
				}
			}

			// Skip if UE not scheduled in this subband
			if (stream_num == -1)
				continue;

			count_subbands++;

			// KEY OPTIMIZATION: Compute SINR ONCE per subband using representative RB
			// Same UEs and precoding are used across all RBs in the subband
			const int rb_idx = schedule_RB;

			// Get channel matrix H and precoding matrix W for representative RB
			const MatrixXcReal& H = ms[ms_idx].H_m[0][rb_idx];  // [NUM_RX_Port x NUM_TX_Port]
			const MatrixXcReal& W_rb = W[rb_idx];  // [NUM_TX_Port x mx_ue_mumimo]

			// Compute effective channel H_bar = sqrt(P) * H * W
			const MatrixXcReal H_bar = sqrt(linear_signal) * H * W_rb;  // [NUM_RX_Port x mx_ue_mumimo]

			// Extract desired user's channel vector
			const MatrixXcReal h_u = H_bar.col(stream_num);  // [NUM_RX_Port x 1]

			// Compute MMSE receiver filter
			// w = h_u^H * (H_bar * H_bar^H + σ²I)^(-1)
			const MatrixXcReal Ryy = H_bar * H_bar.adjoint() + linear_interference * Identity;
			const MatrixXcReal u = h_u.adjoint() * Ryy.inverse();  // [1 x NUM_RX_Port]

			// Desired signal power after receiver combining
			const ComplexReal A_scalar = (u * h_u)(0, 0);
			const Real signal_power = std::norm(A_scalar);  // |A|²

			// MU-MIMO interference power from other co-scheduled streams
			Real interference_power = 0.0;
			for (int stream = 0; stream < mx_ue_mumimo; stream++)
			{
				// Skip if it's the desired stream or no UE scheduled
				if (stream == stream_num || ppschedulewrite[rb_idx][stream].ue_selected == -1)
					continue;

				// Interference from stream i: B = u * h_i
				const MatrixXcReal h_i = H_bar.col(stream);
				const ComplexReal B_scalar = (u * h_i)(0, 0);
				interference_power += std::norm(B_scalar);  // |B|²
			}

			// Add external interference and noise (with MMSE enhancement)
			const Real u_norm_sq = (u * u.adjoint()).norm();
			interference_power += (linear_interference + noise) * u_norm_sq;

			// Compute SINR for this subband (representative RB)
			const Real subband_sinr_linear = signal_power / interference_power;

			// Debug: Check for NaN/Inf at subband level
			if (isnan(subband_sinr_linear) || isinf(subband_sinr_linear) || subband_sinr_linear <= 0) {
				cout << "=== NaN/Inf/Zero in subband_sinr_linear ===" << endl;
				cout << "  sector_idx: " << self_idx << ", ms_idx: " << ms_idx << ", t: " << t << endl;
				cout << "  subband_idx: " << subband_idx << ", rb_idx: " << rb_idx << ", stream_num: " << stream_num << endl;
				cout << "  signal_power: " << signal_power << endl;
				cout << "  interference_power: " << interference_power << endl;
				cout << "  u_norm_sq: " << u_norm_sq << endl;
				cout << "  H norm: " << H.norm() << ", W_rb norm: " << W_rb.norm() << endl;
				cout << "  H_bar norm: " << H_bar.norm() << ", h_u norm: " << h_u.norm() << endl;
				cout << "============================================" << endl;
			}

			sum_sinr_linear += subband_sinr_linear;
			if (subband_sinr_linear < min_sinr_linear) {
				min_sinr_linear = subband_sinr_linear;
			}
		}

		// If UE has scheduled subbands, calculate average or minimum SINR and determine MCS/CQI
		if (count_subbands > 0)
		{
			// Select SINR based on aggregation mode:
			// g_use_min_sinr_for_mcs = 0: Average SINR (default, optimistic)
			// g_use_min_sinr_for_mcs = 1: Minimum SINR (conservative, reduces BLER)
			Real selected_sinr_linear;
			if (g_use_min_sinr_for_mcs == 1) {
				selected_sinr_linear = min_sinr_linear;
			} else {
				selected_sinr_linear = sum_sinr_linear / (Real)count_subbands;
			}
			Real avr_sinr_dB = 10.0 * log10(selected_sinr_linear);

			// Debug: Check for NaN in avr_sinr_dB
			if (isnan(avr_sinr_dB) || isinf(avr_sinr_dB)) {
				cout << "=== NaN/Inf detected in avr_sinr_dB ===" << endl;
				cout << "  sector_idx: " << self_idx << ", ms_idx: " << ms_idx << ", t: " << t << endl;
				cout << "  selected_sinr_linear: " << selected_sinr_linear << endl;
				cout << "  sum_sinr_linear: " << sum_sinr_linear << endl;
				cout << "  min_sinr_linear: " << min_sinr_linear << endl;
				cout << "  count_subbands: " << count_subbands << endl;
				cout << "  g_use_min_sinr_for_mcs: " << g_use_min_sinr_for_mcs << endl;
				cout << "  linear_signal: " << linear_signal << endl;
				cout << "  linear_interference: " << linear_interference << endl;
				cout << "  noise: " << noise << endl;
				cout << "  links[ms_idx].str_signal: " << links[ms_idx].str_signal << endl;
				cout << "  links[ms_idx].interference: " << links[ms_idx].interference << endl;
				cout << "========================================" << endl;
			}
	
			// Apply OLLA (Outer Loop Link Adaptation) offset if enabled
			if (g_olla_enable) {
				avr_sinr_dB += ms[ms_idx].olla_offset;
			}
			// Realized-ESINR feedback correction (measurement-driven, see MS.h)
			if (g_matlab_bler && g_matlab_esinr_fb) {
				avr_sinr_dB += ms[ms_idx].matlab_sinr_corr;
			}

			// Determine MCS and CQI based on adjusted SINR
			int determined_mcs = determine_MCS(avr_sinr_dB);
			int determined_cqi = determine_CQI(avr_sinr_dB);

			// Clip MCS to valid range [0, 27]
			if (determined_mcs > 27)
				determined_mcs = 27;
			else if (determined_mcs < 0)
				determined_mcs = 0;

			// Apply MCS/CQI decision to all scheduled RBs for this UE
			// This is fast because we're just updating existing schedule
			for (int rb_idx = 0; rb_idx < num_rb; rb_idx++)
			{
				for (int stream = 0; stream < mx_ue_mumimo; stream++)
				{
					if (ppschedulewrite[rb_idx][stream].ue_selected == ms_idx)
					{
						ppschedulewrite[rb_idx][stream].mcs_selected = determined_mcs;
						ppschedulewrite[rb_idx][stream].cqi_selected = determined_cqi;
						ppschedulewrite[rb_idx][stream].capacity = selected_sinr_linear;
					}
				}
			}
		}
	}
}


/*===================================================================
FUNCTION: Sector::Initial_Setting( int )

DESCRIPTION:
===================================================================*/
void Sector::Initial_Setting(int rb_idx)
{
	scheduled_ue.clear();
	selected_ue.clear();
	prev_sum_SINR = -pow(10., 10.);
	sum_of_SINR = 0.;

	for (int stream = 0; stream < mx_ue_mumimo; stream++)
	{
		ppschedulewrite[rb_idx][stream].ue_selected = -1;
	}
}


#if 0  // unused: first_max_ue, SUS_procedure, Get_User_Set, Filter_User_Set — replaced by _threadsafe versions
/*===================================================================
FUNCTION: Sector::first_max_ue( int rb_idx )

DESCRIPTION:
  Selects the first UE with maximum metric for SUS algorithm.
  This UE serves as the anchor for semi-orthogonal user selection.

INPUT:  METRIC[ue_idx][rb_idx] - Proportional fairness metric
OUTPUT: selected_ue - First selected UE with highest metric
===================================================================*/
void Sector::first_max_ue(int rb_idx)
{
	selected_ue.clear();

	// Safety check: Ensure METRIC is allocated
	if (METRIC == NULL) {
		printf("[ERROR] first_max_ue: METRIC is NULL at sector %d, rb_idx %d\n", self_idx, rb_idx);
		printf("[ERROR] Read_Ch_Feedback() or Read_CSI_Feedback() must be called first!\n");
		selected_ue.push_back(-1);
		return;
	}

	// Safety check: Ensure rb_idx is valid
	if (rb_idx < 0 || rb_idx >= num_rb) {
		printf("[ERROR] first_max_ue: Invalid rb_idx %d (num_rb=%d) at sector %d\n", rb_idx, num_rb, self_idx);
		selected_ue.push_back(-1);
		return;
	}

	Real MAX_SINR = -pow(10., 10.);
	int    MAX_UE = -1;

	// Find UE with maximum metric (less than previous maximum for fairness)
	for (int ue_idx = 0; ue_idx < (int)ue_in_control.size(); ue_idx++)
	{
		// Bounds check to prevent segmentation fault
		if (ue_idx >= 80) {  // fixed_mem = 80
			printf("[ERROR] first_max_ue: ue_idx %d exceeds fixed_mem (80) at sector %d\n", ue_idx, self_idx);
			break;
		}

		if (METRIC[ue_idx][rb_idx] > MAX_SINR &&  METRIC[ue_idx][rb_idx] < prev_max_sinr)
		{
			MAX_SINR = METRIC[ue_idx][rb_idx];
			MAX_UE = ue_idx;
		}
	}

	selected_ue.push_back(MAX_UE);
	prev_max_sinr = MAX_SINR;
}

/*===================================================================
FUNCTION: Sector::SUS_procedure( int rb_idx )

DESCRIPTION:
  Semi-Orthogonal User Selection (SUS) algorithm for MU-MIMO.
  Iteratively selects UEs that are semi-orthogonal to already selected UEs.

  Algorithm:
  1. First UE selected by first_max_ue() (highest metric)
  2. Find candidate UEs semi-orthogonal to selected UEs (correlation < epsilon)
  3. Among candidates, select UE with highest metric
  4. Repeat until mx_ue_mumimo UEs selected or no more candidates

INPUT:  selected_ue - Initial UE from first_max_ue()
OUTPUT: scheduled_ue - Final set of semi-orthogonal UEs for MU-MIMO
===================================================================*/
void Sector::SUS_procedure(int rb_idx)
{
	// SUS = Semi-Orthogonal User Selection (Optimized incremental filtering)

	// If first UE selection failed, exit
	if (selected_ue[0] == -1)
	{
		return;
	}

	// Initialize user_set with all UEs semi-orthogonal to the first selected UE
	Get_User_Set(rb_idx);

	// Iteratively select additional UEs
	while (selected_ue.size() < mx_ue_mumimo && selected_ue.size() < ue_in_control.size())
	{
		// If user_set is empty, no more semi-orthogonal UEs available
		if (user_set.empty())
		{
			break;
		}

		// Among remaining semi-orthogonal UEs, select one with maximum metric
		Real MAX_SINR = -500;
		int    MAX_UE = -1;
		for (int j = 0; j < (int)user_set.size(); j++)
		{
			if (METRIC[user_set[j]][rb_idx] > MAX_SINR)
			{
				MAX_SINR = METRIC[user_set[j]][rb_idx];
				MAX_UE = user_set[j];
			}
		}

		// If no valid UE found, stop
		if (MAX_UE == -1)
		{
			break;
		}

		// Add selected UE to the list
		selected_ue.push_back(MAX_UE);

		// Remove the selected UE from user_set and filter out UEs that are not
		// semi-orthogonal to the newly selected UE (incremental filtering)
		Filter_User_Set(rb_idx, MAX_UE);
	}

	// Copy to scheduled_ue
	scheduled_ue = selected_ue;
}


/*===================================================================
FUNCTION: Sector::Get_User_Set( int rb_idx )

DESCRIPTION:
  Finds UEs that are semi-orthogonal to all currently selected UEs.

  Semi-orthogonality criterion:
  - For all selected UEs j, correlation |<v_i, v_j>| < epsilon
  - Where v_i, v_j are normalized channel/PMI vectors

  TDD mode: Uses actual channel matrices (CSI_matrix_read)
  FDD mode: Uses quantized PMI vectors (PMI_vector_read)

INPUT:  selected_ue - Currently selected UEs
        epsilon - Orthogonality threshold (global)
OUTPUT: user_set - Set of semi-orthogonal candidate UEs
===================================================================*/
void Sector::Get_User_Set(int rb_idx)
{
	user_set.clear();

	// Safety check: Ensure required data is allocated
	if (METRIC == NULL) {
		printf("[ERROR] Get_User_Set: METRIC is NULL at sector %d, rb_idx %d\n", self_idx, rb_idx);
		return;
	}

	if (TDD_mode == 1 && CSI_matrix_read == NULL) {
		printf("[ERROR] Get_User_Set: CSI_matrix_read is NULL (TDD mode) at sector %d\n", self_idx);
		return;
	}

	if (TDD_mode == 0 && PMI_vector_read == NULL) {
		printf("[ERROR] Get_User_Set: PMI_vector_read is NULL (FDD mode) at sector %d\n", self_idx);
		return;
	}

	int  MIN_SUS_user_idx =  0;
	Real MIN_inner_prod   = -1;

	// Check each UE for semi-orthogonality
	for (int ue_idx = 0; ue_idx < (int)ue_in_control.size(); ue_idx++)
	{
		// Bounds check
		if (ue_idx >= 80) {
			printf("[ERROR] Get_User_Set: ue_idx %d exceeds fixed_mem (80) at sector %d\n", ue_idx, self_idx);
			break;
		}
		bool already_selected = false;

		// Skip if UE is already selected
		for (int _ue_idx = 0; _ue_idx < (int)selected_ue.size(); _ue_idx++)
		{
			if (ue_idx == selected_ue[_ue_idx])
			{
				already_selected = true;
				break;  // Early exit
			}
		}

		bool semi_orthogonal = false;

		if (!already_selected)
		{
			semi_orthogonal = true;  // Assume semi-orthogonal until proven otherwise


			// Check orthogonality with all selected UEs
			for (int j = 0; j < (int)selected_ue.size(); j++)
			{
				// Bounds check for selected_ue[j]
				if (selected_ue[j] < 0 || selected_ue[j] >= 80) {
					printf("[WARNING] Get_User_Set: selected_ue[%d]=%d out of bounds at sector %d\n",
					       j, selected_ue[j], self_idx);
					semi_orthogonal = false;
					break;
				}

				VectorXcReal A, B;

				// TDD mode: Use channel-based vectors, FDD mode: Use PMI vectors
				if (TDD_mode == 1) {
					// Extract effective channel from CSI matrix
					MatrixXcReal H_selected  = CSI_matrix_read[selected_ue[j]][rb_idx];
					MatrixXcReal H_candidate = CSI_matrix_read[        ue_idx][rb_idx];

					// Use dominant eigenvector (right singular vector) for each channel
					// This captures the strongest spatial mode across all receive antennas
					MatrixXcReal HH_selected  = H_selected.adjoint()  * H_selected;   // T×T
					MatrixXcReal HH_candidate = H_candidate.adjoint() * H_candidate;  // T×T

					// Compute eigenvalue decomposition to get dominant eigenvector
					SelfAdjointEigenSolver<MatrixXcReal> es_selected(HH_selected);
					SelfAdjointEigenSolver<MatrixXcReal> es_candidate(HH_candidate);

					// Get the eigenvector corresponding to the largest eigenvalue (last column)
					A = es_selected.eigenvectors().col(es_selected.eigenvectors().cols() - 1);
					B = es_candidate.eigenvectors().col(es_candidate.eigenvectors().cols() - 1);
				} else {
					// FDD mode: Use PMI vectors (already optimal precoding vectors)
					// For multi-layer PMI, use dominant layer (col 0) for pairing correlation.
					A = PMI_vector_read[selected_ue[j]][rb_idx].col(0);
					B = PMI_vector_read[ue_idx        ][rb_idx].col(0);
				}

				// Normalize vectors
				Real norm_A = A.norm();
				Real norm_B = B.norm();

				if (norm_A < 1e-12 || norm_B < 1e-12) {
					// Skip if either vector is nearly zero
					semi_orthogonal = false;
					break;
				}

				A = A / norm_A;
				B = B / norm_B;

				// Compute correlation: |<A, B>|
				ComplexReal inner_prod_complex = (A.adjoint() * B)(0);
				Real inner_prod = std::abs(inner_prod_complex);

				// Check if correlation exceeds threshold (not semi-orthogonal)
				if (inner_prod > epsilon)
				{
					semi_orthogonal = false;
					// Continue to track MIN_inner_prod even if not semi-orthogonal
				}

				// Track minimum inner product for debugging/fallback
				if (MIN_inner_prod < 0 || MIN_inner_prod > inner_prod)
				{
					MIN_inner_prod = inner_prod;
					MIN_SUS_user_idx = ue_idx;
				}
			}

		}

		if (semi_orthogonal)
		{
			user_set.push_back(ue_idx);
		}
	}
}


/*===================================================================
FUNCTION: Sector::Filter_User_Set( int rb_idx, int newly_selected_ue )

DESCRIPTION:
  Incrementally filters user_set to remove:
  1. The newly selected UE itself
  2. UEs that are NOT semi-orthogonal to the newly selected UE

  This is more efficient than rebuilding the entire user_set from scratch,
  as it only checks against the one newly added UE.

INPUT:  rb_idx - Resource block index
        newly_selected_ue - The UE just added to selected_ue
        user_set - Current set of candidate UEs (modified in-place)
OUTPUT: user_set - Filtered to remain semi-orthogonal to all selected UEs
===================================================================*/
void Sector::Filter_User_Set(int rb_idx, int newly_selected_ue)
{
	// Safety check
	if (newly_selected_ue < 0 || newly_selected_ue >= 80) {
		printf("[WARNING] Filter_User_Set: newly_selected_ue=%d out of bounds\n", newly_selected_ue);
		user_set.clear();
		return;
	}

	// Extract the channel/PMI vector for the newly selected UE
	VectorXcReal A_new;

	if (TDD_mode == 1) {
		MatrixXcReal H_new = CSI_matrix_read[newly_selected_ue][rb_idx];
		MatrixXcReal HH_new = H_new.adjoint() * H_new;
		SelfAdjointEigenSolver<MatrixXcReal> es_new(HH_new);
		A_new = es_new.eigenvectors().col(es_new.eigenvectors().cols() - 1);
	} else {
		A_new = PMI_vector_read[newly_selected_ue][rb_idx].col(0);  // dominant layer
	}

	// Normalize
	Real norm_A = A_new.norm();
	if (norm_A < 1e-12) {
		// If newly selected UE has invalid channel, clear user_set
		user_set.clear();
		return;
	}
	A_new = A_new / norm_A;

	// Filter user_set: keep only UEs that are semi-orthogonal to the newly selected UE
	std::vector<int> filtered_set;
	for (int i = 0; i < (int)user_set.size(); i++)
	{
		int candidate_ue = user_set[i];

		// Remove the newly selected UE itself
		if (candidate_ue == newly_selected_ue) {
			continue;
		}

		// Bounds check
		if (candidate_ue < 0 || candidate_ue >= 80) {
			continue;
		}

		// Check semi-orthogonality with newly selected UE
		VectorXcReal B;

		if (TDD_mode == 1) {
			MatrixXcReal H_candidate = CSI_matrix_read[candidate_ue][rb_idx];
			MatrixXcReal HH_candidate = H_candidate.adjoint() * H_candidate;
			SelfAdjointEigenSolver<MatrixXcReal> es_candidate(HH_candidate);
			B = es_candidate.eigenvectors().col(es_candidate.eigenvectors().cols() - 1);
		} else {
			B = PMI_vector_read[candidate_ue][rb_idx].col(0);  // dominant layer
		}

		// Normalize
		Real norm_B = B.norm();
		if (norm_B < 1e-12) {
			continue;  // Skip invalid channels
		}
		B = B / norm_B;

		// Compute correlation
		ComplexReal inner_prod_complex = (A_new.adjoint() * B)(0);
		Real inner_prod = std::abs(inner_prod_complex);

		// Keep only if semi-orthogonal
		if (inner_prod <= epsilon) {
			filtered_set.push_back(candidate_ue);
		}
	}

	// Replace user_set with filtered version
	user_set = filtered_set;
}
#endif  // unused: first_max_ue, SUS_procedure, Get_User_Set, Filter_User_Set


/*===================================================================
Thread-safe versions of scheduling functions for parallel subband processing
These functions use explicit local variables instead of member variables
===================================================================*/

void Sector::first_max_ue_threadsafe(int rb_idx, vector<int>& local_selected_ue)
{
	local_selected_ue.clear();

	if (METRIC == NULL || rb_idx < 0 || rb_idx >= num_rb) {
		local_selected_ue.push_back(-1);
		return;
	}

	Real MAX_SINR = -pow(10., 10.);
	int MAX_UE = -1;

	// Simply find the UE with maximum METRIC (no prev_max_sinr needed)
	for (int ue_idx = 0; ue_idx < (int)ue_in_control.size(); ue_idx++)
	{
		if (ue_idx >= 80) break;

		if (METRIC[ue_idx][rb_idx] > MAX_SINR)
		{
			MAX_SINR = METRIC[ue_idx][rb_idx];
			MAX_UE = ue_idx;
		}
	}

	local_selected_ue.push_back(MAX_UE);
}


void Sector::Get_User_Set_threadsafe(int rb_idx, vector<int>& local_selected_ue, vector<int>& local_user_set)
{
	local_user_set.clear();

	if (METRIC == NULL) return;
	if (TDD_mode == 1 && CSI_matrix_read == NULL) return;
	if (TDD_mode == 0 && PMI_vector_read == NULL) return;

	for (int ue_idx = 0; ue_idx < (int)ue_in_control.size(); ue_idx++)
	{
		if (ue_idx >= 80) break;

		bool already_selected = false;
		for (int _ue_idx = 0; _ue_idx < (int)local_selected_ue.size(); _ue_idx++)
		{
			if (ue_idx == local_selected_ue[_ue_idx])
			{
				already_selected = true;
				break;
			}
		}

		if (!already_selected)
		{
			bool semi_orthogonal = true;

			for (int j = 0; j < (int)local_selected_ue.size(); j++)
			{
				if (local_selected_ue[j] < 0 || local_selected_ue[j] >= 80) {
					semi_orthogonal = false;
					break;
				}

				// Use pre-computed vectors (both TDD and FDD)
				// TDD: PMI_vector_read contains dominant eigenvectors (from Store_CSI_for_TDD)
				// FDD: PMI_vector_read contains quantized PMI vectors (from Quantization_of_Ch)
				// For multi-layer PMI, pairing is based on dominant layer (col 0).
				VectorXcReal A = PMI_vector_read[local_selected_ue[j]][rb_idx].col(0);
				VectorXcReal B = PMI_vector_read[ue_idx][rb_idx].col(0);

				Real norm_A = A.norm();
				Real norm_B = B.norm();

				if (norm_A < 1e-12 || norm_B < 1e-12) {
					semi_orthogonal = false;
					break;
				}

				// Normalize (in case not already normalized)
				A = A / norm_A;
				B = B / norm_B;

				ComplexReal inner_prod_complex = (A.adjoint() * B)(0);
				Real inner_prod = std::abs(inner_prod_complex);

				if (inner_prod > epsilon)
				{
					semi_orthogonal = false;
				}
			}

			if (semi_orthogonal)
			{
				local_user_set.push_back(ue_idx);
			}
		}
	}
}


void Sector::Filter_User_Set_threadsafe(int rb_idx, int newly_selected_ue, vector<int>& local_user_set)
{
	if (newly_selected_ue < 0 || newly_selected_ue >= 80) {
		local_user_set.clear();
		return;
	}

	// Use pre-computed vectors (both TDD and FDD)
	// TDD: PMI_vector_read contains dominant eigenvectors (from Store_CSI_for_TDD)
	// FDD: PMI_vector_read contains quantized PMI vectors (from Quantization_of_Ch)
	// For multi-layer PMI, use dominant layer (col 0) for pairing correlation.
	VectorXcReal A_new = PMI_vector_read[newly_selected_ue][rb_idx].col(0);

	Real norm_A = A_new.norm();
	if (norm_A < 1e-12) {
		local_user_set.clear();
		return;
	}
	A_new = A_new / norm_A;

	std::vector<int> filtered_set;
	for (int i = 0; i < (int)local_user_set.size(); i++)
	{
		int candidate_ue = local_user_set[i];

		if (candidate_ue == newly_selected_ue) continue;
		if (candidate_ue < 0 || candidate_ue >= 80) continue;

		// Use pre-computed vectors (both TDD and FDD) — dominant layer for pairing
		VectorXcReal B = PMI_vector_read[candidate_ue][rb_idx].col(0);

		Real norm_B = B.norm();
		if (norm_B < 1e-12) continue;
		B = B / norm_B;

		ComplexReal inner_prod_complex = (A_new.adjoint() * B)(0);
		Real inner_prod = std::abs(inner_prod_complex);

		if (inner_prod <= epsilon) {
			filtered_set.push_back(candidate_ue);
		}
	}

	local_user_set = filtered_set;
}


void Sector::SUS_procedure_threadsafe(int rb_idx, vector<int>& local_selected_ue, vector<int>& local_scheduled_ue, vector<int>& local_user_set)
{
	if (local_selected_ue[0] == -1)
	{
		return;
	}

	Get_User_Set_threadsafe(rb_idx, local_selected_ue, local_user_set);

	while (local_selected_ue.size() < mx_ue_mumimo && local_selected_ue.size() < ue_in_control.size())
	{
		if (local_user_set.empty())
		{
			break;
		}

		Real MAX_SINR = -500;
		int MAX_UE = -1;
		for (int j = 0; j < (int)local_user_set.size(); j++)
		{
			if (METRIC[local_user_set[j]][rb_idx] > MAX_SINR)
			{
				MAX_SINR = METRIC[local_user_set[j]][rb_idx];
				MAX_UE = local_user_set[j];
			}
		}

		if (MAX_UE == -1)
		{
			break;
		}

		local_selected_ue.push_back(MAX_UE);
		Filter_User_Set_threadsafe(rb_idx, MAX_UE, local_user_set);
	}

	local_scheduled_ue = local_selected_ue;
}


/*===================================================================
OPTIMIZED versions with pre-computed eigenvector cache
These avoid redundant eigenvalue decompositions
===================================================================*/

void Sector::Get_User_Set_threadsafe_cached(int rb_idx, vector<int>& local_selected_ue, vector<int>& local_user_set,
                                              const vector<VectorXcReal>& eigenvector_cache)
{
	local_user_set.clear();

	if (METRIC == NULL) return;

	for (int ue_idx = 0; ue_idx < (int)ue_in_control.size(); ue_idx++)
	{
		if (ue_idx >= 80) break;

		bool already_selected = false;
		for (int _ue_idx = 0; _ue_idx < (int)local_selected_ue.size(); _ue_idx++)
		{
			if (ue_idx == local_selected_ue[_ue_idx])
			{
				already_selected = true;
				break;
			}
		}

		if (!already_selected)
		{
			bool semi_orthogonal = true;

			for (int j = 0; j < (int)local_selected_ue.size(); j++)
			{
				if (local_selected_ue[j] < 0 || local_selected_ue[j] >= 80) {
					semi_orthogonal = false;
					break;
				}

				// Use pre-computed eigenvectors (already normalized)
				const VectorXcReal& A = eigenvector_cache[local_selected_ue[j]];
				const VectorXcReal& B = eigenvector_cache[ue_idx];

				// Check validity
				if (A.norm() < 1e-12 || B.norm() < 1e-12) {
					semi_orthogonal = false;
					break;
				}

				ComplexReal inner_prod_complex = (A.adjoint() * B)(0);
				Real inner_prod = std::abs(inner_prod_complex);

				if (inner_prod > epsilon)
				{
					semi_orthogonal = false;
				}
			}

			if (semi_orthogonal)
			{
				local_user_set.push_back(ue_idx);
			}
		}
	}
}


void Sector::Filter_User_Set_threadsafe_cached(int rb_idx, int newly_selected_ue, vector<int>& local_user_set,
                                                 const vector<VectorXcReal>& eigenvector_cache)
{
	if (newly_selected_ue < 0 || newly_selected_ue >= 80) {
		local_user_set.clear();
		return;
	}

	// Use pre-computed eigenvector (already normalized)
	const VectorXcReal& A_new = eigenvector_cache[newly_selected_ue];

	if (A_new.norm() < 1e-12) {
		local_user_set.clear();
		return;
	}

	std::vector<int> filtered_set;
	for (int i = 0; i < (int)local_user_set.size(); i++)
	{
		int candidate_ue = local_user_set[i];

		if (candidate_ue == newly_selected_ue) continue;
		if (candidate_ue < 0 || candidate_ue >= 80) continue;

		const VectorXcReal& B = eigenvector_cache[candidate_ue];

		if (B.norm() < 1e-12) continue;

		ComplexReal inner_prod_complex = (A_new.adjoint() * B)(0);
		Real inner_prod = std::abs(inner_prod_complex);

		if (inner_prod <= epsilon) {
			filtered_set.push_back(candidate_ue);
		}
	}

	local_user_set = filtered_set;
}


void Sector::SUS_procedure_threadsafe_cached(int rb_idx, vector<int>& local_selected_ue, vector<int>& local_scheduled_ue,
                                               vector<int>& local_user_set, const vector<VectorXcReal>& eigenvector_cache)
{
	if (local_selected_ue[0] == -1)
	{
		return;
	}

	Get_User_Set_threadsafe_cached(rb_idx, local_selected_ue, local_user_set, eigenvector_cache);

	while (local_selected_ue.size() < mx_ue_mumimo && local_selected_ue.size() < ue_in_control.size())
	{
		if (local_user_set.empty())
		{
			break;
		}

		Real MAX_SINR = -500;
		int MAX_UE = -1;
		for (int j = 0; j < (int)local_user_set.size(); j++)
		{
			if (METRIC[local_user_set[j]][rb_idx] > MAX_SINR)
			{
				MAX_SINR = METRIC[local_user_set[j]][rb_idx];
				MAX_UE = local_user_set[j];
			}
		}

		if (MAX_UE == -1)
		{
			break;
		}

		local_selected_ue.push_back(MAX_UE);
		Filter_User_Set_threadsafe_cached(rb_idx, MAX_UE, local_user_set, eigenvector_cache);
	}

	local_scheduled_ue = local_selected_ue;
}


/*===================================================================
FUNCTION: Sector::Transmit_Precoding( int )

DESCRIPTION:
===================================================================*/
void Sector::Transmit_Precoding(int rb_idx, const vector<int>& local_scheduled_ue)
{
	MatrixXcReal Identity = Matrix2cReal::Identity();

	// Per-UE rank. For eType II with rank-adaptive, each UE contributes self_RI
	// stream-layers; otherwise the global g_type2_rank. Cap at mx_ue_mumimo.
	auto per_ue_rank = [&](int ms_idx) -> int {
		if (g_codebook_type != 3) return 1;
		int R = (g_rank_adaptive == 1) ? ms[ms_idx].self_RI : g_type2_rank;
		if (R < 1) R = 1;
		return R;
	};

	// Count total valid streams capped at mx_ue_mumimo
	int num_valid_streams = 0;
	for (int i = 0; i < (int)local_scheduled_ue.size(); i++) {
		if (local_scheduled_ue[i] > -1) {
			int ms_idx = ue_in_control[local_scheduled_ue[i]];
			num_valid_streams += per_ue_rank(ms_idx);
			if (num_valid_streams >= mx_ue_mumimo) {
				num_valid_streams = mx_ue_mumimo;
				break;
			}
		}
	}

	if (num_valid_streams == 0) {
		W[rb_idx] = MatrixXcReal::Zero(NUM_TX_Port, mx_ue_mumimo);
		return;
	}

	// Create Lambda matrix — one row per (UE, layer) stream
	MatrixXcReal Lambda = MatrixXcReal::Zero(num_valid_streams, NUM_TX_Port);

	// Fill Lambda: for each scheduled UE, stack R_ue rows (col 0..R_ue-1 of PMI_W)
	int stream_row = 0;
	for (int i = 0; i < (int)local_scheduled_ue.size() && stream_row < num_valid_streams; i++)
	{
		if (local_scheduled_ue[i] < 0) continue;
		int ms_idx = ue_in_control[local_scheduled_ue[i]];
		int R_ue  = per_ue_rank(ms_idx);
		const MatrixXcReal& PMIm = PMI_vector_read[local_scheduled_ue[i]][rb_idx];
		int R_avail = (int)PMIm.cols();
		int R_use   = (R_avail < R_ue) ? R_avail : R_ue;
		for (int r = 0; r < R_use && stream_row < num_valid_streams; r++)
		{
			Lambda.row(stream_row) = PMIm.col(r).transpose();
			stream_row++;
		}
	}
	num_valid_streams = stream_row;
	if (num_valid_streams == 0) {
		W[rb_idx] = MatrixXcReal::Zero(NUM_TX_Port, mx_ue_mumimo);
		return;
	}
	Lambda.conservativeResize(num_valid_streams, Eigen::NoChange);
	int num_valid_ue = num_valid_streams;

	// Regularized Zero-Forcing precoder for the selected UEs
	if ( num_valid_ue > 0 )
	{
		// Compute Lambda * Lambda^H
		MatrixXcReal LambdaLambdaH = Lambda * Lambda.adjoint();
		Real det = LambdaLambdaH.determinant().real();
		const Real EPSILON = 1e-10;

		// Add Tikhonov regularization (Regularized Zero-Forcing, MU-MIMO)
		Real alpha_reg = (std::abs(det) < EPSILON) ? 1e-4 : 1e-6;
		MatrixXcReal Identity_reg = MatrixXcReal::Identity(LambdaLambdaH.rows(), LambdaLambdaH.cols());
		LambdaLambdaH += alpha_reg * Identity_reg;

		// RZF precoder: (NUM_TX_Port × num_valid_ue)
		MatrixXcReal W_valid = Lambda.adjoint() * LambdaLambdaH.inverse();

		// Per-column power normalization
		VectorXReal sqrt_p = VectorXReal::Zero(num_valid_ue);
		for (int i = 0; i < num_valid_ue; i++)
		{
			Real col_norm = W_valid.col(i).norm();
			sqrt_p(i) = (col_norm < EPSILON || std::isnan(col_norm) || std::isinf(col_norm))
			            ? 1.0 : (1.0 / col_norm);
		}
		W_valid = W_valid * sqrt_p.asDiagonal();

		// Pad to full mx_ue_mumimo columns (consumer code expects uniform width).
		// With multi-layer UEs the stream layout matches ppschedulewrite:
		//   stream 0..R-1 = UE0's R layers, stream R..2R-1 = UE1's R layers, ...
		// W_valid columns are in that exact order already (as they came from Lambda rows).
		MatrixXcReal W_full = MatrixXcReal::Zero(NUM_TX_Port, mx_ue_mumimo);
		int ncopy = (num_valid_ue < mx_ue_mumimo) ? num_valid_ue : mx_ue_mumimo;
		W_full.leftCols(ncopy) = W_valid.leftCols(ncopy);
		W[rb_idx] = W_full;
	}
}


/*===================================================================
FUNCTION: Sector::Transmit_Precoding_TDD( int rb_idx )

DESCRIPTION:
  TDD reciprocity-based precoding using full channel matrices.
  Computes Zero-Forcing (ZF) or Regularized Zero-Forcing (RZF) precoding
  directly from channel knowledge without codebook quantization.

  For MU-MIMO with K users and M TX antennas:
  - H_combined: K × M matrix (each row is one UE's channel)
  - W_ZF = H^H (H H^H)^-1  (Zero-Forcing)
  - W_RZF = H^H (H H^H + αI)^-1  (Regularized ZF / MMSE)

INPUT:  CSI_matrix_read[ue_idx][rb_idx] - Channel matrices from scheduled UEs
OUTPUT: W[rb_idx] - Precoding matrix (NUM_TX_Port × num_valid_ue)
===================================================================*/
void Sector::Transmit_Precoding_TDD(int rb_idx, const vector<int>& local_scheduled_ue)
{
	MatrixXcReal Identity = Matrix2cReal::Identity();

	// Count valid scheduled UEs (those with index > -1)
	int num_valid_ue = 0;
	vector<int> valid_ue_indices;  // Keep track of which UEs are valid

	for (int i = 0; i < (int)local_scheduled_ue.size(); i++) {
		if (local_scheduled_ue[i] > -1) {
			num_valid_ue++;
			valid_ue_indices.push_back(i);
		}
	}

	// Early exit if no valid UEs scheduled
	if (num_valid_ue == 0) {
		return;
	}

	// Build combined channel matrix H_combined (num_valid_ue × NUM_TX_Port)
	// Each row is one UE's effective channel (using dominant eigenmode or first RX antenna)
	MatrixXcReal H_combined = MatrixXcReal::Zero(num_valid_ue, NUM_TX_Port);

	for (int i = 0; i < num_valid_ue; i++)
	{
		int ue_control_idx = local_scheduled_ue[valid_ue_indices[i]];

		// Get full channel matrix from CSI feedback (NUM_RX_Port × NUM_TX_Port)
		MatrixXcReal H_full = CSI_matrix_read[ue_control_idx][rb_idx];

		// Option 1: Use first receive antenna's channel (SISO equivalent)
		// This is simplest and works well for single-stream per user
		H_combined.row(i) = H_full.row(0);

		// Option 2 (Alternative): Use dominant eigenvector via SVD
		// Uncomment below to use SVD-based approach
		/*
		JacobiSVD<MatrixXcReal> svd(H_full, ComputeThinU | ComputeThinV);
		RowVectorXcReal dominant_channel = svd.matrixV().col(0).transpose();
		H_combined.row(i) = dominant_channel;
		*/
	}

	// Compute H * H^H for ZF precoding
	MatrixXcReal H_HH = H_combined * H_combined.adjoint();

	// Check for near-singular matrix
	Real det = H_HH.determinant().real();
	const Real EPSILON = 1e-10;

	if (std::abs(det) < EPSILON) {
		printf("[WARNING] TDD: Nearly singular H*H^H (det=%e) at rb_idx=%d, sector=%d\n",
		       det, rb_idx, self_idx);

		// Print channel correlation information
		printf("[DEBUG] TDD: num_valid_ue=%d, H_combined size: %lldx%lld\n",
		       num_valid_ue, (long long)H_combined.rows(), (long long)H_combined.cols());

		for (int i = 0; i < num_valid_ue; i++) {
			printf("  UE %d: channel_norm = %e\n", i, H_combined.row(i).norm());
		}
	}

	// Use Regularized ZF (RZF) / MMSE precoding for better stability
	// RZF: W = H^H (H H^H + αI)^-1
	Real alpha = 1e-6;  // Regularization parameter
	if (std::abs(det) < EPSILON) {
		alpha = 1e-4;  // Increase regularization for ill-conditioned channels
		printf("[INFO] TDD: Using increased regularization alpha=%e\n", alpha);
	}

	MatrixXcReal regularization = alpha * MatrixXcReal::Identity(num_valid_ue, num_valid_ue);
	MatrixXcReal H_HH_reg = H_HH + regularization;

	// Compute precoding matrix: W = H^H (H H^H + αI)^-1
	W[rb_idx] = H_combined.adjoint() * H_HH_reg.inverse();

	// Power normalization: normalize each column to have unit norm
	VectorXReal sqrt_p = VectorXReal::Zero(num_valid_ue);

	for (int i = 0; i < num_valid_ue; i++)
	{
		Real col_norm = W[rb_idx].col(i).norm();
		if (col_norm < EPSILON || std::isnan(col_norm) || std::isinf(col_norm)) {
			printf("[ERROR] TDD: W[%d].col(%d).norm() = %e is invalid\n", rb_idx, i, col_norm);
			sqrt_p(i) = 1.0;  // Fallback value
		} else {
			sqrt_p(i) = 1.0 / col_norm;
		}
	}

	W[rb_idx] = W[rb_idx] * sqrt_p.asDiagonal();

	// Debug: Check for NaN in final precoding matrix
	bool has_nan = false;
	for (int i = 0; i < W[rb_idx].rows(); i++) {
		for (int j = 0; j < W[rb_idx].cols(); j++) {
			if (std::isnan(W[rb_idx](i,j).real()) || std::isnan(W[rb_idx](i,j).imag())) {
				printf("[ERROR] TDD: FINAL W[%d](%d,%d) has NaN\n", rb_idx, i, j);
				has_nan = true;
			}
		}
	}

	if (!has_nan && std::abs(det) < EPSILON) {
		printf("[SUCCESS] TDD: Precoding matrix computed successfully with regularization\n");
	}
}


/*===================================================================
FUNCTION: Sector::Set_Min_Cqi( void )

AUTHOR: SJ, KR, DW

DESCRIPTION: Check the cqi for every allocated resource block, and Set it with minimum Cqi.

NOTES: Updated 6 OCT 08
===================================================================*/
void Sector::Set_Min_Cqi(void)
{
	for (int ue_idx = 0; ue_idx < (int)ue_in_control.size(); ue_idx++)
	{
		int min_cqi = NUM_CQI_MCS_TYPES;
		// MIN CQI within UE setting for each frequency
		for (int freq_idx = 0; freq_idx < num_rb; freq_idx++)
		{
			for (int stream = 0; stream < mx_ue_mumimo; stream++)
			if (ppschedulewrite[freq_idx][stream].ue_selected == ue_in_control[ue_idx])
			{
				if (min_cqi > MCS_decision[ue_idx][freq_idx])
				{
					min_cqi = MCS_decision[ue_idx][freq_idx];
				}
			}
		}

		// Set MIN CQI Scheduling for each frequency
		for (int freq_idx = 0; freq_idx < num_rb; freq_idx++)
		{
			for (int stream = 0; stream < mx_ue_mumimo; stream++)
			if (ppschedulewrite[freq_idx][stream].ue_selected == ue_in_control[ue_idx])
			{
				ppschedulewrite[freq_idx][stream].mcs_selected = min_cqi;
			}
		}
	}
}



void Sector::Scheduling_algorithm_module_MU_MIMO(void)
{
	/*===================================================================
	Subband-based User selection algorithm with Multi-threading

	Process each subband independently in parallel to exploit:
	1. Frequency-selective fading (different UEs optimal at different frequencies)
	2. Multi-core CPU resources (parallel scheduling computation)

	OUTPUT : scheduled_ue (per subband, stored in ppschedulewrite)
	===================================================================*/

	// Subband configuration
	const int subband_size = 16;  // Number of RBs per subband (typical for 5G NR)
	int num_subband = (int)(num_rb / subband_size);
	if (num_rb % subband_size != 0) {
		num_subband += 1;  // Add one more subband for remaining RBs
	}

	// Process each subband independently with OpenMP parallel for
	// Dynamic scheduling ensures load balancing (last subband may be smaller)
	#if ENABLE_MULTITHREADING
	#pragma omp parallel for schedule(dynamic)
	#endif
	for (int subband_idx = 0; subband_idx < num_subband; subband_idx++)
	{
		// Thread-local variables for independent scheduling state
		vector<int> local_selected_ue;
		vector<int> local_scheduled_ue;
		vector<int> local_user_set;

		// Determine RB range for this subband
		int rb_start = subband_idx * subband_size;
		int rb_end = rb_start + subband_size;

		// Handle the last subband (may have fewer RBs)
		if (rb_end > num_rb) {
			rb_end = num_rb;
		}

		// Select a representative RB in this subband for scheduling decision
		// (typically the middle RB to get average channel condition)
		int schedule_RB = rb_start + (rb_end - rb_start) / 2;

		// Perform user selection for this subband based on selected algorithm
		if (g_mumimo_scheduling_algorithm == 1) {
			// Chordal Distance based scheduling
			ChordalDistance_Scheduling_threadsafe(schedule_RB, local_selected_ue, local_scheduled_ue);
		} else {
			// Default: SUS (Semi-Orthogonal User Selection)
			first_max_ue_threadsafe(schedule_RB, local_selected_ue);

			// OPTIMIZATION: TDD and FDD now use the SAME optimized code path!
			// TDD: PMI_vector_read contains pre-computed dominant eigenvectors (from Store_CSI_for_TDD)
			// FDD: PMI_vector_read contains quantized PMI vectors (from Quantization_of_Ch)
			// Both are already normalized in MS, so no eigenvector computation needed here!
			SUS_procedure_threadsafe(schedule_RB, local_selected_ue, local_scheduled_ue, local_user_set);
		}

		/*===================================================================
		SU-vs-MU decision (g_su_fallback=1): compare best single-UE SU metric
		against the MU group metric. If SU wins, override local_scheduled_ue
		with just that single UE.

		SU metric per UE u: M_SU(u) = self_RI(u) · METRIC[u]
			(dominant-layer PF metric times its optimal rank approximates
			 the per-UE rate at R layers, divided by long-term avg rate)
		MU metric: sum of METRIC[u] across local_scheduled_ue.
		===================================================================*/
		if (g_su_fallback == 1 && g_codebook_type == 3 && TDD_mode == 0) {
			// MU group metric (current selection)
			Real mu_metric_sum = 0.0;
			for (int ui = 0; ui < (int)local_scheduled_ue.size(); ui++) {
				int sel = local_scheduled_ue[ui];
				if (sel < 0) continue;
				mu_metric_sum += METRIC[sel][schedule_RB];
			}

			// Best SU candidate across ue_in_control
			Real best_su_metric = 0.0;
			int  best_su_ue    = -1;
			for (int ue_idx = 0; ue_idx < (int)ue_in_control.size(); ue_idx++) {
				int ms_idx = ue_in_control[ue_idx];
				int R_ue = (g_rank_adaptive == 1) ? ms[ms_idx].self_RI : g_type2_rank;
				if (R_ue < 1) R_ue = 1;
				Real su_m = (Real)R_ue * METRIC[ue_idx][schedule_RB];
				if (su_m > best_su_metric) {
					best_su_metric = su_m;
					best_su_ue    = ue_idx;
				}
			}

			if (best_su_ue >= 0 && best_su_metric > mu_metric_sum) {
				// Override: schedule the best single UE alone
				local_scheduled_ue.clear();
				local_scheduled_ue.push_back(best_su_ue);
			}
		}

		/*===================================================================
		Write selected UEs to schedule map for all RBs in this subband.
		Each UE occupies R_ue consecutive streams (R_ue = ms[u].self_RI when
		rank-adaptive, else g_type2_rank for eType II, else 1).
		Stream budget capped at mx_ue_mumimo.
		OUTPUT : ppschedulewrite[rb_idx][stream].ue_selected / .layer_idx
		===================================================================*/
		for (int rb_idx = rb_start; rb_idx < rb_end; rb_idx++)
		{
			// Initialize schedule for this RB
			for (int stream = 0; stream < mx_ue_mumimo; stream++)
			{
				ppschedulewrite[rb_idx][stream].ue_selected = -1;
				ppschedulewrite[rb_idx][stream].layer_idx   = 0;
			}

			// Write selected UEs — each UE occupies R_ue consecutive slots
			int stream_w = 0;
			for (int ue_idx = 0; ue_idx < (int)local_scheduled_ue.size(); ue_idx++)
			{
				if (local_scheduled_ue[ue_idx] < 0) continue;
				int ms_idx = ue_in_control[local_scheduled_ue[ue_idx]];

				int R_ue;
				if (TDD_mode == 0 && g_codebook_type == 3) {
					R_ue = (g_rank_adaptive == 1) ? ms[ms_idx].self_RI : g_type2_rank;
				} else {
					R_ue = 1;
				}
				if (R_ue < 1) R_ue = 1;

				for (int layer = 0; layer < R_ue && stream_w < mx_ue_mumimo; layer++, stream_w++)
				{
					ppschedulewrite[rb_idx][stream_w].ue_selected  = ms_idx;
					ppschedulewrite[rb_idx][stream_w].layer_idx    = layer;
					ppschedulewrite[rb_idx][stream_w].mcs_selected = MCS_decision[local_scheduled_ue[ue_idx]][rb_idx];
					ppschedulewrite[rb_idx][stream_w].cqi_selected = CQIndex_decision[local_scheduled_ue[ue_idx]][rb_idx];
				}
				if (stream_w >= mx_ue_mumimo) break;
			}

			// TDD mode: Use channel-based precoding, FDD mode: Use PMI-based precoding
			if (TDD_mode == 1) {
				Transmit_Precoding_TDD(rb_idx, local_scheduled_ue);
			} else {
				Transmit_Precoding(rb_idx, local_scheduled_ue);
			}
		}
	}
}


void Sector::AVR_Cqi_update(void)   ////
{
	for (int ue_idx = 0; ue_idx < ue_in_control.size(); ue_idx++)
	{
		for (int rb_idx = 0; rb_idx < num_rb; rb_idx++)
		{
			if ((TYPE == 12) && (Configuration_Type == 1 || Configuration_Type == 3)) // 
			{
				bool unscheduled = 1;
				ms[ue_in_control[ue_idx]].unscheduled_stack += 1;

				for (int stream = 0; stream < mx_ue_mumimo; stream++)
				{
					if (ppschedulewrite[rb_idx][stream].ue_selected == ue_in_control[ue_idx])
					{
						unscheduled = 0;
						ms[ue_in_control[ue_idx]].unscheduled_stack -= 1;
					}
				}

				if (CQI_AVR[ue_idx][rb_idx] && unscheduled)
				{
					CQI_AVR[ue_idx][rb_idx] -= (1. / N_pf) * log2(1 + CQI_read[ue_idx][rb_idx]);
				}
			}
			else
			{

				if (CQI_AVR[ue_idx][rb_idx])
				{
					CQI_AVR[ue_idx][rb_idx] -= (1. / _N_pf) * log2(1 + CQI_read[ue_idx][rb_idx]);
				}
			}
		}
	}
}


void Sector::Scheduling_algorithm_module_RoundRobin(void)
{
	/*===================================================================
	User selection algorithm

	OUTPUT : selected_ue
	===================================================================*/
	int selected_ue = (int)(ue_in_control.size() * randnum.u());


	/*===================================================================
	Write selected ue to schedule map

	OUTPUT : ppschedulewrite[rb_idx][stream].ue_selected
	===================================================================*/
	for (int rb_idx = 0; rb_idx < num_rb; rb_idx++)
	{
		for (int stream = 0; stream < mx_ue_mumimo; stream++)
		{

			ppschedulewrite[rb_idx][stream].ue_selected = ue_in_control[selected_ue];

		}

	}
}


/*===================================================================
FUNCTION: Sector::Compute_Chordal_Distance

DESCRIPTION:
  Computes the Chordal Distance between the projection matrix of
  previously selected users and a candidate user's channel.

  Chordal Distance formula (from paper):
  dc(A, B) = (1/sqrt(2)) * ||A_o * A_o^H - B_o * B_o^H||_F

  where A_o and B_o are orthonormal bases for subspaces A and B.

INPUT:
  P_projection - Projection matrix of previously selected channels (P_o * P_o^H)
  H_candidate  - Candidate user's subspace basis source, N_tx x k (columns span
                 the candidate subspace). FDD: PMI precoder (N_tx x R). TDD:
                 the caller passes the transposed channel (N_tx x N_rx).
                 NOTE: callers must hand in the N_tx-row form; this function does
                 NOT transpose internally (keeps it consistent with P_projection,
                 which is also built from N_tx-row bases).

OUTPUT:
  Returns the Chordal Distance (Real)
===================================================================*/
Real Sector::Compute_Chordal_Distance(const MatrixXcReal& P_projection, const MatrixXcReal& H_candidate)
{
	// QR decomposition to get an orthonormal basis for the candidate subspace.
	// H_candidate is already N_tx x k (columns = subspace), so use it directly.
	Eigen::HouseholderQR<MatrixXcReal> qr(H_candidate);
	int rank = min((int)H_candidate.rows(), (int)H_candidate.cols());
	MatrixXcReal Q = qr.householderQ() * MatrixXcReal::Identity(H_candidate.rows(), rank);

	// Compute projection matrix for candidate: G_{k,o} * G_{k,o}^H
	MatrixXcReal G_projection = Q * Q.adjoint();

	// Compute Chordal Distance: (1/sqrt(2)) * ||P_projection - G_projection||_F
	MatrixXcReal diff = P_projection - G_projection;
	Real chordal_dist = diff.norm() / sqrt(2.0);

	return chordal_dist;
}


/*===================================================================
FUNCTION: Sector::ChordalDistance_Scheduling_threadsafe

DESCRIPTION:
  Chordal Distance based MU-MIMO user selection algorithm.

  Algorithm (from IEEE TCOM 2012 paper by Ko & Lee):

  Step I: Initialization
    - Select user with maximum Frobenius norm as first user

  Step II: Loop (i = 2 to K)
    - Update orthonormal basis of selected channels
    - For each candidate: metric = C_k^alpha * ||P_projection - G_k_projection||_F
    - Select user with maximum metric
    - Early termination if capacity doesn't improve

INPUT:
  rb_idx - Resource block index for scheduling

OUTPUT:
  local_selected_ue  - Indices of selected users (in ue_in_control)
  local_scheduled_ue - Final scheduled users
===================================================================*/
void Sector::ChordalDistance_Scheduling_threadsafe(int rb_idx, vector<int>& local_selected_ue,
                                                    vector<int>& local_scheduled_ue)
{
	local_selected_ue.clear();
	local_scheduled_ue.clear();

	int num_ue = ue_in_control.size();
	if (num_ue == 0) return;

	// Safety checks
	if (TDD_mode == 1 && CSI_matrix_read == NULL) return;
	if (TDD_mode == 0 && PMI_vector_read == NULL) return;

	// ===================================================================
	// Step I: Initialization - Select user with maximum PF-weighted metric
	// PF metric: METRIC[ue_idx][rb_idx] = instantaneous_rate^alpha / average_rate^beta
	// Channel quality: CQI_read[ue_idx][rb_idx] = SINR (includes pathloss, interference, noise)
	// Combined metric: PF_metric * SINR^chordal_alpha
	// ===================================================================
	vector<Real> channel_quality(num_ue);  // SINR-based channel quality
	vector<bool> available(num_ue, true);

	Real max_metric = -1.0;
	int first_user = -1;

	for (int ue_idx = 0; ue_idx < num_ue; ue_idx++) {
		if (ue_idx >= 80) break;  // Bounds check

		// Use SINR (CQI_read) which already includes:
		// - Large scale fading (pathloss, shadowing)
		// - Small scale fading (channel coefficients)
		// - Inter-cell interference
		// - Noise
		// This is more accurate than raw Frobenius norm
		Real sinr_linear = CQI_read[ue_idx][rb_idx];
		if (sinr_linear < 1e-10) sinr_linear = 1e-10;  // Avoid zero/negative

		channel_quality[ue_idx] = sinr_linear;

		// PF-weighted initial selection metric:
		// PF_metric * SINR^alpha (channel quality weighted by PF fairness)
		Real pf_metric = METRIC[ue_idx][rb_idx];
		Real combined_metric = pf_metric * pow(sinr_linear, g_chordal_alpha);

		if (combined_metric > max_metric) {
			max_metric = combined_metric;
			first_user = ue_idx;
		}
	}

	if (first_user < 0) return;

	// Select first user
	local_selected_ue.push_back(first_user);
	available[first_user] = false;

	// Get channel dimension for projection matrix size
	int N_tx = NUM_TX_Port;  // Number of transmit antennas

	// Initialize combined channel matrix for orthonormalization
	// P_o starts with first user's channel (orthonormalized)
	MatrixXcReal P_combined;

	if (TDD_mode == 1) {
		MatrixXcReal H_first = CSI_matrix_read[first_user][rb_idx];
		P_combined = H_first.transpose();  // N_tx x N_rx
	} else {
		// For multi-layer PMI, use all R layers as the subspace representation
		P_combined = PMI_vector_read[first_user][rb_idx];  // N_tx x R
	}

	// ===================================================================
	// Step II: Iteratively select users with maximum Chordal Distance metric
	// ===================================================================
	while ((int)local_selected_ue.size() < mx_ue_mumimo && (int)local_selected_ue.size() < num_ue)
	{
		// Orthonormalize combined channel matrix using QR decomposition
		Eigen::HouseholderQR<MatrixXcReal> qr(P_combined);
		int rank = min((int)P_combined.rows(), (int)P_combined.cols());
		MatrixXcReal P_o = qr.householderQ() * MatrixXcReal::Identity(P_combined.rows(), rank);

		// Compute projection matrix: A = P_o * P_o^H
		MatrixXcReal P_projection = P_o * P_o.adjoint();

		// Find user with maximum metric: PF_metric * C_k^alpha * chordal_distance
		// This combines Proportional Fairness with Chordal Distance based selection
		Real max_combined_metric = -1.0;
		int best_user = -1;

		for (int ue_idx = 0; ue_idx < num_ue; ue_idx++) {
			if (!available[ue_idx]) continue;
			if (ue_idx >= 80) continue;

			// Get candidate channel
			MatrixXcReal H_candidate;
			if (TDD_mode == 1) {
				// Transpose to N_tx x N_rx so it matches the N_tx-row form that
				// Compute_Chordal_Distance and P_combined expect.
				H_candidate = CSI_matrix_read[ue_idx][rb_idx].transpose();
			} else {
				H_candidate = PMI_vector_read[ue_idx][rb_idx];  // N_tx x R
			}

			// Compute Chordal Distance
			Real chordal_dist = Compute_Chordal_Distance(P_projection, H_candidate);

			// Compute PF-weighted Chordal Distance metric:
			// PF_metric * SINR^alpha * chordal_distance
			// - PF_metric ensures fairness among users (prevents starvation)
			// - SINR^alpha (channel quality including interference) ensures good link quality
			// - chordal_distance ensures spatial compatibility (orthogonality)
			Real pf_metric = METRIC[ue_idx][rb_idx];
			Real metric = pf_metric * pow(channel_quality[ue_idx], g_chordal_alpha) * chordal_dist;

			if (metric > max_combined_metric) {
				max_combined_metric = metric;
				best_user = ue_idx;
			}
		}

		// If no valid user found, stop
		if (best_user < 0) break;

		// Add selected user
		local_selected_ue.push_back(best_user);
		available[best_user] = false;

		// Update combined channel matrix by adding new user's channel
		if (TDD_mode == 1) {
			MatrixXcReal H_new = CSI_matrix_read[best_user][rb_idx];
			MatrixXcReal H_new_T = H_new.transpose();  // N_tx x N_rx

			// Horizontally concatenate
			MatrixXcReal P_new(P_combined.rows(), P_combined.cols() + H_new_T.cols());
			P_new << P_combined, H_new_T;
			P_combined = P_new;
		} else {
			MatrixXcReal v_new = PMI_vector_read[best_user][rb_idx];  // N_tx x R

			// Horizontally concatenate
			MatrixXcReal P_new(P_combined.rows(), P_combined.cols() + v_new.cols());
			P_new << P_combined, v_new;
			P_combined = P_new;
		}
	}

	// Copy to scheduled_ue
	local_scheduled_ue = local_selected_ue;
}


#if 0  // unused: ChordalDistance_Scheduling — replaced by ChordalDistance_Scheduling_threadsafe
/*===================================================================
FUNCTION: Sector::ChordalDistance_Scheduling (non-threadsafe version)

DESCRIPTION:
  Wrapper for non-threadsafe context. Uses class member variables.
===================================================================*/
void Sector::ChordalDistance_Scheduling(int rb_idx)
{
	ChordalDistance_Scheduling_threadsafe(rb_idx, selected_ue, scheduled_ue);
}
#endif  // unused: ChordalDistance_Scheduling