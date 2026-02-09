#include "const.h"
#include "common.h"
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

#include <cstring>
#include <cstdlib>


void MS::Configuration(int _self_idx,int _self_sector_idx)
{
	self_idx = _self_idx;
	self_sector_idx = _self_sector_idx;

	total_esinr  = 0;
	sum_cqi      = 0;
	sum_mcs_type = 0;

	links[_self_idx].total_num_scheduled = 0;
	links[_self_idx].total_num_tx = 0;
	links[_self_idx].failed = 0;
	links[_self_idx].total_num_re_tx = 0;
	links[_self_idx].re_tx_failed = 0;

	// Cache noise value in linear scale (constant for entire simulation)
	noise_linear_cached = dBm2linear(thermal_noise + 10.0 * log10(bandwidth) + MS_noisefig);

	Set_MS_Location();
}

void MS::Set_MS_Location()
{
	//////////// Hex grid MS configuration //////////////////////////
	/*
	INPUT:
	OUTPUT: MS location
	Notes: 
	*/

	Real angle;
	Real R;
	Real cell_radius = inter_site_distance * (1 / sqrt(3.));

	do
	{
		Real x_tmp = cell_radius*randnum.u();
		Real y_tmp = cell_radius*sqrt(3.)/2.*randnum.u();

		if ( y_tmp > -sqrt(3)*(x_tmp-cell_radius))
		{
			Real y_tmp2 = cell_radius - x_tmp;
			Real x_tmp2 = y_tmp;

			loc.x = x_tmp2*cos(-pi/6) - y_tmp2*sin(-pi/6);
			loc.y = x_tmp2*sin(-pi/6) + y_tmp2*cos(-pi/6);
		} 
		else 
		{
			loc.x = x_tmp;
			loc.y = y_tmp;
		}
	}
	while (sqrt(loc.x * loc.x + loc.y * loc.y) < min_distance);

	angle = atan(loc.y / loc.x) * 180./pi;
	R = sqrt(loc.x*loc.x + loc.y*loc.y);

	//////////////////////////////////////////// 3 Sector /////////////
	if (self_sector_idx == 1)
	{
		angle = angle + 120;
		loc.x = R*cos(angle*pi / 180.);
		loc.y = R*sin(angle*pi / 180.);
	}
	else if (self_sector_idx == 2)
	{
		angle = angle + 240;
		loc.x = R*cos(angle*pi / 180.);
		loc.y = R*sin(angle*pi / 180.);
	}

	////////////////////////////////////////////////// Indoor or Outdoor /////////////////////////////////
	
	if (TYPE == 12)  //// IMT-2020 , Dense Urban eMBB
	{
		// 230721 by jhnoh 
		// Assume all outdoor user
		// 251011 by jhnoh 
		// Recover the original code for indoor user
		
		if (randnum.u() < 0.8) // indoor
		{
			Indoor = true;
		}
		else
		{
			Indoor = false;
		}
	
		//Indoor = false;
		if (carrier_freq > 20000000000.)   /// 30GHz scenario -> all outdoor
		{
			Indoor = false;
		}
		else
		{
			if (randnum.u() < 0.8) // indoor
			{
				Indoor = true;
			}
			else
			{
				Indoor = false;
			}
		}
	}
	else if (TYPE == 13)  //// IMT-2020 , Rural eMBB
	{
		if (Configuration_Type == 2)    ///// configuration C -> 40% indoor, 20% outdoor pedestrian, 20% outdoor in car
		{
			if (randnum.u() < 0.4) // indoor 40%
			{
				Indoor = true;
			}
			else
			{
				Indoor = false;
			}
		}
		else
		{
			if (randnum.u() < 0.5) // indoor
			{
				Indoor = true;
			}
			else
			{
				Indoor = false;
			}
		}
		//Indoor = false;  // 251011: Commented out - this was overriding the indoor/outdoor setting above
	}
	else
	{
		Indoor = false;
	}
}

void MS::Allocate_memory()
{
	TLSPs  = new ArrayXReal[num_BS + num_mTRP];
	LSPs   = new ArrayXReal[num_BS + num_mTRP];
}

void MS::Delete_memory()
{
	delete [] TLSPs;
	delete []  LSPs;
}

void MS::Reset2Default()
{
		unscheduled_stack = 0;
		self_idx          = 0;
		self_sector_idx   = 0;
		nearest_bs_idx    = 0;
		loc               = {0,0};
		//d_rx[2][4][2][1][2];//  antenna element location vector [M][N][P][Mg][Ng] (max num)
		for (int i_0 = 0; i_0 < 2; i_0++)
			for (int i_1 = 0; i_1 < 4; i_1++)
				for (int i_2 = 0; i_2 < 2; i_2++)
					for (int i_3 = 0; i_3 < 1; i_3++)
						for (int i_4 = 0; i_4 < 2; i_4++)
							d_rx[i_0][i_1][i_2][i_3][i_4] = {0,0,0};


		// array antenna rotate angle
		alpha             = 0;
		beta              = 0;
		gamma             = 0;

		speed             = 0;
		Indoor            = false;

		//receive_downlink
		receive_self_ms_idx = 0;
		received_sinr_avg   = 0;
		ESINR_linear        = 0;
		sum_cqi             = 0;
		sum_mcs_type        = 0;
		total_esinr         = 0;
		total_estimate_sinr = 0;
		_info_bits          = 0;	
		_mcs_idx            = 0;
		_cqi_idx            = 0;
		_avr_sinr           = 0;
		_mod_type           = 0;
		num_rx_rb           = 0;
		_num_traffic        = 0;
		_code_rate          = 0;
		_snr_scale_factor   = 0;
 		num_scheduled       = 0;
		num_comp_tx         = 0;
			
		// cqi_offset
		CQI_offset          = 0;
		num_added           = 0;
		ue_BLER_Value       = 0;

		//int *propagation_condition;
		floor_idx           = 0;
		ms_d_in             = 0;
		MS_HEIGHT_FINAL     = 0;
		pos_indoor          = {0,0};

		comp_mode_rx_flag   = 0;
		num_rx_rbs          = 0;
		re_tx_failed_cnt    = 0;
		cqi_offset          = 0;

		// OLLA: Initialize with configured initial offset and moving window
		olla_offset         = g_olla_init_offset;
		olla_history.clear();
		olla_history.resize(g_olla_window_size, true);  // Initialize with all ACKs (optimistic start)
		olla_history_idx    = 0;
		olla_nack_count     = 0;  // All ACKs initially, so NACK count is 0

		pmi_l               = 0;
		pmi_m               = 0;
		pmi_n               = 0;

		//CQI             = NULL;
		//PMI             = NULL;
		//PMI_vector      = NULL;

		H_m             = NULL;

		ppSchedulerRead = NULL;

		for(int bs_idx = 0; bs_idx < num_BS; bs_idx++) 
		{
			TLSPs[bs_idx] = ArrayXReal::Zero(0);
			LSPs [bs_idx] = ArrayXReal::Zero(0);
		}

		rbs_rx.clear(); 
		rb_indices.clear();
}


/*===================================================================
FUNCTION: UE::return_throughput( void )



===================================================================*/

Real MS::Return_Throughput(int _ms_idx) const
{
	return (links[_ms_idx]._throughput);
}

void MS::Kill_Memory(void)
{
	
#if 0
	delete[] H_m;


	for (int i = 0; i < num_compute_coef - 1; i++)
	{
		delete[] HI_m[i];
	}

	delete[] HI_m;
#endif

}
