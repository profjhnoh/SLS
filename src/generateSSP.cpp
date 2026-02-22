#include "common.h"


#include <cstring>
#include <cstdlib>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

void Generate_SSP();
void Compute_Channel_Coef();
void Update_Channel_Coef();
Real Get_distance(LOCATION, LOCATION);


void Generate_SSP()
{
	
	cout << "Generate SSP......................"<<std::flush;// << endl;

	if (TYPE == 11)   //// Get all BS-MS channel coefficient
	{
		for (int bs_idx = 0; bs_idx < num_BS; bs_idx++)
		{
			for (int ms_idx = 0; ms_idx < num_MS; ms_idx++)
			{
				channel[bs_idx][ms_idx].Set_SmallScaleParameter(bs_idx, ms_idx);
			}
		}
	}
	else           ////// Get MS - adjacent_BS channel coefficient
	{
		for (int bs_idx = 0; bs_idx < num_BS; bs_idx++)
		{
			for (int ms_idx = 0; ms_idx < num_MS; ms_idx++)
			{
				channel[bs_idx][ms_idx].Set_SmallScaleParameter(bs_idx, ms_idx);
			}
		}
	}
	cout << "DONE" << endl;
}



void Compute_Channel_Coef()
{
	if (CH_CAL == 1)  /// use runtime calculate
	{
		cout << "Channel Matrix...................."<<std::flush;

		if (TYPE == 11 && num_Indoor_TRxP == 1)   //// Get all BS-MS channel coefficient
		{
			for (int ms_idx = 0; ms_idx < num_MS; ms_idx++)
			{
				int adj_sector_num_to_BS = (int)(links[ms_idx]._sector_in_control);
				int sec_idx = adj_sector_num_to_BS;  // InH 1TRxP: bs_idx == sector
				channel[adj_sector_num_to_BS][ms_idx].GetChannelImpulseResponse(adj_sector_num_to_BS, ms_idx, sec_idx);
			}
		}
		else  // Get MS - adjacent_BS channel coefficient
		{
			for (int coeff_idx = 0; coeff_idx < num_compute_coef; coeff_idx++)
			{
				#if ENABLE_MULTITHREADING
				#pragma omp parallel num_threads(num_of_threads)
				{
					#pragma omp for
				#endif
					for (int ms_idx = 0; ms_idx < num_MS; ms_idx++)
					{
						int adj_sector_num_to_BS;
						int adj_sector;
						adj_sector_num_to_BS = (int)(links[ms_idx].adj_sector[coeff_idx] / 3);
						adj_sector = links[ms_idx].adj_sector[coeff_idx];
						int sec_idx = adj_sector % 3;
						channel[adj_sector_num_to_BS][ms_idx].GetChannelImpulseResponse(adj_sector_num_to_BS, ms_idx, sec_idx);
					}
				#if ENABLE_MULTITHREADING
				}
				#endif
			}
		}
		cout <<"DONE"<<endl<<endl;;
	}
	else
	{
		//
	}

}

void Update_Channel_Coef()
{
	if (CH_CAL == 1)  /// use runtime calculate
	{
		if (TYPE == 11 && num_Indoor_TRxP == 1)   //// Get all BS-MS channel coefficient
		{
			for (int ms_idx = 0; ms_idx < num_MS; ms_idx++)
			{
				int adj_sector_num_to_BS = (int)(links[ms_idx]._sector_in_control);
				int sec_idx = adj_sector_num_to_BS;
				channel[adj_sector_num_to_BS][ms_idx].Update_H_usn_per_time(t, ms_idx, sec_idx);
			}
		}
		else  ////// Get MS - adjacent_BS channel coefficient
		{
			for (int coeff_idx = 0; coeff_idx < num_compute_coef; coeff_idx++)
			{
				#if ENABLE_MULTITHREADING
				#pragma omp parallel num_threads(num_of_threads)
				{
					#pragma omp for
				#endif
					for (int ms_idx = 0; ms_idx < num_MS; ms_idx++)
					{
						int adj_sector_num_to_BS = (int)(links[ms_idx].adj_sector[coeff_idx] / 3);
						int adj_sector = links[ms_idx].adj_sector[coeff_idx];
						int sec_idx = adj_sector % 3;
						channel[adj_sector_num_to_BS][ms_idx].Update_H_usn_per_time(t, ms_idx, sec_idx);
					}
				#if ENABLE_MULTITHREADING
				}
				#endif
			}
		}
	}
	else
	{
		//
	}

}