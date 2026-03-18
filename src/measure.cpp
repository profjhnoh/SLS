#include <iostream>
#include <fstream>

using namespace std;

#include <cstring>
#include <cstdlib>

#include "common.h"
#include "const.h"

Real FER(Real SINR, int MCS);
void Averaging_Sinr(int ue_idx, Real sinr);

typedef pair<int,unsigned long int> throughput_pair;
bool comparator_throughput ( const throughput_pair& l, const throughput_pair& r) 
{
    return l.second < r.second;
}

void Measure()
{
	#if ENABLE_MULTITHREADING
	#pragma omp parallel num_threads(num_of_threads)
	{
		#pragma omp for
	#endif
		for (int idx = 0; idx < num_MS; idx++)
		{
			ms[idx].Receive_DL_mTRP();
		}
	#if ENABLE_MULTITHREADING
	}
	#endif

	throughput_pair * per_ue_thru = new throughput_pair[num_MS];
	for (int idx = 0; idx < num_MS; idx++)
	{
		if (ms[idx].rb_indices.size() > 0) {
			if (ms[idx].received_sinr_avg < 0)
				cout << " something wrong with SINR 3" << endl;

			Averaging_Sinr(idx, ms[idx].received_sinr_avg ); 
		}
	}

	
	if ( t > 50 )
	{
		for (int idx = 0; idx < num_MS; idx++)
		{
			per_ue_thru[idx].first  = idx;
			per_ue_thru[idx].second = ms[idx].Return_Throughput(idx);
		}

		sort(per_ue_thru, per_ue_thru + num_MS, comparator_throughput);
		for (int idx = 0; idx < num_MS; idx++)
		{
			int ue_idx   = per_ue_thru[idx].first;
			int sec_idx  = links[ue_idx]._sector_in_control;
			auto it = find(sector_indices_schedule.begin(),
			                sector_indices_schedule.end(),
							sec_idx);
			if ( it == sector_indices_schedule.end() )
				sector_indices_schedule.push_back(sec_idx);
		}		

	}
	//cout << " Receive_DL... " << endl;
}

#if 0  // unused: FER — never called
Real FER(Real SINR, int MCS)
{
	Real Frame_Error_Rate = 1;
	Real ESINR_L = SINR; // linear value
	switch (MCS)
	{
	case -1:
		Frame_Error_Rate = 0.0;
		break;

	case 0:
		if (ESINR_L < pow(10.0, -7.737 / 10.0))
			Frame_Error_Rate = 1.0;
		else
			Frame_Error_Rate = exp(69.109 - 410.4 * ESINR_L);
		break;
	case 1:
		if (ESINR_L < pow(10.0, -5.852 / 10.0))
			Frame_Error_Rate = 1.0;
		else
			Frame_Error_Rate = exp(70.072 - 269.63 * ESINR_L);

		break;

	case 2:
		if (ESINR_L < pow(10.0, -3.737 / 10.0))
			Frame_Error_Rate = 1.0;
		else
			Frame_Error_Rate = exp(69.19 - 163.42 * ESINR_L);

		break;

	case 3:
		if (ESINR_L < pow(10.0, -1.718 / 10.0))
			Frame_Error_Rate = 1.0;
		else
			Frame_Error_Rate = exp(70.491 - 104.69 * ESINR_L);

		break;

	case 4:
		if (ESINR_L < pow(10.0, 0.3206 / 10.0))
			Frame_Error_Rate = 1.0;
		else
			Frame_Error_Rate = exp(77.142 - 71.651* ESINR_L);

		break;

	case 5:
		if (ESINR_L < pow(10.0, 2.14 / 10.0))
			Frame_Error_Rate = 1.0;
		else
			Frame_Error_Rate = exp(47.496 - 29.018 * ESINR_L);

		break;

	case 6:
		if (ESINR_L < pow(10.0, 4.096 / 10.0))
			Frame_Error_Rate = 1.0;
		else
			Frame_Error_Rate = exp(45.163 - 17.585 * ESINR_L);

		break;

	case 7:
		if (ESINR_L < pow(10.0, 6.052 / 10.0))
			Frame_Error_Rate = 1.0;
		else
			Frame_Error_Rate = exp(45.936 - 11.402 * ESINR_L);

		break;

	case 8:
		if (ESINR_L < pow(10.0, 8.03 / 10.0))
			Frame_Error_Rate = 1.0;
		else
			Frame_Error_Rate = exp(44.625 - 7.0239 * ESINR_L);

		break;

	case 9:
		if (ESINR_L < pow(10.0, 10.03 / 10.0))
			Frame_Error_Rate = 1.0;
		else
			Frame_Error_Rate = exp(46.766 - 4.6398 * ESINR_L);

		break;

	case 10:
		if (ESINR_L < pow(10.0, 11.82 / 10.0))
			Frame_Error_Rate = 1.0;
		else
			Frame_Error_Rate = exp(41.834 - 2.7525 * ESINR_L);

		break;
	case 11:
		if (ESINR_L < pow(10.0, 13.65 / 10.0))
			Frame_Error_Rate = 1.0;
		else
			Frame_Error_Rate = exp(33.542 - 1.4485 * ESINR_L);

		break;
	case 12:
		if (ESINR_L < pow(10.0, 15.69 / 10.0))
			Frame_Error_Rate = 1.0;
		else
			Frame_Error_Rate = exp(29.392 - 0.79247 * ESINR_L);

		break;
	case 13:
		if (ESINR_L < pow(10.0, 17.5 / 10.0))
			Frame_Error_Rate = 1.0;
		else
			Frame_Error_Rate = exp(29.525 - 0.52504 * ESINR_L);

		break;
	case 14:
		if (ESINR_L < pow(10.0, 18.98 / 10.0))
			Frame_Error_Rate = 1.0;
		else
			Frame_Error_Rate = exp(8.154 - 0.10315 * ESINR_L);

		break;
	}
	return Frame_Error_Rate;
}
#endif  // unused: FER