
/*
////////////////////////
5G DownLink System Level Simulator
////////////////////////
*/

#include <stdio.h>
#include <time.h>
#include <omp.h>
#include <thread>
#include <chrono>

#include "const.h"
#include "common.h"
#include <ctime>

//#include "indicators.hpp"

// Global variable definition for ray-level Doppler precision mode
int USE_RAY_LEVEL_DOPPLER = 0;  // Default: use cluster-average Doppler (fast)

// Global variable definitions for precoding-based SINR calculation
int USE_PRECODING_BASED_SINR = 0;  // Default: use feedback CQI
Real INTERCELL_INTERFERENCE_MARGIN_DB = 0.0;  // Default: 0 dB margin (assumes perfect inter-cell interference knowledge)

void Set_simul_param(int argc, char *argv[]);
void Set_Parameter(int scenario);
void StandardInitialization();
void InitializeSystem();
void Set_fastfading_param();
void Allocate_memory(void);
void Initialdrop();
void Initialize_CHANNEL();
void Print_Location(int drop_idx);
void Parameter_initialization();
void Generate_SSP();
void Compute_Channel_Coef();
void Update_Channel_Coef();
void ClearMAPs2StartADrop();
void PerDropStatistics(int);
void Measure();
void Generate_LSP(void);
void Link_configuration();
void Get_CouplingLoss();
void Coupling_Loss_CDF(void);
void Collect_LSP_from_ServingCell();
void LSP_DS_CDF();
void LSP_ASD_CDF();
void LSP_ASA_CDF();
void LSP_ZSD_CDF();
void LSP_ZSA_CDF();
void Geometry_CDF(void);
void Wideband_SIR_CDF(void);
void Delete_memory(void);
//void Delete_CHIR_memory(void);
void loading(void);
#ifdef ENABLE_PROGRESSBAR
void Progressbar(int);
void ResetProgressbar(int);
#endif
void Scheduling(void);
void scheduling_statistics();
void measure_statistics();
void ClearDrop();

void Print_Calib_Debug_Info();

// Singular Value CDF collection (multithreaded)
void Init_Singular_Value_Collection();
void Collect_Singular_Values_All();
void Singular_Value_CDF();

// Precoding-based metrics CDF collection (TX+RX digital beamforming)
void Init_Precoding_Metrics();
void Finalize_Precoding_Metrics();  // Compute time-averaged values before CDF
void Precoding_SINR_CDF();
void Precoding_Coupling_Loss_CDF();

int main(int argc, char *argv[])
{
    clock_t start, finish;
    Real duration;

	// Enable nested parallelism for OpenMP
	// This allows parallel regions inside parallel regions
	// Outer level: sector-level parallelism
	// Inner level: subband-level parallelism
	omp_set_nested(1);
	omp_set_max_active_levels(2);

	CH_CAL = 1;                  //// for runtime calculate
	Set_simul_param(argc, argv);
	Set_Parameter(scenario);
	StandardInitialization(); // Table Setting
	InitializeSystem(); 
	Set_fastfading_param();   // Correlatio matrix setting for LSP
	Allocate_memory();

	// Initialize singular value collection
	Init_Singular_Value_Collection();

	// Initialize precoding metrics collection
	Init_Precoding_Metrics();

	for (drop_idx = 0; drop_idx < num_drop; drop_idx++)
	{
		cout << "Drop # " << drop_idx+1 << "/" << num_drop << " ing... " << endl;
		cout << "  " << endl;
		
		Initialdrop();
		Initialize_CHANNEL();
		Generate_LSP();
		Generate_SSP();
		Print_Location(drop_idx);
		start = clock();
		Link_configuration();
		finish = clock();

		duration = (Real)(finish - start) / CLOCKS_PER_SEC;
    	cout << duration << "초" << endl;

		if (Calibration_mode == 0)  // normal simulation
		{
			ClearMAPs2StartADrop();
		}

		Get_CouplingLoss();
		Collect_LSP_from_ServingCell();
		Print_Calib_Debug_Info();

		int old_time_percent = -1;
		if (Calibration_mode == 0)  // normal simulation
		{
			Parameter_initialization();
			Compute_Channel_Coef();

#ifdef ENABLE_PROGRESSBAR
			ResetProgressbar( drop_idx );
#endif
			for (t = 0; t < run_times + SCHEDULE_DELAY; t++) // N_pf +
			{
#ifdef ENABLE_PROGRESSBAR
				Progressbar( t );
#endif
				Update_Channel_Coef();

				// Channel Update with progress indicator
				int processed_ms = 0;
				#if ENABLE_MULTITHREADING
				#pragma omp parallel num_threads(num_of_threads)
				{
					#pragma omp for
				#endif
					for (int idx = 0; idx < num_MS; idx++)
					{
						ms[idx].Channel_Update_MIMO(idx);
					}
				#if ENABLE_MULTITHREADING
				}
				#endif

				if (t > SCHEDULE_DELAY-1)  
				{
					Scheduling();
					scheduling_statistics();
				}
				

				if (t >= SCHEDULE_DELAY-1)
				{
					Measure();
					measure_statistics();
					// Note: Precoding metrics are collected inside Receive_DL_mTRP()
					// to avoid redundant SINR computation
				}
				// Collect singular values from channel matrices (multithreaded)
				Collect_Singular_Values_All();				
			}
			PerDropStatistics(drop_idx);
			ClearDrop();
			cout << "DROP finish.. " << endl;
		}
		else  
		{
			///// calibration mode -> only coupling loss, geometry
		}
	}
	Geometry_CDF();
	Wideband_SIR_CDF();
	Coupling_Loss_CDF();
	LSP_DS_CDF();
	LSP_ASD_CDF();
	LSP_ASA_CDF();
	LSP_ZSD_CDF();
	LSP_ZSA_CDF();
	Singular_Value_CDF();

	// Output precoding-based CDFs (TX+RX digital beamforming, time-averaged per UE)
	Finalize_Precoding_Metrics();  // Compute time-averaged values
	Precoding_SINR_CDF();
	Precoding_Coupling_Loss_CDF();

	return 0;
}




