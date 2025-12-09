
#include "common.h"
//#include "MS.h"

void Allocate_memory(void)
{
	for (int ue_idx = 0; ue_idx < num_MS; ue_idx++)
	{
		ms[ue_idx].Allocate_memory();
	}

	for (int bs_idx = 0; bs_idx < num_BS + num_mTRP; bs_idx++)
	{
		for (int ue_idx = 0; ue_idx < num_MS; ue_idx++)
		{
			channel[bs_idx][ue_idx].Allocate_memory();
		}
	}

	TX_LOS_gain_theta = new ComplexReal[NUM_TX_Port];
	TX_LOS_gain_pi    = new ComplexReal[NUM_TX_Port];

	TX_NLOS_gain_theta = new ComplexReal **[NUM_TX_Port];
	TX_NLOS_gain_pi    = new ComplexReal **[NUM_TX_Port];

	for (int p = 0; p < NUM_TX_Port; p++)
	{
		TX_NLOS_gain_theta[p] = new ComplexReal *[MAX_NUM_CLUSTERS];
		TX_NLOS_gain_pi   [p] = new ComplexReal *[MAX_NUM_CLUSTERS];
		 
		for (int i = 0; i < MAX_NUM_CLUSTERS; i++)
		{
			TX_NLOS_gain_theta[p][i] = new ComplexReal[MAX_NUM_RAYS];
			TX_NLOS_gain_pi   [p][i] = new ComplexReal[MAX_NUM_RAYS];
		}
	}

	RX_LOS_gain_theta = new ComplexReal[NUM_RX_Port];
	RX_LOS_gain_pi    = new ComplexReal[NUM_RX_Port];

	RX_NLOS_gain_theta = new ComplexReal **[NUM_RX_Port];
	RX_NLOS_gain_pi    = new ComplexReal **[NUM_RX_Port];

	for (int p = 0; p < NUM_RX_Port; p++)
	{
		RX_NLOS_gain_theta[p] = new ComplexReal *[MAX_NUM_CLUSTERS];
		RX_NLOS_gain_pi[p]    = new ComplexReal *[MAX_NUM_CLUSTERS];

		for (int i = 0; i < MAX_NUM_CLUSTERS; i++)
		{
			RX_NLOS_gain_theta[p][i] = new ComplexReal[MAX_NUM_RAYS];
			RX_NLOS_gain_pi[p][i]    = new ComplexReal[MAX_NUM_RAYS];
		}
	}

	comp_mode = new int [num_MS];

}

void Delete_memory(void)
{
	// CHIR 메모리 해제 (할당된 것만)
	for (int bs_idx = 0; bs_idx < num_BS + num_mTRP; bs_idx++)
	{
		for (int ue_idx = 0; ue_idx < num_MS; ue_idx++)
		{
			if (channel[bs_idx][ue_idx].CHIR_allocated)
			{
				channel[bs_idx][ue_idx].Delete_CHIR_memory();
			}
			channel[bs_idx][ue_idx].Delete_memory();
		}
		delete[] channel[bs_idx];
	}
	delete[] channel;

	// MS 메모리 해제
	for (int ue_idx = 0; ue_idx < num_MS; ue_idx++)
	{
		ms[ue_idx].Delete_memory();
	}

	// TX/RX antenna gain 메모리 해제
	delete[] TX_LOS_gain_theta;
	delete[] TX_LOS_gain_pi;

	for (int p = 0; p < NUM_TX_Port; p++)
	{
		for (int i = 0; i < MAX_NUM_CLUSTERS; i++)
		{
			delete[] TX_NLOS_gain_theta[p][i];
			delete[] TX_NLOS_gain_pi[p][i];
		}
		delete[] TX_NLOS_gain_theta[p];
		delete[] TX_NLOS_gain_pi[p];
	}
	delete[] TX_NLOS_gain_theta;
	delete[] TX_NLOS_gain_pi;

	delete[] RX_LOS_gain_theta;
	delete[] RX_LOS_gain_pi;

	for (int p = 0; p < NUM_RX_Port; p++)
	{
		for (int i = 0; i < MAX_NUM_CLUSTERS; i++)
		{
			delete[] RX_NLOS_gain_theta[p][i];
			delete[] RX_NLOS_gain_pi[p][i];
		}
		delete[] RX_NLOS_gain_theta[p];
		delete[] RX_NLOS_gain_pi[p];
	}
	delete[] RX_NLOS_gain_theta;
	delete[] RX_NLOS_gain_pi;

	delete[] comp_mode;
}