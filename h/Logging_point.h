#include <complex>
#include <vector>
#ifdef ENABLE_PROGRESSBAR
#include "indicators.hpp"
#endif

using namespace std;

void VoIP_Statistics();
void Averaging_Sinr(int ue_idx, Real sinr);

// Singular Value CDF collection (3GPP compliant, multithreaded)
// 3GPP: PRB singular values = eigenvalues of mean covariance matrix
// Output: absolute SVs and ratio of (2nd, 3rd, ..., xth) to 1st (largest) SV in dB
// Collected at t=0 (initial time)
void Init_Singular_Value_Collection();
void Collect_Singular_Values_All();  // Multithreaded collection for all UEs
void Singular_Value_CDF();

extern vector<vector<Real>> singular_value_samples;       // Absolute SV values
extern vector<vector<Real>> singular_value_ratio_samples; // SV ratios in dB (3GPP: 10*log10(SV_n/SV_1))
extern int num_singular_values;

// Precoding-based metrics CDF collection (TX+RX digital beamforming)
// Time-averaged per UE: accumulate over all time slots, then average for CDF
void Init_Precoding_Metrics();
void Accumulate_Precoding_Metrics(int ue_idx, Real wideband_sinr_linear, Real coupling_loss_dB);  // Called from Receive_DL
void Finalize_Precoding_Metrics();     // Compute time-averaged values before CDF output
void Precoding_SINR_CDF();
void Precoding_Coupling_Loss_CDF();

extern vector<Real> precoding_wideband_sinr_samples;    // Precoding-applied wideband SINR (time-averaged per UE)
extern vector<Real> precoding_coupling_loss_samples;    // Precoding-applied coupling loss (time-averaged per UE)

int num_ues_in_outage;

Real * throuhput_ue;

Real     * avr_sinr;
int        * num_of_added;

Real     * geometry;
Real     * wideband_sir;
Real     * coupling_loss;
Real     * RMS_delay_spread;
Real     * circular_angle_spread_AOA;
Real     * circular_angle_spread_AOD;

Real     * LSP_DS;
Real     * LSP_ASD;
Real     * LSP_ASA;
Real     * LSP_ZSD;
Real     * LSP_ZSA;

void Collect_LSP_from_ServingCell();
void Get_LSP_DS(Real val, int ue_idx);
void Get_LSP_ASD(Real val, int ue_idx);
void Get_LSP_ASA(Real val, int ue_idx);
void Get_LSP_ZSD(Real val, int ue_idx);
void Get_LSP_ZSA(Real val, int ue_idx);
void LSP_DS_CDF();
void LSP_ASD_CDF();
void LSP_ASA_CDF();
void LSP_ZSD_CDF();
void LSP_ZSA_CDF();

#ifdef ENABLE_PROGRESSBAR
using namespace indicators;
ProgressBar * bar;
#endif
/*
ProgressBar bar{
	option::BarWidth{50},
	option::Start{"["},
	option::Fill{"■"},
	option::Lead{"■"},
	option::Remainder{"-"},
	option::End{" ]"},
	option::PrefixText{"Drop"},
	option::ShowElapsedTime{true},
	option::ShowRemainingTime{true},    
	option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}
};
*/
#ifdef ENABLE_PROGRESSBAR
int old_time_percent = -1;
#endif


/*===================================================================
FUNCTION: STASTICS()

AUTHOR: Noh, jung-hoon

DESCRIPTION: Structure for getting cdf , pdf histograms
PDF will be added

NOTES: 2009-08-12

===================================================================*/
struct STASTICS
{
	Real * samples;
	int         size;
	char * name;

	vector < Real > cdf;
	// setting up basic parameter 

	void setup(Real * _samples, int _size, char * _name)
	{
		samples = _samples;
		size = _size;
		name = _name;
		cdf.clear();
	}
	
	// get cdf ( percetile type 0 ~100% )
	void get_cdf(void)
	{
		
		sort(samples, samples + size); 
		
		cdf.push_back(samples[0]);

		for (int percentile = 1; percentile < 101; percentile++)
		{
			//int idx = (int)itpp::round(size*percentile*0.01);
			int idx = (int)((size-1)*percentile*0.01);
			cdf.push_back(samples[idx]);
		}
		
	}

	// Print Out CDF (with date/time in filename)
	void print_cdf(void)
	{
		ofstream print_out_cdf;
		char cdf_path[200];

		// Include date/time in filename: {name}_{YYYY-MM-DD_HH-MM}_cdf.dat
		sprintf(cdf_path, "./%s/%s_%04d-%02d-%02d_%02d-%02d_cdf.dat",
				folder_name, name, g_year, g_mon, g_day, g_hour, g_min);

		print_out_cdf.open(cdf_path, ios::out);

		for (int percentile = 0; percentile < 101; percentile++)
		{
			print_out_cdf << percentile << "\t" << cdf[percentile] << endl;
		}
		print_out_cdf.close();
	}
};






