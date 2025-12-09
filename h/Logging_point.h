#include <complex>
#include <vector>
#include "indicators.hpp"

using namespace std;

void VoIP_Statistics();
void Averaging_Sinr(int ue_idx, Real sinr);


int num_ues_in_outage;

Real * throuhput_ue;

Real     * avr_sinr;
int        * num_of_added;

Real     * geometry;
Real     * coupling_loss;
Real     * RMS_delay_spread;
Real     * circular_angle_spread_AOA;
Real     * circular_angle_spread_AOD;

using namespace indicators;
ProgressBar * bar;
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
int old_time_percent = -1;


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

	// Print Out CDF
	void print_cdf(void)
	{
		/*
		ofstream print_out_cdf;
		char cdf_path[50];
		strcpy(cdf_path, path);
		strcat(cdf_path, "/");
		strcat(cdf_path, name);
		strcat(cdf_path, "_cdf.dat");

		print_out_cdf.open(cdf_path, ios::out);

		for (int percentile = 0; percentile < 101; percentile++)
		{
			print_out_cdf << percentile << "\t" << cdf[percentile] << endl;
		}
		print_out_cdf.close();

		*/

		ofstream print_out_cdf;
		char cdf_path[100];
		sprintf(cdf_path,"./%s/%s_cdf.dat",folder_name,name);

		//strcpy(cdf_path, name);
		//strcat(cdf_path, "_cdf.dat");

		print_out_cdf.open(cdf_path, ios::out);

		for (int percentile = 0; percentile < 101; percentile++)
		{
			print_out_cdf << percentile << "\t" << cdf[percentile] << endl;
			
			//cout << cdf[percentile] << endl;
		}
		print_out_cdf.close();
	}
};






