#include "common.h"


void Print_Location(int drop_idx)
{
	if (drop_idx == (num_drop - 1))
	{
		ofstream Location("BS_MS_Location.m");
		Location << "clc;clear all;" << endl;

		Location << "BS= [";
		for (int bs_idx = 0; bs_idx < num_BS + num_mTRP; bs_idx++)
		{
			Location << bs[bs_idx].loc.x << "," << bs[bs_idx].loc.y << "," << bs_height  << ";\t";
		}
		Location << "];" << endl;

		Location << "MS= [";
		for (int ms_idx = 0; ms_idx < num_MS; ms_idx++)
		{
			Location << ms[ms_idx].loc.x << "," << ms[ms_idx].loc.y << "," << ms[ms_idx].MS_HEIGHT_FINAL << ";\t";
		}
		Location << "];" << endl;

		Location << "figure()" << endl;
		Location << "stem3(BS(:,1),BS(:,2),BS(:,3),'ro', 'filled');hold on;grid on;grid minor;" << endl;
		Location << "stem3(MS(:,1),MS(:,2),MS(:,3),'x', 'filled');" << endl;
		Location << "title('BS and MS Location');" << endl;
		Location << "legend('BS','MS');" << endl;
		Location.close();
	}
	

}

Real Link_SINR[200];
Real Link_SINR_PDF[200];
int total_interest_links;


void CDF(int drop_idx)
{
	Real SINR_sample;

	if (drop_idx == 0)
	{
		total_interest_links = 0;
		for (int idx = 0; idx < 200; idx++)
		{
			Link_SINR[idx] = 0;
		}
	}

	for (int link_idx = 0; link_idx < num_LINK; link_idx++)
	{

		SINR_sample = links[link_idx].SINR;
		total_interest_links++;

		for (int idx = 0; idx < 200; idx++)
		{
			if (SINR_sample <= -15. + idx*0.25)
			{
				Link_SINR[idx] ++;
			}
		}
	}


	if (drop_idx == (num_drop - 1))
	{
		
		ofstream  SINR;
		ofstream SINR_excel;
		SINR.open("SINR_DL.m", ios::out);
		SINR_excel.open("SINR_DL.xls", ios::out);
		SINR << "clc;clear all;" << endl;
		SINR << "sinr= [";
		for (int idx = 0; idx < 200; idx++)
		{
			Link_SINR_PDF[idx] = Link_SINR[idx] / (Real)(total_interest_links);
			SINR << -15. + idx*0.25 << "," << Link_SINR_PDF[idx] << ";\t";
			SINR_excel << -15. + idx*0.25 << "\t" << Link_SINR_PDF[idx] << endl;
		}
		SINR << "];" << endl;
		SINR << "figure()" << endl;
		SINR << "plot(sinr(:,1),sinr(:,2));" << endl;
		SINR << "xlabel('SINR(dB)');ylabel('CDF');"<<endl;
		SINR.close();
		SINR_excel.close();
	}
}