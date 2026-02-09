#ifndef _BS_H_
#define _BS_H_

#include "const.h"

class BS
{
    public:
		BS();
		~BS();

	    LOCATION loc;

		LOCATION3D d_tx[3][32][64][2][2][2]; // sector(rx) [num_sector][M][N][P][Mg][Ng] - increased to max values (M=32, N=64)

		Real ant[3][3];  /// [sector num][0=alpha or 1=beta or 2=gamma]

};


#endif