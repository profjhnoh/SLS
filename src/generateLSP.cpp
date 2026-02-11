#include "common.h"
#include "const.h"
#include "rv_generate.h"
#include "paramFastFading.h"
#include <iostream>
#include <Eigen/Dense>
#include <vector> 

Real normal(Real mu, Real sigma);

void new_memory();
void Set_grid_size();
void Generate_TLSPs_RN_map();

void delete_memory();
void Set_auto_corr_matrix();
void Get_TLSPs_by_interpolation();
void Transform_TLSP_to_LSP();
void Check_corr_var();

Real h_m(Real);

VectorXReal *****TLSPs_map; // [BS][floor][condition][x][y] (LSP)
VectorXReal **corr_RN_vecter; //[UE][BS] (LSP)

MatrixXReal **R; // [condition][LSP]
MatrixXReal **inv_L_tilde;  // [condition][LSP]
RowVectorXReal **lamda_n_T;  // [condition][LSP]

int size_x_grid, size_y_grid;
LOCATION grid_origin_point;                  
Real max_x, min_x, max_y, min_y;

void Generate_LSP()
{
	Set_grid_size();
	new_memory();
	Set_auto_corr_matrix();
	cout << "Generate_TLSPs_RN_map............."<<std::flush;
	Generate_TLSPs_RN_map();
	cout << "DONE" << endl;
	cout << "Get_TLSPs_by_interpolation........"<<std::flush; 
	Get_TLSPs_by_interpolation();
	cout <<"DONE" << endl;

	Transform_TLSP_to_LSP();

	//Check_corr_var();
	delete_memory();
}

void Set_grid_size()
{
	if (test_environment<3)
	{
		size_x_grid = int(120. / grid_interval + 1) + 2;
		size_y_grid = int(50 / grid_interval + 1) + 2;
	}
	else
	{
		size_x_grid = int(inter_site_distance * 8 / sqrt(3.) / grid_interval + 1) + 2;
		size_y_grid = int(inter_site_distance * 5 / grid_interval + 1) + 2;
	}


	//// grid origin point : the highest, leftmost grid point
	if (test_environment<3)
	{
		grid_origin_point.x = -1 * (((size_x_grid - 1) * grid_interval) / 2);
		grid_origin_point.y = (((size_y_grid - 1) * grid_interval) / 2) + 25;

	}
	else
	{
		grid_origin_point.x = -1 * (((size_x_grid - 1) * grid_interval) / 2);
		grid_origin_point.y = ((size_y_grid - 1) * grid_interval) / 2;
	}

	
}


void new_memory()
{
	cout << "TLSPs_map........................."<<std::flush;
	TLSPs_map = new VectorXReal ****[num_BS + num_mTRP ];

	for (int bs_idx = 0; bs_idx < num_BS + num_mTRP; bs_idx++)
	{
		TLSPs_map[bs_idx] = new VectorXReal ***[num_floor];

		for (int floor_idx = 0; floor_idx < num_floor; floor_idx++)
		{
			TLSPs_map[bs_idx][floor_idx] = new VectorXReal**[num_propagation_condition];

			for (int idx = 0; idx < num_propagation_condition; idx++)
			{
				TLSPs_map[bs_idx][floor_idx][idx] = new VectorXReal*[size_x_grid];

				for (int x = 0; x < size_x_grid; x++)
				{
					TLSPs_map[bs_idx][floor_idx][idx][x] = new VectorXReal[size_y_grid];

					for (int y = 0; y < size_y_grid; y++)
					{
						if (idx == LOS_propagation)
						{
							TLSPs_map[bs_idx][floor_idx][idx][x][y] = VectorXReal(7);
						}
						else
						{
							TLSPs_map[bs_idx][floor_idx][idx][x][y] = VectorXReal(6);
						}
					}
				}
			}
		}
	}
	cout<<"DONE"<< endl;
    cout <<"corr_RN_vector...................."<<std::flush;
	corr_RN_vecter = new VectorXReal *[num_MS];

	for (int ue_idx = 0; ue_idx < num_MS; ue_idx++)
	{
		corr_RN_vecter[ue_idx] = new VectorXReal[num_BS + num_mTRP];

		for (int bs_idx = 0; bs_idx < num_BS + num_mTRP; bs_idx++)
		{
			if (channel[bs_idx][ue_idx].Propagation == LOS_propagation)  ///// LOS propagation = LOS = 1
			{
				corr_RN_vecter[ue_idx][bs_idx] = VectorXReal(7);
			}
			else
			{
				corr_RN_vecter[ue_idx][bs_idx] = VectorXReal(6);
			}
		}
	}
	cout << "DONE" << endl;
	cout <<"TLSPs............................."<<std::flush;

	for (int ue_idx = 0; ue_idx < num_MS; ue_idx++)
	{
		for (int bs_idx = 0; bs_idx < num_BS + num_mTRP; bs_idx++)
		{
			if (channel[bs_idx][ue_idx].Propagation == LOS_propagation)    ///// LOS propagation = LOS = 1
			{
				ms[ue_idx].TLSPs[bs_idx] = ArrayXReal(7);
				ms[ue_idx].LSPs[bs_idx] = ArrayXReal(7);
			}
			else
			{
				ms[ue_idx].TLSPs[bs_idx] = ArrayXReal(6);
				ms[ue_idx].LSPs[bs_idx] = ArrayXReal(6);
			}

		}
	}
	cout << "DONE" << endl;

	R = new MatrixXReal*[num_propagation_condition];
	inv_L_tilde = new MatrixXReal*[num_propagation_condition];
	lamda_n_T = new RowVectorXReal*[num_propagation_condition];

	for (int idx = 0; idx < num_propagation_condition; idx++)
	{
		int num_LSP;

		if (idx == LOS_propagation) { num_LSP = 7; }
		else                        { num_LSP = 6; }

		R[idx] = new MatrixXReal[num_LSP];
		inv_L_tilde[idx] = new MatrixXReal[num_LSP];
		lamda_n_T[idx] = new RowVectorXReal[num_LSP];

		for (int LSP_idx = 0; LSP_idx < num_LSP; LSP_idx++)
		{
			R[idx][LSP_idx] = MatrixXReal(num_neighbor + 1, num_neighbor + 1);
			inv_L_tilde[idx][LSP_idx] = MatrixXReal(num_neighbor, num_neighbor);
		}

	}

}

void delete_memory()
{
	for (int bs_idx = 0; bs_idx < num_BS + num_mTRP; bs_idx++)
	{
		for (int floor_idx = 0; floor_idx < num_floor; floor_idx++)
		{
			for (int idx = 0; idx < num_propagation_condition; idx++)
			{
				for (int x = 0; x < size_x_grid; x++)
				{
					delete[]TLSPs_map[bs_idx][floor_idx][idx][x];
				}

				delete[]TLSPs_map[bs_idx][floor_idx][idx];

			}
			delete[]TLSPs_map[bs_idx][floor_idx];
		}
		delete[]TLSPs_map[bs_idx];
	}
	delete[]TLSPs_map;


	for (int ue_idx = 0; ue_idx < num_MS; ue_idx++)
	{
		delete[]corr_RN_vecter[ue_idx];
	}

	delete[]corr_RN_vecter;


	for (int idx = 0; idx < num_propagation_condition; idx++)
	{
		delete[]R[idx];
		delete[]lamda_n_T[idx];
		delete[]inv_L_tilde[idx];
	}

	delete[]R;
	delete[]lamda_n_T;
	delete[]inv_L_tilde;
}


void Set_auto_corr_matrix()
{
	Real alpha[7];
	int num_LSP;

	MatrixXReal L = MatrixXReal(num_neighbor + 1, num_neighbor + 1);
	MatrixXReal L_tilde = MatrixXReal(num_neighbor, num_neighbor);

	for (int propagation_condition = 0; propagation_condition < num_propagation_condition; propagation_condition++)
	{
		if (propagation_condition == LOS_propagation)	{ num_LSP = 7; }
		else  		                                    { num_LSP = 6; }

		if (propagation_condition == LOS_propagation)  // LOS case 7 Large-scale parameters else 6 LSP
		{
			for (int LSP_idx = 0; LSP_idx < num_LSP; LSP_idx++)
			{
				if (LSP_idx == 0) { alpha[0] = corr_dist_DS[propagation_condition]; }
				else if (LSP_idx == 1) { alpha[1] = corr_dist_ASD[propagation_condition]; }
				else if (LSP_idx == 2) { alpha[2] = corr_dist_ASA[propagation_condition]; }
				else if (LSP_idx == 3) { alpha[3] = corr_dist_SF[propagation_condition]; }
				else if (LSP_idx == 4) { alpha[4] = corr_dist_K[propagation_condition]; }
				else if (LSP_idx == 5) { alpha[5] = corr_dist_ZSD[propagation_condition]; }
				else if (LSP_idx == 6) { alpha[6] = corr_dist_ZSA[propagation_condition]; }
			}
		}
		else
		{
			for (int LSP_idx = 0; LSP_idx < num_LSP; LSP_idx++)
			{
				if (LSP_idx == 0) { alpha[0] = corr_dist_DS[propagation_condition]; }
				else if (LSP_idx == 1) { alpha[1] = corr_dist_ASD[propagation_condition]; }
				else if (LSP_idx == 2) { alpha[2] = corr_dist_ASA[propagation_condition]; }
				else if (LSP_idx == 3) { alpha[3] = corr_dist_SF[propagation_condition]; }
				else if (LSP_idx == 4) { alpha[4] = corr_dist_ZSD[propagation_condition]; }
				else if (LSP_idx == 5) { alpha[5] = corr_dist_ZSA[propagation_condition]; }
			}
		}


		for (int LSP_idx = 0; LSP_idx < num_LSP; LSP_idx++)
		{
			Real d = grid_interval / 20;

			if (num_neighbor == 4)
			{
				R[propagation_condition][LSP_idx] <<	1,				h_m(d),			h_m(2 * d),		h_m(d),			h_m(sqrt(2)*d),
														h_m(d),			1,				h_m(d),			h_m(sqrt(2)*d), h_m(d),
														h_m(2 * d),		h_m(d),			1,				h_m(sqrt(5)*d), h_m(sqrt(2)*d),
														h_m(d),			h_m(sqrt(2)*d), h_m(sqrt(5)*d), 1,				h_m(d),
														h_m(sqrt(2)*d), h_m(d),			h_m(sqrt(2)*d), h_m(d),			1;
			}
			else if (num_neighbor == 8)
			{
				R[propagation_condition][LSP_idx] <<	1,				h_m(d),			h_m(2 * d),		h_m(d),			 h_m(d),		  h_m(sqrt(5)*d),  h_m(d),			h_m(3 * d),		 h_m(sqrt(2)*d),
														h_m(d),			1,				h_m(d),			h_m(sqrt(2)*d),  h_m(sqrt(2)*d),  h_m(sqrt(2)*d),  h_m(2 * d),		h_m(2 * d),		 h_m(d),
														h_m(2 * d),		h_m(d),			1,				h_m(sqrt(5)*d),  h_m(sqrt(5)*d),  h_m(d),		   h_m(3 * d),		h_m(d),			 h_m(sqrt(2)*d),
														h_m(d),			h_m(sqrt(2)*d), h_m(sqrt(5)*d), 1,				 h_m(2 * d),	  h_m(sqrt(8)*d),  h_m(sqrt(2)*d),  h_m(sqrt(10)*d), h_m(d),
														h_m(d),			h_m(sqrt(2)*d), h_m(sqrt(5)*d), h_m(2 * d),		 1,				  h_m(2 * d),	   h_m(sqrt(2)*d),  h_m(sqrt(10)*d), h_m(sqrt(5)*d),
														h_m(sqrt(5)*d), h_m(sqrt(2)*d), h_m(d),			h_m(sqrt(8)*d),  h_m(2 * d),	  1,			   h_m(sqrt(10)*d), h_m(sqrt(2)*d),  h_m(sqrt(5)*d),
														h_m(d),			h_m(2 * d),		h_m(3 * d),		h_m(sqrt(2)*d),  h_m(sqrt(2)*d),  h_m(sqrt(10)*d), 1,				h_m(4 * d),		 h_m(sqrt(5)*d),
														h_m(3 * d),		h_m(2 * d),		h_m(d),			h_m(sqrt(10)*d), h_m(sqrt(10)*d), h_m(sqrt(2)*d),  h_m(4 * d),		1,				 h_m(sqrt(5)*d),
														h_m(sqrt(2)*d), h_m(d),			h_m(sqrt(2)*d), h_m(d),			 h_m(sqrt(5)*d),  h_m(sqrt(5)*d),  h_m(sqrt(5)*d),  h_m(sqrt(5)*d),  1;
			}
			else
			{
				cout << "check num_neighbor = 4 or 8" << endl;
			}

		}

		for (int LSP_idx = 0; LSP_idx < num_LSP; LSP_idx++)
		{
			LLT<MatrixXReal> LLT_R(R[propagation_condition][LSP_idx]);
			L = LLT_R.matrixL();
			L_tilde = L.block(0, 0, num_neighbor, num_neighbor); //   L_tilde = LLT_R_tilde= L(size_row-1, size_col-1)

			inv_L_tilde[propagation_condition][LSP_idx] = L_tilde.inverse();
			lamda_n_T[propagation_condition][LSP_idx] = L.row(num_neighbor);
		}
	}

}


void Generate_TLSPs_RN_map()
{
	int num_LSP;
	MatrixXReal iid_RN_map(size_x_grid, size_y_grid); // gaussian iid Random number grid

	VectorXReal s_tilde(num_neighbor);
	VectorXReal a(num_neighbor + 1);
	LOCATION *relative_position_to_s_n;
	relative_position_to_s_n = new LOCATION[num_neighbor];

	vector <int> **available_neighbor;
	available_neighbor = new vector <int>*[size_x_grid];


	for (int x = 0; x<size_x_grid; x++)
	{
		available_neighbor[x] = new vector <int>[size_y_grid];
	}

	if (num_neighbor == 4)
	{
		// s3 sn
		// s0 s1 s2
		relative_position_to_s_n[0].x = -1;
		relative_position_to_s_n[0].y = -1;
		relative_position_to_s_n[1].x = 0;
		relative_position_to_s_n[1].y = -1;
		relative_position_to_s_n[2].x = 1;
		relative_position_to_s_n[2].y = -1;
		relative_position_to_s_n[3].x = -1;
		relative_position_to_s_n[3].y = 0;

	}
	else if (num_neighbor == 8)
	{
		//
		//   s3 sn
		//s6 s0 s1 s2 s7
		//   s4    s5
		relative_position_to_s_n[0].x = -1;
		relative_position_to_s_n[0].y = -1;
		relative_position_to_s_n[1].x = 0;
		relative_position_to_s_n[1].y = -1;
		relative_position_to_s_n[2].x = 1;
		relative_position_to_s_n[2].y = -1;
		relative_position_to_s_n[3].x = -1;
		relative_position_to_s_n[3].y = 0;
		relative_position_to_s_n[4].x = -1;
		relative_position_to_s_n[4].y = -2;
		relative_position_to_s_n[5].x = 1;
		relative_position_to_s_n[5].y = -2;
		relative_position_to_s_n[6].x = -2;
		relative_position_to_s_n[6].y = -1;
		relative_position_to_s_n[7].x = 2;
		relative_position_to_s_n[7].y = -1;
	}
	else
	{
		cout << "check num_neighbor = 4 or 8" << endl;
	}

	for (int y = 0; y < size_y_grid; y++)
	{
		for (int x = 0; x < size_x_grid; x++)
		{
			for (int neighbor_idx = 0; neighbor_idx < num_neighbor; neighbor_idx++)
			{
				LOCATION neighbor_location;

				neighbor_location.x = x + relative_position_to_s_n[neighbor_idx].x;
				neighbor_location.y = y + relative_position_to_s_n[neighbor_idx].y;
				if (0 <= neighbor_location.x  && neighbor_location.x < size_x_grid && 0 <= neighbor_location.y  && neighbor_location.y < size_y_grid)
				{
					available_neighbor[x][y].push_back(neighbor_idx);
				}
			}
			available_neighbor[x][y].push_back(num_neighbor); // self position
		}
	}

	for (int bs_idx = 0; bs_idx < num_BS + num_mTRP; bs_idx++)
	{
		//cout << "bs_idx : " << bs_idx + 1 << "/" << num_BS << endl;

		for (int floor_idx = 0; floor_idx < num_floor; floor_idx++)
		{
			for (int propagation_condition = 0; propagation_condition < num_propagation_condition; propagation_condition++)
			{
				if (propagation_condition == LOS_propagation)	{ num_LSP = 7; }
				else 								            { num_LSP = 6; }

				VectorXReal **corr_RN_map;

				corr_RN_map = new VectorXReal *[size_x_grid];

				for (int x = 0; x< size_x_grid; x++)
				{
					corr_RN_map[x] = new VectorXReal[size_y_grid];
					for (int y = 0; y< size_y_grid; y++)
					{
						corr_RN_map[x][y] = VectorXReal(num_LSP);
					}
				}


				for (int LSP_idx = 0; LSP_idx < num_LSP; LSP_idx++)
				{

					RowVectorXReal _lamda_n_T(num_neighbor + 1);
					Real a_n;
					Real s_n;
					_lamda_n_T << lamda_n_T[propagation_condition][LSP_idx];


					for (int x = 0; x < size_x_grid; x++)
					{
						for (int y = 0; y < size_y_grid; y++)
						{
							iid_RN_map(x, y) = randnum.n();
						}
					}


					for (int y = 0; y < size_y_grid; y++)
					{

						for (int x = 0; x < size_x_grid; x++)
						{
							a_n = iid_RN_map(x, y);

							if (available_neighbor[x][y].size() == num_neighbor + 1)
							{
								for (int s_idx = 0; s_idx < num_neighbor; s_idx++)
								{
									int neighbor_pos_x = x + int(relative_position_to_s_n[s_idx].x);
									int neighbor_pos_y = y + int(relative_position_to_s_n[s_idx].y);

									s_tilde(s_idx) = corr_RN_map[neighbor_pos_x][neighbor_pos_y](LSP_idx);
								}

								a << inv_L_tilde[propagation_condition][LSP_idx] * s_tilde, a_n;
								s_n = _lamda_n_T * a;							
							}
							else
							{
								int num_available_neighbor = available_neighbor[x][y].size() - 1;

								MatrixXReal _R(num_available_neighbor + 1, num_available_neighbor + 1);
								MatrixXReal _L(num_available_neighbor + 1, num_available_neighbor + 1);
								MatrixXReal _L_tilde(num_available_neighbor, num_available_neighbor);
								MatrixXReal _inv_L_tilde(num_available_neighbor, num_available_neighbor);
								RowVectorXReal _lamda_n_T(num_available_neighbor + 1);
								VectorXReal _a(num_available_neighbor + 1);
								VectorXReal _s_tilde(num_available_neighbor);


								if (num_available_neighbor > 0)
								{

									for (int available_neighbor_idx = 0; available_neighbor_idx < num_available_neighbor; available_neighbor_idx++)
									{
										int neighbor_pos_x = x + int(relative_position_to_s_n[available_neighbor[x][y][available_neighbor_idx]].x);
										int neighbor_pos_y = y + int(relative_position_to_s_n[available_neighbor[x][y][available_neighbor_idx]].y);

										_s_tilde(available_neighbor_idx) = corr_RN_map[neighbor_pos_x][neighbor_pos_y](LSP_idx);

									}
								}

								for (int R_prime_row_idx = 0; R_prime_row_idx < num_available_neighbor + 1; R_prime_row_idx++)
								{
									for (int R_prime_column_idx = 0; R_prime_column_idx < num_available_neighbor + 1; R_prime_column_idx++)
									{
										_R(R_prime_row_idx, R_prime_column_idx) = R[propagation_condition][LSP_idx](available_neighbor[x][y][R_prime_row_idx], available_neighbor[x][y][R_prime_column_idx]);
									}

								}

								if (num_available_neighbor > 0)
								{
									LLT<MatrixXReal> LLT_R(_R);
									_L = LLT_R.matrixL();
									_L_tilde = _L.block(0, 0, num_available_neighbor, num_available_neighbor); //   L_tilde = LLT_R_tilde= L(size_row-1, size_col-1)
									_inv_L_tilde = _L_tilde.inverse();
									_lamda_n_T = _L.row(num_available_neighbor);
									_a << _inv_L_tilde * _s_tilde, a_n;
									s_n = _lamda_n_T * _a;
								}
								else
								{
									s_n = a_n;
								}
							}

							corr_RN_map[x][y](LSP_idx) = s_n;
						}
					}

				}


				for (int x = 0; x< size_x_grid; x++)
				{
					for (int y = 0; y< size_y_grid; y++)
					{
						switch (propagation_condition)
						{
						case LOS_propagation:
							TLSPs_map[bs_idx][floor_idx][propagation_condition][x][y] = sqrt_corr_matrx_LOS * corr_RN_map[x][y];
							break;
						case NLOS_propagation:
							TLSPs_map[bs_idx][floor_idx][propagation_condition][x][y] = sqrt_corr_matrx_NLOS* corr_RN_map[x][y];
							break;
						case OUT2IN_propagation:
							TLSPs_map[bs_idx][floor_idx][propagation_condition][x][y] = sqrt_corr_matrx_OUT2IN* corr_RN_map[x][y];
							break;
						default:
							cout << "check propagation_condition " << endl;
						}
					}
				}
				delete[]corr_RN_map;
			}
		}
	}


	for (int x = 0; x<size_x_grid; x++)
	{
		delete[]available_neighbor[x];
	}
	delete[] available_neighbor;
	delete[] relative_position_to_s_n;

}
/*===================================================================
FUNCTION: Get_TLSPs_by_interpolation()

AUTHOR: G.J

DESCRIPTION:
-

NOTES:

===================================================================*/

void Get_TLSPs_by_interpolation()
{
	int x_[4];
	int y_[4];
	Real x_pos;
	Real y_pos;
	Real D_corr = grid_interval;

	VectorXReal (6);
	VectorXReal radom_NLOS(7);
	for (int bs_idx = 0; bs_idx < num_BS + num_mTRP ; bs_idx++)
	{

		for (int ue_idx = 0; ue_idx < num_MS; ue_idx++)
		{
			ArrayXReal random_V;
			int propagation_condition = channel[bs_idx][ue_idx].Propagation;
			if (propagation_condition == LOS_propagation)
			{
				random_V = ArrayXReal(7);
			
				for (int i = 0; i < 7; i++)
				{
					random_V(i) = randnum.n();
				}
			}
			else
			{
				random_V = ArrayXReal(6);
				for (int i = 0; i < 6; i++)
				{
					random_V(i) = randnum.n();
				}
			
			}

			//cout << "bs_idx = "<< bs_idx << endl;
			switch (propagation_condition)
			{
				case LOS_propagation:
					ms[ue_idx].TLSPs[bs_idx] = sqrt_corr_matrx_LOS * random_V.matrix();
					break;
				case NLOS_propagation:
					ms[ue_idx].TLSPs[bs_idx] = sqrt_corr_matrx_NLOS* random_V.matrix();
					break;
				case OUT2IN_propagation:
					ms[ue_idx].TLSPs[bs_idx] = sqrt_corr_matrx_OUT2IN* random_V.matrix();
					break;
				default:
					cout << "check propagation_condition " << endl;
			}
		}
	}
}


Real h_m(Real x)
{
	return exp(-x);
}

void Check_corr_var()
{
	MatrixXReal corr_var_LOS = MatrixXReal::Zero(7, 7);
	MatrixXReal corr_var_NLOS = MatrixXReal::Zero(6, 6);
	MatrixXReal corr_var_OUT2IN = MatrixXReal::Zero(6, 6);

	int num_LSP;

	for (int bs_idx = 0; bs_idx < num_BS + num_mTRP; bs_idx++)
	{
		for (int floor_idx = 0; floor_idx < num_floor; floor_idx++)
		{
			for (int propagation_condition = 0; propagation_condition < num_propagation_condition; propagation_condition++)
			{
				if (propagation_condition == LOS_propagation) { num_LSP = 7; }
				else							  { num_LSP = 6; }
				for (int x = 0; x < size_x_grid; x++)
				{
					for (int y = 0; y < size_y_grid; y++)
					{
						switch (propagation_condition)
						{
						case LOS_propagation:
							corr_var_LOS += TLSPs_map[bs_idx][floor_idx][propagation_condition][x][y].matrix()*TLSPs_map[bs_idx][floor_idx][propagation_condition][x][y].matrix().transpose();
							break;
						case NLOS_propagation:
							corr_var_NLOS += TLSPs_map[bs_idx][floor_idx][propagation_condition][x][y].matrix()*TLSPs_map[bs_idx][floor_idx][propagation_condition][x][y].matrix().transpose();
							break;
						case OUT2IN_propagation:
							corr_var_OUT2IN += TLSPs_map[bs_idx][floor_idx][propagation_condition][x][y].matrix()*TLSPs_map[bs_idx][floor_idx][propagation_condition][x][y].matrix().transpose();
							break;
						default:
							cout << "check propagation_condition " << endl;
						}
					}
				}
			}
		}
	}

	cout << "--------------------cross correlatin at grid --------------------------" << endl;
	cout << "var_LOS  " << endl;
	cout << corr_var_LOS / num_BS / num_floor / size_x_grid / size_y_grid << endl;
	cout << endl;

	cout << "var_NLOS  " << endl;
	cout << corr_var_NLOS / num_BS / num_floor / size_x_grid / size_y_grid << endl;

	cout << endl;


	cout << "var_OUT2IN " << endl;
	cout << corr_var_OUT2IN / num_BS / num_floor / size_x_grid / size_y_grid << endl;
	cout << endl;

	cout << "---------------------------------------------------------------------------" << endl;

	getchar();




	corr_var_LOS = MatrixXReal::Zero(7, 7);
	corr_var_NLOS = MatrixXReal::Zero(6, 6);
	corr_var_OUT2IN = MatrixXReal::Zero(6, 6);
	Real num_LOS = 0;
	Real num_NLOS = 0;
	Real num_OUT2IN = 0;

	cout << "--------------------ue cross correlatin  ---------------------------------------" << endl;
	for (int bs_idx = 0; bs_idx < num_BS + num_mTRP; bs_idx++)
	{
		for (int ue_idx = 0; ue_idx < num_MS; ue_idx++)
		{
			if (channel[bs_idx][ue_idx].Propagation == LOS_propagation)
			{
				corr_var_LOS += ms[ue_idx].TLSPs[bs_idx].matrix() * ms[ue_idx].TLSPs[bs_idx].matrix().transpose();
				num_LOS++;
			}
			else if (channel[bs_idx][ue_idx].Propagation == NLOS_propagation)
			{
				corr_var_NLOS += ms[ue_idx].TLSPs[bs_idx].matrix() * ms[ue_idx].TLSPs[bs_idx].matrix().transpose();
				num_NLOS++;
			}
			else if (channel[bs_idx][ue_idx].Propagation == OUT2IN_propagation)
			{
				corr_var_OUT2IN += ms[ue_idx].TLSPs[bs_idx].matrix() * ms[ue_idx].TLSPs[bs_idx].matrix().transpose();
				num_OUT2IN++;
			}
			else
			{
				cout << "check Check_corr_var(), ue[ue_idx].propagation_condition[bs_idx]" << endl;
			}
		}
	}

	cout << "num_LOS : " << num_LOS << endl;
	cout << "var_LOS  " << endl;
	cout << corr_var_LOS / num_LOS << endl;
	cout << endl;
	cout << "num_NLOS : " << num_NLOS << endl;
	cout << "var_NLOS  " << endl;
	cout << endl;
	cout << corr_var_NLOS / num_NLOS << endl;
	cout << endl;
	if (num_OUT2IN != 0)
	{
		cout << "var_OUT2IN : " << num_OUT2IN << endl;
		cout << "var_OUT2IN " << endl;
		cout << corr_var_OUT2IN / num_OUT2IN << endl;
		cout << endl;
	}

	cout << "---------------------------------------------------------------------------" << endl;

	getchar();

}

void Transform_TLSP_to_LSP()
{
	for (int bs_idx = 0; bs_idx < num_BS + num_mTRP; bs_idx++)
	{
		for (int ue_idx = 0; ue_idx < num_MS; ue_idx++)
		{
			ms[ue_idx].LSPs[bs_idx](0) = channel[bs_idx][ue_idx].sigma_SF * ms[ue_idx].TLSPs[bs_idx](0); // Generate SF[dB]
			if (channel[bs_idx][ue_idx].Propagation == LOS_propagation)
			{
				ms[ue_idx].LSPs[bs_idx](LOS_K) = channel[bs_idx][ue_idx].sigma_K_factor * ms[ue_idx].TLSPs[bs_idx](LOS_K) + channel[bs_idx][ue_idx].mu_K_factor;// [dB]
				ms[ue_idx].LSPs[bs_idx](LOS_DS) = pow(10, channel[bs_idx][ue_idx].sigma_DS * ms[ue_idx].TLSPs[bs_idx](LOS_DS) + channel[bs_idx][ue_idx].mu_DS);
				//ms[ue_idx].LSPs[bs_idx](LOS_DS) = pow(10, normal(channel[bs_idx][ue_idx].mu_DS, channel[bs_idx][ue_idx].sigma_DS));
				ms[ue_idx].LSPs[bs_idx](LOS_ASD) = MIN(pow(10., channel[bs_idx][ue_idx].sigma_ASD * ms[ue_idx].TLSPs[bs_idx](LOS_ASD) + channel[bs_idx][ue_idx].mu_ASD), 104.);
				ms[ue_idx].LSPs[bs_idx](LOS_ASA) = MIN(pow(10., channel[bs_idx][ue_idx].sigma_ASA * ms[ue_idx].TLSPs[bs_idx](LOS_ASA) + channel[bs_idx][ue_idx].mu_ASA), 104.);
				ms[ue_idx].LSPs[bs_idx](LOS_ZSD) = MIN(pow(10., channel[bs_idx][ue_idx].sigma_ZSD * ms[ue_idx].TLSPs[bs_idx](LOS_ZSD) + channel[bs_idx][ue_idx].mu_ZSD), 52.);
				ms[ue_idx].LSPs[bs_idx](LOS_ZSA) = MIN(pow(10., channel[bs_idx][ue_idx].sigma_ZSA * ms[ue_idx].TLSPs[bs_idx](LOS_ZSA) + channel[bs_idx][ue_idx].mu_ZSA), 52.);
			}
			else if (channel[bs_idx][ue_idx].Propagation == NLOS_propagation)
			{
				ms[ue_idx].LSPs[bs_idx](NLOS_DS) = pow(10., channel[bs_idx][ue_idx].sigma_DS * ms[ue_idx].TLSPs[bs_idx](NLOS_DS) + channel[bs_idx][ue_idx].mu_DS);
				//ms[ue_idx].LSPs[bs_idx](NLOS_DS)  = pow(10, normal(channel[bs_idx][ue_idx].mu_DS, channel[bs_idx][ue_idx].sigma_DS));
				ms[ue_idx].LSPs[bs_idx](NLOS_ASD) = MIN(pow(10., channel[bs_idx][ue_idx].sigma_ASD * ms[ue_idx].TLSPs[bs_idx](NLOS_ASD) + channel[bs_idx][ue_idx].mu_ASD), 104.);
				ms[ue_idx].LSPs[bs_idx](NLOS_ASA) = MIN(pow(10., channel[bs_idx][ue_idx].sigma_ASA * ms[ue_idx].TLSPs[bs_idx](NLOS_ASA) + channel[bs_idx][ue_idx].mu_ASA), 104.);
				ms[ue_idx].LSPs[bs_idx](NLOS_ZSD) = MIN(pow(10., channel[bs_idx][ue_idx].sigma_ZSD * ms[ue_idx].TLSPs[bs_idx](NLOS_ZSD) + channel[bs_idx][ue_idx].mu_ZSD), 52.);
				ms[ue_idx].LSPs[bs_idx](NLOS_ZSA) = MIN(pow(10., channel[bs_idx][ue_idx].sigma_ZSA * ms[ue_idx].TLSPs[bs_idx](NLOS_ZSA) + channel[bs_idx][ue_idx].mu_ZSA), 52.);
				//cout << "Check 6" << endl;
			}
			else if (channel[bs_idx][ue_idx].Propagation == OUT2IN_propagation)
			{
				ms[ue_idx].LSPs[bs_idx](OUT2IN_DS) = pow(10., channel[bs_idx][ue_idx].sigma_DS * ms[ue_idx].TLSPs[bs_idx](OUT2IN_DS) + channel[bs_idx][ue_idx].mu_DS);
				//ms[ue_idx].LSPs[bs_idx](OUT2IN_DS) = pow(10, normal(channel[bs_idx][ue_idx].mu_DS, channel[bs_idx][ue_idx].sigma_DS));
				ms[ue_idx].LSPs[bs_idx](OUT2IN_ASD) = MIN(pow(10., channel[bs_idx][ue_idx].sigma_ASD * ms[ue_idx].TLSPs[bs_idx](OUT2IN_ASD) + channel[bs_idx][ue_idx].mu_ASD), 104.);
				ms[ue_idx].LSPs[bs_idx](OUT2IN_ASA) = MIN(pow(10., channel[bs_idx][ue_idx].sigma_ASA * ms[ue_idx].TLSPs[bs_idx](OUT2IN_ASA) + channel[bs_idx][ue_idx].mu_ASA), 104.);
				ms[ue_idx].LSPs[bs_idx](OUT2IN_ZSD) = MIN(pow(10., channel[bs_idx][ue_idx].sigma_ZSD * ms[ue_idx].TLSPs[bs_idx](OUT2IN_ZSD) + channel[bs_idx][ue_idx].mu_ZSD), 52.);
				ms[ue_idx].LSPs[bs_idx](OUT2IN_ZSA) = MIN(pow(10., channel[bs_idx][ue_idx].sigma_ZSA * ms[ue_idx].TLSPs[bs_idx](OUT2IN_ZSA) + channel[bs_idx][ue_idx].mu_ZSA), 52.);
			}
			else
			{
				cout << "check Transform_TLSP_to_LSP" << endl;
				getchar();
			}
		}
	}
}
