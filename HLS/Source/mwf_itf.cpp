#include "mainProject.hpp"

void mwf_itf(float gama, float beta, std::complex<float> initial_coeffs[12][256], std::complex<float> rv[NUM_CHAN][NUM_CHAN][256], std::complex<float> ry[NUM_CHAN][NUM_CHAN][256], std::complex<float> rx[NUM_CHAN][NUM_CHAN][256], int m, float qL[NUM_CHAN], float qR[NUM_CHAN])
{
	
	int ch = 6;
	int size_w_column = 256;
	std::complex<float> w_return[12][256];
	//~ std::vector<std::vector<std::complex<float>>> w_return_aux;
	
	std::complex<float> wL[6][256];
	std::complex<float> wR[6][256];
	
	std::complex<float> q[12];
	std::complex<float> ccL, ccR, itFin, c;
	std::complex<float> cL[6];
	std::complex<float> cR[6];
	std::complex<float> dL[6];
	std::complex<float> dR[6];
	float a,b,d,e;

	std::complex<float> w_ITF[12];
	
	for(int i=0; i<ch;i++)
		for(int j=0; j<size_w_column; j++)
			wL[i][j] = initial_coeffs[i][j];
	
	for(int i=0; i<ch;i++)
		for(int j=0; j<size_w_column; j++)
			wR[i][j] = initial_coeffs[i+ch][j];
	
	
	for(int i=0; i<ch;i++)
		q[i] = qL[i];
	
	for(int i=0; i<ch;i++)
		q[i+ch] = qR[i];
	
	std::complex<float> iMM[12][12];


	for (int i = 0; i < 12 ; ++i)
		for (int j = 0; j < 12 ; ++j)
			iMM[i][j] = 0;
	
	for (int i = 0; i < 2*ch ; ++i)
		iMM[i][i] = 1;
	
	std::complex<float> ry_y[6][6];
	std::complex<float> rx_x[6][6];
	std::complex<float> rv_v[6][6];
	
	std::complex<float> rxx[12][12];
	std::complex<float> ryy[12][12];
	
	for (int bin=0; bin<(m/2+1); bin++)
	{
		for (int i=0; i<ch; i++)
			for (int j=0; j<ch; j++)
			{
				ry_y[i][j]=ry[i][j][bin];
				rx_x[i][j]=rx[i][j][bin];
				rv_v[i][j]=rv[i][j][bin];
			}
		
		for (int i=0; i<12; i++)
			for (int j=0; j<12; j++)
			{
				if(i < ch)
				{
					if(j < ch)
					{
						rxx[i][j] = rx_x[i][j];
						ryy[i][j] = ry_y[i][j];
					}
					else
					{
						rxx[i][j] = 0;
						ryy[i][j] = 0;
					}
				}
				else
				{
					if(j < ch)
					{
						rxx[i][j] = 0;
						ryy[i][j] = 0;
					}
					else
					{
						rxx[i][j] = rx_x[i-ch][j-ch];
						ryy[i][j] = ry_y[i-ch][j-ch];
					}
				}
			}
		std::complex<float> v1[12];
		for (int i=0; i<12; i++)
		{
			std::complex<float> soma_linha(0,0);
			for (int j=0; j<12; j++)
				soma_linha = soma_linha + (beta * rxx[i][j] * q[j]);
				
			v1[i] = soma_linha;
		}
		
		std::complex<float> m1[12][12];
		for (int i=0; i<12; i++)
			for (int j=0; j<12; j++)
				m1[i][j] = iMM[i][j] - (beta*ryy[i][j]);
		
		std::complex<float> w_MWF[12];
		for (int i=0; i<12; i++)
		{
			std::complex<float> soma_linha(0,0);
			for (int j=0; j<12; j++)
				soma_linha = soma_linha + (m1[i][j] * initial_coeffs[j][bin]);
				
			soma_linha = soma_linha + v1[i];
			w_MWF[i] = soma_linha;
		}
		
		if (gama == 0)
		{
			for (int i=0; i<12; i++)
				w_return[i][bin] = w_MWF[i];
		}
		else
		{
			// ccL
			std::complex<float> res0[6];
			for (int j=0; j<6; j++)
			{
				std::complex<float> soma_linha(0,0);
				for (int i=0; i<ch; i++)
					soma_linha = soma_linha + (rv_v[j][i]*qR[i]);

				res0[j] = soma_linha;
			}

			std::complex<float> soma_linha_0(0,0);
			for (int j=0; j<ch; j++)
				soma_linha_0 = soma_linha_0 + (res0[j] * qL[j]);
			
			ccL = soma_linha_0;
			// Fim ccL

			//ccR

			std::complex<float> soma_linha_1(0,0);
			for (int j=0; j<ch; j++)
				soma_linha_1 = soma_linha_1 + (res0[j] * qR[j]);

			ccR = soma_linha_1;

			// Fim ccR

			itFin = ccL / ccR;


			for (int j=0; j<ch; j++)
			{
				std::complex<float> soma_linhaL(0,0);
				std::complex<float> soma_linhaR(0,0);
				for (int i=0; i<ch; i++)
				{
					soma_linhaL = soma_linhaL + (rv_v[j][i]*wL[i][bin]);
					soma_linhaR = soma_linhaR + (rv_v[j][i]*wR[i][bin]);
				}

				cL[j] = soma_linhaL;
				cR[j] = soma_linhaR;
			}

			std::complex<float> soma_linha_2_R(0,0);
			std::complex<float> soma_linha_2_L(0,0);
			for (int j=0; j<ch; j++)
			{
				soma_linha_2_L = soma_linha_2_L + (cL[j] * std::conj(wL[j][bin]));
				soma_linha_2_R = soma_linha_2_R + (cR[j] * std::conj(wR[j][bin]));
			}

			a = soma_linha_2_L.real();
			b = soma_linha_2_R.real();

			std::complex<float> soma_linha_3_c(0,0);
			for (int j=0; j<ch; j++)
				soma_linha_3_c = soma_linha_3_c + (cL[j] * std::conj(wR[j][bin]));

			c = soma_linha_3_c;

			d = 1/b;
			d = d * -d;

			for (int i = 0; i < ch; i++)
				dL[i] = std::complex<float>(d,0) * std::complex<float>(b,0) * (cL[i] - (std::conj(itFin)*cR[i]));

			std::complex<float> resitFin = itFin * c;
			e = 2 * (resitFin.real()) - a;

			for (int i = 0; i < ch; i++)
				dR[i] = std::complex<float>(d,0) * ((std::complex<float>(e,0)*cR[i]) - (std::complex<float>(b,0)*itFin*cL[i]));

			for (int i = 0; i < 12; i++)
			{
				if (i < ch)
					w_ITF[i] = dL[i];
				else
					w_ITF[i] = dR[i-ch];
			}

			for (int i=0; i<12; i++)
				w_return[i][bin] = w_MWF[i] + (std::complex<float>(gama,0)* w_ITF[i]);

		}

		if(bin > 0 && bin < m/2 )
		{
			for (int i=0; i<12; i++)
				w_return[i][m-bin] = std::conj(w_return[i][bin]);
		}
	}
	
	for(int i=0; i<12;i++)
		for(int j=0; j<size_w_column; j++)
			initial_coeffs[i][j] = w_return[i][j];	
}
