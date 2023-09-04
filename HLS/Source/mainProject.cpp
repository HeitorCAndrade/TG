#include "mainProject.hpp"

//void mainFilter(float gama, float beta, int & ct_frame, std::complex<float> wm, float weighting[256], int estimation_mode, float frame_in_tmp_in[256][6], int frame_vad[256], float & conty, float & contv, std::complex<float> (&rv)[6][6][256], std::complex<float> (&ry)[6][6][256], std::complex<float> (&rx)[6][6][256], float qL[6], float qR[6], std::complex<float> (&frame_out_tmp_in_aux)[256][2], std::complex<float> (&w)[12][256])
void mainFilter(float gama, float beta, float lambda, int BINS, int frameEndPriorEst, int & cont_frame, std::complex<float> wm, float weighting[256], int estimation_mode, float frame_in_tmp_in[256][NUM_CHAN], int frame_vad[256], int & cont_y, int & cont_v, std::complex<float> (&rv)[NUM_CHAN][NUM_CHAN][256], std::complex<float> (&ry)[NUM_CHAN][NUM_CHAN][256], std::complex<float> (&rx)[NUM_CHAN][NUM_CHAN][256], float qL[NUM_CHAN], float qR[NUM_CHAN], std::complex<float> (&frame_out_tmp_in_aux)[256][2], std::complex<float> (&w)[2*NUM_CHAN][256])
{


	//int length_fft = 256;
	////////////////////////////////
	//int frameEndPriorEst=1001;


	int channel = 6;
	float ctFrameTeste[256];
	std::complex<float> ctFrameTesteCompl[BINS];
	std::complex<float> phaseMod[256];
	std::complex<float> conjPhaseMod[256];
	cont_frame++;
    //teste 123...

	bool ovflo_all = false;
	static cmpxDataIn xn_input[256];
	static cmpxDataIn xn_input2[256];
	static cmpxDataOut xk_output[256];
	static cmpxDataOut xk_output2[256];
	//int frameEndPriorEst = 1170;
	//int frameEndPriorEst = 0;
	float frameMultPeso[NUM_CHAN][256];
	std::complex<float> frame_in_frq_in[256][NUM_CHAN];     //pode ser enviado do Testbench
	std::complex<float> frame_out_frq_in[256][2];
	std::complex<float> frame_out_frq_tmp[2][256];
	int sumVad = 0;
	std::complex<float> zeroComplex(0,0);
	//float lambda=0.9998;
	std::complex<float> normalization = std::complex<float>(256, 0);


	for (int i = 0; i < 256; i++)
		ctFrameTeste[i] = cont_frame*64*i;

	for (int i = 0; i < 256; i++){
		std::complex<float> numComplex = wm;
		ctFrameTesteCompl[i] = pow(wm,(int)-ctFrameTeste[i]);

		if (((i+1) % 2 == 0) && (cont_frame % 2 != 0))
			ctFrameTesteCompl[i] = ctFrameTesteCompl[i] * std::complex<float>(-1,0);
	}
   //Cálculo maluco de Phasemod
	for (int i = 0; i < 256; i++)
	{
		std::complex<float> numb = pow(-1,i);
		std::complex<float> numbb = numb*ctFrameTesteCompl[i];
		float real = numbb.real();
		float imag = numbb.imag();
		if (real > 0)
		{
			if (real > 0.5)
				real = 1.00000000000000;
			else
				real = 0.00000000000000;
		}
		else
		{
			if (real > -0.5)
				real = 0.00000000000000;
			else
				real = -1.00000000000000;
		}

		if (imag > 0)
		{
			if (imag > 0.5)
				imag = 1.00000000000000;
			else
				imag = 0.00000000000000;
		}
		else
		{
			if (imag > -0.5)
				imag = 0.00000000000000;
			else
				imag = -1.00000000000000;
		}

		phaseMod[i] = std::complex<float>(real,imag);
	}
////////////////////////////////////////////////////////////////////////////
	for (int j = 0; j < BINS; j++)
		for (int i = 0; i < NUM_CHAN; i++)
			frameMultPeso[i][j] = weighting[j] * frame_in_tmp_in[j][i];

////////////////////////////// Laço da FFT //////////////////////////////////

	for(int ch = 0; ch < NUM_CHAN; ch++)
	{
		for (int i = 0; i < BINS; i++)
			xn_input[i] = cmpxDataIn((float)frameMultPeso[ch][i], 0);

		bool ovflo;

		fft_top(1, xn_input, xk_output, &ovflo);

		ovflo_all |= ovflo;
		for (int i=0; i<BINS; i++)
			frame_in_frq_in[i][ch] = xk_output[i] * phaseMod[i];
	}

////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////

	if (cont_frame <= frameEndPriorEst)  //Só entra aqui no período de estimação PRIOR
	{
		std::complex<float> sumFrqIn = 0;
		for (int i=0; i<BINS; i++)
			for (int j=0; j<NUM_CHAN; j++)
				sumFrqIn = sumFrqIn + frame_in_frq_in[i][j];
			
		if (sumFrqIn.real() != zeroComplex.real() || sumFrqIn.imag() != zeroComplex.imag())  // Gustavo usou operador lógico OR
		{
			sumVad = 0;
			for (int i=0; i<BINS; i++)
				sumVad = sumVad + frame_vad[i];  //soma todos os resultados do VAD dentro do Frame
			
			corr_matrix_estimation(0, BINS, NUM_CHAN, frame_in_frq_in, rv, ry, rx, sumVad, cont_y, cont_v, lambda);
		}
		
		for (int i=0; i<BINS; i++)
		{
			std::complex<float> soma_linha0(0,0);
			std::complex<float> soma_linha1(0,0);
			for (int j=0; j<NUM_CHAN; j++)
			{
				soma_linha0 = soma_linha0 + (frame_in_frq_in[i][j] * qL[j]);
				soma_linha1 = soma_linha1 + (frame_in_frq_in[i][j] * qR[j]);
			}
			frame_out_frq_in[i][0] = soma_linha0;
			frame_out_frq_in[i][1] = soma_linha1;
		}
		
		if (cont_frame == frameEndPriorEst)
		{
			for(int bin = 0; bin < NUM_CHAN; bin++)
				for(int i = 0; i < NUM_CHAN; ++i)
					for(int j = 0; j < BINS; ++j)
					{
						rv[bin][i][j] = rv[bin][i][j] / (float)cont_v;
						ry[bin][i][j] = ry[bin][i][j] / (float)cont_y;
						rx[bin][i][j] = ry[bin][i][j] - rv[bin][i][j];
					}
		}
	}			
	else //Aqui entra para começar a filtragem com coeficientes verdadeiros
	{
		sumVad = 0;
		for (int i=0; i<BINS; i++)
			sumVad = sumVad + frame_vad[i];

		corr_matrix_estimation(1, BINS, NUM_CHAN, frame_in_frq_in, rv, ry, rx, sumVad, cont_y, cont_v, lambda);
		
		mwf_itf(gama, beta, w, rv, ry, rx, 256, qL, qR);
		
		for (int i=0; i<BINS; i++)
		{
			std::complex<float> soma_linha0(0,0);
			std::complex<float> soma_linha1(0,0);
			for (int j=0; j<NUM_CHAN; j++)
			{
				std::complex<float> a0 = w[j][i];
				std::complex<float> a1 = w[j+NUM_CHAN][i];
				std::complex<float> b0 = frame_in_frq_in[i][j];

				std::complex<float> res0 = frame_in_frq_in[i][j] * std::conj(w[j][i]);
				soma_linha0 = soma_linha0 + res0;
				std::complex<float> res1 = frame_in_frq_in[i][j] * std::conj(w[j+NUM_CHAN][i]);
				soma_linha1 = soma_linha1 + res1;
			}
			frame_out_frq_in[i][0] = soma_linha0;
			frame_out_frq_in[i][1] = soma_linha1;
		}
	}


    for (int i = 0; i < BINS; i++)
    {
    	std::complex<float> numbb = std::conj(phaseMod[i]);
    	float real = numbb.real();
    	float imag = numbb.imag();
    	if (real > 0)
    	{
    		if (real > 0.5)
    			real = 1.00000000000000;
    		else
    			real = 0.00000000000000;
    	}
    	else
    	{
    		if (real > -0.5)
				real = 0.00000000000000;
			else
				real = -1.00000000000000;
    	}

    	if (imag > 0)
    	{
    		if (imag > 0.5)
    			imag = 1.00000000000000;
    		else
    			imag = 0.00000000000000;
    	}
    	else
    	{
    		if (imag > -0.5)
    			imag = 0.00000000000000;
			else
				imag = -1.00000000000000;
    	}

    	conjPhaseMod[i] = std::complex<float>(real,imag);
    }



    for (int i=0; i<BINS; i++)
    	for (int j=0; j<2; j++)
    		frame_out_frq_tmp[j][i] = conjPhaseMod[i] * frame_out_frq_in[i][j];


	for (int j=0; j<2; j++)
	{

		bool ovflo;

		for (int i=0; i<BINS; i++)
			xn_input2[i] = cmpxDataIn(frame_out_frq_tmp[j][i].real(),frame_out_frq_tmp[j][i].imag());

		fft_top(0, xn_input2, xk_output2, &ovflo);    //FFT inversa

		for (int i=0; i<BINS; i++)
			xk_output2[i] = xk_output2[i] / normalization;

		for (int i=0; i<BINS; i++)
		{
			std::complex<float> numB = std::complex<float>((float)weighting[i],0);
			std::complex<float> mult = xk_output2[i] * numB;
			frame_out_tmp_in_aux[i][j] = mult;
		}
	}
}


void mainFilterForResults(int cont_frame , int frameEndPriorEst, std::complex<float> wm, float weighting[256], float frame_in_tmp_in[256][NUM_CHAN], float qL[NUM_CHAN], float qR[NUM_CHAN], std::complex<float> (&frame_out_frq_out_puro)[256][2], std::complex<float> (&frame_out_frq_out_filtrado)[256][2], std::complex<float> w[12][256])
{
	int length_fft = 256;
	int channel = 6;
	float ctFrameTeste[256];
	std::complex<float> ctFrameTesteCompl[256];
	std::complex<float> phaseMod[256];
	std::complex<float> conjPhaseMod[256];


	bool ovflo_all = false;
	static cmpxDataIn xn_input[256];
	static cmpxDataIn xn_input2[256];
	static cmpxDataOut xk_output[256];
	static cmpxDataOut xk_output2[256];
	//int frameEndPriorEst = 1501;
	float frameMultPeso[6][256];
	std::complex<float> frame_in_frq_in[256][6];
	std::complex<float> frame_out_frq_tmp[2][256];
	int sumVad = 0;
	std::complex<float> zeroComplex(0,0);
	float lambda=0.9998;
	std::complex<float> normalization = std::complex<float>(256, 0);

	for (int i = 0; i < 256; i++)
		ctFrameTeste[i] = cont_frame*64*i;

	for (int i = 0; i < 256; i++){
		std::complex<float> numComplex = wm;
		ctFrameTesteCompl[i] = pow(wm,(int)-ctFrameTeste[i]);

		if (((i+1) % 2 == 0) && (cont_frame % 2 != 0))
			ctFrameTesteCompl[i] = ctFrameTesteCompl[i] * std::complex<float>(-1,0);
	}

	for (int i = 0; i < 256; i++)
	{
		std::complex<float> numb = pow(-1,i);
		std::complex<float> numbb = numb*ctFrameTesteCompl[i];
		float real = numbb.real();
		float imag = numbb.imag();
		if (real > 0)
		{
			if (real > 0.5)
				real = 1.00000000000000;
			else
				real = 0.00000000000000;
		}
		else
		{
			if (real > -0.5)
				real = 0.00000000000000;
			else
				real = -1.00000000000000;
		}

		if (imag > 0)
		{
			if (imag > 0.5)
				imag = 1.00000000000000;
			else
				imag = 0.00000000000000;
		}
		else
		{
			if (imag > -0.5)
				imag = 0.00000000000000;
			else
				imag = -1.00000000000000;
		}

		phaseMod[i] = std::complex<float>(real,imag);
	}

	for (int j = 0; j < 256; j++)
		for (int i = 0; i < 6; i++)
			frameMultPeso[i][j] = weighting[j] * frame_in_tmp_in[j][i];

	for(int ch = 0; ch < 6; ch++)
	{
		for (int i = 0; i < 256; i++)
			xn_input[i] = cmpxDataIn((float)frameMultPeso[ch][i], 0);

		bool ovflo;

		fft_top(1, xn_input, xk_output, &ovflo);

		ovflo_all |= ovflo;
		for (int i=0; i<length_fft; i++)
			frame_in_frq_in[i][ch] = xk_output[i] * phaseMod[i];
	}

	for (int i=0; i<length_fft; i++)
	{
		std::complex<float> soma_linha0(0,0);
		std::complex<float> soma_linha1(0,0);
		for (int j=0; j<channel; j++)
		{
			soma_linha0 = soma_linha0 + (frame_in_frq_in[i][j] * qL[j]);
			soma_linha1 = soma_linha1 + (frame_in_frq_in[i][j] * qR[j]);
		}
		frame_out_frq_out_puro[i][0] = soma_linha0;
		frame_out_frq_out_puro[i][1] = soma_linha1;
	}

	for (int i=0; i<length_fft; i++)
	{
		std::complex<float> soma_linha0(0,0);
		std::complex<float> soma_linha1(0,0);
		for (int j=0; j<channel; j++)
		{
			std::complex<float> res0 = frame_in_frq_in[i][j] * std::conj(w[j][i]);
			soma_linha0 = soma_linha0 + res0;
			std::complex<float> res1 = frame_in_frq_in[i][j] * std::conj(w[j+6][i]);
			soma_linha1 = soma_linha1 + res1;
		}
		frame_out_frq_out_filtrado[i][0] = soma_linha0;
		frame_out_frq_out_filtrado[i][1] = soma_linha1;
	}
}
