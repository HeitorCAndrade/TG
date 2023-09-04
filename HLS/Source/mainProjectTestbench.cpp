/*******************************************************************************
Vendor: Xilinx 
Associated Filename: fft_tb.cpp
Purpose: Xilinx FFT IP-XACT IP in Vivado HLS
Revision History: September 26, 2013 - initial release
                                                
*******************************************************************************
#-  (c) Copyright 2011-2018 Xilinx, Inc. All rights reserved.
#-
#-  This file contains confidential and proprietary information
#-  of Xilinx, Inc. and is protected under U.S. and
#-  international copyright and other intellectual property
#-  laws.
#-
#-  DISCLAIMER
#-  This disclaimer is not a license and does not grant any
#-  rights to the materials distributed herewith. Except as
#-  otherwise provided in a valid license issued to you by
#-  Xilinx, and to the maximum extent permitted by applicable
#-  law: (1) THESE MATERIALS ARE MADE AVAILABLE "AS IS" AND
#-  WITH ALL FAULTS, AND XILINX HEREBY DISCLAIMS ALL WARRANTIES
#-  AND CONDITIONS, EXPRESS, IMPLIED, OR STATUTORY, INCLUDING
#-  BUT NOT LIMITED TO WARRANTIES OF MERCHANTABILITY, NON-
#-  INFRINGEMENT, OR FITNESS FOR ANY PARTICULAR PURPOSE; and
#-  (2) Xilinx shall not be liable (whether in contract or tort,
#-  including negligence, or under any other theory of
#-  liability) for any loss or damage of any kind or nature
#-  related to, arising under or in connection with these
#-  materials, including for any direct, or any indirect,
#-  special, incidental, or consequential loss or damage
#-  (including loss of data, profits, goodwill, or any type of
#-  loss or damage suffered as a result of any action brought
#-  by a third party) even if such damage or loss was
#-  reasonably foreseeable or Xilinx had been advised of the
#-  possibility of the same.
#-
#-  CRITICAL APPLICATIONS
#-  Xilinx products are not designed or intended to be fail-
#-  safe, or for use in any application requiring fail-safe
#-  performance, such as life-support or safety devices or
#-  systems, Class III medical devices, nuclear facilities,
#-  applications related to the deployment of airbags, or any
#-  other applications that could lead to death, personal
#-  injury, or severe property or environmental damage
#-  (individually and collectively, "Critical
#-  Applications"). Customer assumes the sole risk and
#-  liability of any use of Xilinx products in Critical
#-  Applications, subject only to applicable laws and
#-  regulations governing limitations on product liability.
#-
#-  THIS COPYRIGHT NOTICE AND DISCLAIMER MUST BE RETAINED AS
#-  PART OF THIS FILE AT ALL TIMES. 
#- ************************************************************************


This file contains confidential and proprietary information of Xilinx, Inc. and 
is protected under U.S. and international copyright and other intellectual 
property laws.

DISCLAIMER
This disclaimer is not a license and does not grant any rights to the materials 
distributed herewith. Except as otherwise provided in a valid license issued to 
you by Xilinx, and to the maximum extent permitted by applicable law: 
(1) THESE MATERIALS ARE MADE AVAILABLE "AS IS" AND WITH ALL FAULTS, AND XILINX 
HEREBY DISCLAIMS ALL WARRANTIES AND CONDITIONS, EXPRESS, IMPLIED, OR STATUTORY, 
INCLUDING BUT NOT LIMITED TO WARRANTIES OF MERCHANTABILITY, NON-INFRINGEMENT, OR 
FITNESS FOR ANY PARTICULAR PURPOSE; and (2) Xilinx shall not be liable (whether 
in contract or tort, including negligence, or under any other theory of 
liability) for any loss or damage of any kind or nature related to, arising under 
or in connection with these materials, including for any direct, or any indirect, 
special, incidental, or consequential loss or damage (including loss of data, 
profits, goodwill, or any type of loss or damage suffered as a result of any 
action brought by a third party) even if such damage or loss was reasonably 
foreseeable or Xilinx had been advised of the possibility of the same.

CRITICAL APPLICATIONS
Xilinx products are not designed or intended to be fail-safe, or for use in any 
application requiring fail-safe performance, such as life-support or safety 
devices or systems, Class III medical devices, nuclear facilities, applications 
related to the deployment of airbags, or any other applications that could lead 
to death, personal injury, or severe property or environmental damage 
(individually and collectively, "Critical Applications"). Customer assumes the 
sole risk and liability of any use of Xilinx products in Critical Applications, 
subject only to applicable laws and regulations governing limitations on product 
liability. 

THIS COPYRIGHT NOTICE AND DISCLAIMER MUST BE RETAINED AS PART OF THIS FILE AT 
ALL TIMES.

*******************************************************************************/

# define M_PI           3.14159265358979323846  /* pi */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <sstream>
//#include "fft_top.h"
#include "mainProject.hpp"
#include "Experiment_Constraints.hpp"
#include <stdio.h>

#ifdef _WIN32
#include <Windows.h>
#else
#include <unistd.h>
#endif

#include <fstream>
#include <string>
#include <sstream>

//#include <math.h>
using namespace std;



//Escolhe se pretende gerar as métricas
//#define METRICAS

float sciToDub(const string& str) {

   stringstream ss(str);
   float d = 0;
   ss >> d;

   if (ss.fail()) {
      string s = "Unable to format ";
      s += str;
      s += " as a number!";
      throw (s);
   }

   return (d);
}

//bool compare_float(float x, float y, float epsilon = 0.001f){
//   if(fabs(x - y) < epsilon)
//      return true; //they are same
//      return false; //they are not same
//}


int main()
{

	std::cout.precision(20);

    //Variáveis declaradas aqui
	std::string line, token;
	std::ifstream infilePeso("peso.txt");

    //Ponteiros para arquivos de entradas
	ifstream* inputAudio = new ifstream;
	ifstream* inputAudioVAD = new ifstream;

	struct Data_Experiments Experiment;

	std::complex<float> wm(0.999698818696204,0.024541228522912);

	#ifdef METRICAS


	std::complex<float> frame_out_frq_out_puroFala[256][2];         //Colocar essas variáveis no .hpp   - frame_in_frq_sp
	std::complex<float> frame_out_frq_out_filtradoFala[256][2];     //Colocar essas variáveis no .hpp   - frame_out_frq_sp
	std::complex<float> frame_out_frq_out_puroRuido[256][2];        //Colocar essas variáveis no .hpp   - frame_in_frq_no
	std::complex<float> frame_out_frq_out_filtradoRuido[256][2];    //Colocar essas variáveis no .hpp   - frame_out_frq_no

	//Original Input Speech and Noise audio files
	std::ifstream inputSpeechAudio(Input_Speech_audioFile);
	std::ifstream inputNoiseAudio(Input_Noise_audioFile);


	//vetores de média de N frames passados (usado nas métricas)
	std::complex<float> frames_speech_in_L[BINS][NS];
	std::complex<float> frames_speech_in_R[BINS][NS];
	std::complex<float> frames_noise_in_L[BINS][NS];
	std::complex<float> frames_noise_in_R[BINS][NS];

	std::complex<float> frames_speech_out_L[BINS][NS];
    std::complex<float> frames_speech_out_R[BINS][NS];
	std::complex<float> frames_noise_out_L[BINS][NS];
	std::complex<float> frames_noise_out_R[BINS][NS];

	std::complex<float> pvL_k[BINS];
	std::complex<float> pvR_k[BINS];
	std::complex<float> pvLR_k[BINS];

	std::complex<float> pxL_k[BINS];
	std::complex<float> pxR_k[BINS];
	std::complex<float> pxLR_k[BINS];

	std::complex<float> pzvL_k[BINS];
	std::complex<float> pzvR_k[BINS];
	std::complex<float> pzvLR_k[BINS];

	std::complex<float> pzxL_k[BINS];
	std::complex<float> pzxR_k[BINS];
	std::complex<float> pzxLR_k[BINS];

	//SNR Data
	    float snrPuroL = 0;
	    float snrPuroR = 0;
	    float snrFiltradoL = 0;
	    float snrFiltradoR = 0;
	    float Delta_SNRg[BINS];
	    float Delta_ILDv[BINS];
	    float Delta_ILDx[BINS];
	    float Delta_IPDv[BINS];
	    float Delta_IPDx[BINS];



	    //ILD
	    float ild_sp = 0;
	    float ild_no = 0;
	    float Delta_ild_sp = 0;
	    float Delta_ild_no = 0;
    #endif

	int countline = 0;
    //Abre o arquivo com os pesos da WOLA
	while (std::getline(infilePeso, line)){
		std::istringstream iss(line);
		if(countline>BINS-1)
			exit(0);

		std::getline(iss, token, '\n');
		float numberDub = sciToDub(token);
		weighting[countline] = numberDub;
		countline++;
	}

	/////////////////////////////////////////////////////////////////////////////////////////
	// Carrega arquivos de entrada

	for (int numTeste = 0; numTeste < number_of_experiments; numTeste++){
		for (int countGama = 0; countGama < sizeGama; countGama++){

			//Preparação dos arquivos de saída
			ofstream output_signal; //Classe para escrita de arquivos (Arquivo de saída)
			output_signal.open (Output_File_Name[countGama]);
            #ifdef METRICAS
				//std::ofstream deltaSNR;
				std::ofstream power_vectors;
				std::ofstream debugFFT;
				std::ofstream GLOBAL_MEASURES;

				//deltaSNR.open (Delta_SNR[countGama]);
				//deltaILD.open(Delta_ILD[countGama]);
				power_vectors.open("Power_vectors.txt");
				debugFFT.open(debug_FFT[countGama]);
				GLOBAL_MEASURES.open(MEASURES[countGama]);

            #endif
			//ofstream wSaida; //Classe para escrita de arquivos
			//wSaida.open (Output_Coefficients[countGama]);  //Arquivo de coeficientes
			////////////////////////////////////////////////////
			//Teste salva Rv no bin 10
			//ofstream Rv_Matrix; //Classe para escrita de arquivos
			//Rv_Matrix.open("Rv_Vitis.txt");
			//////////////////////////////////////////////////////

            ///////////////////////////////////////////////////////////////
			//Teste frame_in_tmp_in
			//ofstream frame_entrada_tempo;
			//frame_entrada_tempo.open("Frame_Input_Tempo.txt");

		    inputAudio->open(testBench_audioFile);  //Carrega arquivo com áudio de entrada .hpp
			inputAudioVAD->open(testBench_VAD);     //Carrega arquivo com resultado do VAD .hpp


		// INICIALIZACAO

            // Inicializa as matrizes de estimação com ZERO
			for(int x = 0; x < NUM_MICS; x++)
			{
				for(int y = 0; y < NUM_MICS; y++)
				{
					for(int z = 0; z < BINS; z++)
					{
						std::complex<float> entry(0,0);
						Rv[x][y][z] = entry;
						Rx[x][y][z] = entry;
						Ry[x][y][z] = entry;
					}
				}
			}
			//////////////////////////////////////////////////////////////////////////////
			// FIM INICIALIZACAO

			// IMPORTANTES
			float input[NUM_MICS];
			float realPart=0.0, imagPart=0.0;

			// FIM MEU CODIGO
			int estimation_mode = 0;


			//COEFICIENTS FOR TRIVIAL FILTER (ALL PASS FILTER) APPLIED ON PRIOR_ESTIMATION_TIME

			//test
			std::complex<float> w[12][256];
			///

			for (int i = 0; i < 256; i++)
			{
				 w[0][i] = 1;
				 w[9][i] = 1;
				 //w[11][i] = 5;

			}


			int countVAD = 0;


			//Início do carregamento das amostras

			for (cont_samples = 0; cont_samples < tamSign; cont_samples++)
			{
				for (int k =0; k < BINS-1; k++){
					for (int g =0; g < NUM_MICS; g++){
						frame_in_tmp_in[k][g] = frame_in_tmp_in[k+1][g];
                        #ifdef METRICAS
				        	frame_in_tmp_sp[k][g] = frame_in_tmp_sp[k+1][g];   // Monta frame de speech
				        	frame_in_tmp_no[k][g] = frame_in_tmp_no[k+1][g];   // Monta frame de noise
                        #endif
					}}
////////////////////////////////   INPUT   /////////////////////////////////////////////////////////////
				std::getline(*inputAudio, line);
				std::istringstream iss(line);
				for (int i = 0; i < NUM_MICS; i++) {
					std::getline(iss, token, '\t');
					float numberDub = sciToDub(token);
					input[i] = numberDub;
				}
				for (int g = 0; g < NUM_MICS; g++)
					frame_in_tmp_in[BINS-1][g] = input[g];

				for (int k =0; k < BINS-1; k++)          //Esse buffer é usado para a IFFT
					for (int g =0; g < 2; g++)
						frame_out_tmp_in[k][g] = frame_out_tmp_in[k+1][g];

				for (int g =0; g < 2; g++)
					frame_out_tmp_in[BINS-1][g] = 0;
//////////////////////////////////  SPEECH    //////////////////////////////////////////////////////////
             #ifdef METRICAS
				std::getline(inputSpeechAudio, line); //Load data from Speech File
				std::istringstream iss2(line);
				for (int i = 0; i < NUM_MICS; i++) {
					std::getline(iss2, token, '\t');
					float numberDub = sciToDub(token);
				    input[i] = numberDub; }
				for (int g = 0; g < NUM_MICS; g++)
					frame_in_tmp_sp[BINS-1][g] = input[g]; //New SPEECH sample get in (BINS-1) position

//////////////////////////////////  NOISE    //////////////////////////////////////////////////////////
				std::getline(inputNoiseAudio, line);  //Load data from Noise File
				std::istringstream iss3(line);
				for (int i = 0; i < NUM_MICS; i++) {
					std::getline(iss3, token, '\t');
					float numberDub = sciToDub(token);
					input[i] = numberDub; }
				for (int g = 0; g < NUM_MICS; g++)
					frame_in_tmp_no[BINS-1][g] = input[g]; //New NOISE sample get in (BINS-1) position

            #endif
////////////////////////////////  VAD FILE  ///////////////////////////////////////////////////////////////
				//Início do carregamento do arquivo de VAD
				for(int i = BINS-1; i > 0; i--) {
					frame_vad[i] = frame_vad[i-1];
				}

				std::getline(*inputAudioVAD, line);
				std::istringstream iss4(line);
				std::getline(iss4, token, '\n');
				int numberDubInt = (int)sciToDub(token);
				frame_vad[0] = numberDubInt;
////////////////////////////////////////////////////////////////////////////////////////////////////////////

				if (((cont_samples+1) % BLOCK) == 0)
				{

					mainFilter(gama[countGama], beta, lambda, BINS, frameEndPriorEst, cont_frames, wm, weighting, estimation_mode, frame_in_tmp_in, frame_vad, cont_y, cont_v, Rv, Ry, Rx, qL, qR, frame_out_tmp_in_aux, w);


                    #ifdef METRICAS

					  //Não carregar coeficientes de arquivo, pegar direto o w calculado por mainfilter

						 //Calcula números do ruído
					 	mainFilterForResults(cont_frames,frameEndPriorEst, wm, weighting, frame_in_tmp_no, qL, qR, frame_out_frq_out_puroRuido, frame_out_frq_out_filtradoRuido, w);


					     //Calcula números da fala
						mainFilterForResults(cont_frames,frameEndPriorEst, wm, weighting, frame_in_tmp_sp, qL, qR, frame_out_frq_out_puroFala, frame_out_frq_out_filtradoFala, w);

					for (int i=0; i<NS-1; i++){
					 for (int k=0; k<BINS; k++){
						 frames_speech_in_L[k][i]=frames_speech_in_L[k][i+1];
						 frames_speech_in_R[k][i]=frames_speech_in_R[k][i+1];
						 frames_noise_in_L[k][i]=frames_noise_in_L[k][i+1];
						 frames_noise_in_R[k][i]=frames_noise_in_R[k][i+1];

						 frames_speech_out_L[k][i]=frames_speech_out_L[k][i+1];
						 frames_speech_out_R[k][i]=frames_speech_out_R[k][i+1];
						 frames_noise_out_L[k][i]=frames_noise_out_L[k][i+1];
						 frames_noise_out_R[k][i]=frames_noise_out_R[k][i+1];
					 }}
					for (int g = 0; g < BINS; g++){
						frames_speech_in_L[g][NS-1] = frame_out_frq_out_puroFala[g][0]; //New NOISE sample get in (BINS-1) position
						frames_speech_in_R[g][NS-1] =frame_out_frq_out_puroFala[g][1];
						frames_noise_in_L[g][NS-1] =frame_out_frq_out_puroRuido[g][0];
						frames_noise_in_R[g][NS-1] =frame_out_frq_out_puroRuido[g][1];

						frames_speech_out_L[g][NS-1] = frame_out_frq_out_filtradoFala[g][0]; //New NOISE sample get in (BINS-1) position
						frames_speech_out_R[g][NS-1] =frame_out_frq_out_filtradoFala[g][1];
						frames_noise_out_L[g][NS-1] =frame_out_frq_out_filtradoRuido[g][0];
						frames_noise_out_R[g][NS-1] =frame_out_frq_out_filtradoRuido[g][1];
					}
					if (cont_frames >= frameEndPriorEst ){  //começa contabilizar métricas

						/*for (int wCount = 0; wCount < 256; wCount++){ //Salva coeficientes
								for (int wCountMic = 0; wCountMic < 12; wCountMic++){
									wSaida <<  std::fixed << std::setprecision(16) << w[wCountMic][wCount].real() << " " << w[wCountMic][wCount].imag() << "\t";}
								wSaida << "\n";
							}
						//int f=5;
						for(int linha=0; linha<=5; linha ++){
							for(int contador=0; contador<12; contador ++){
								Rv_Matrix << std::fixed << std::setprecision(16) << Rv[linha][contador][5].real() << " " << Rv[linha][contador][5].imag() << "\t";

							}
							Rv_Matrix << "\n";
						}
                        */
						//for(int linha=0; linha<=255; linha ++){
						//	for(int coluna=0; coluna<=5; coluna ++){
						//		frame_entrada_tempo <<  std::fixed << std::setprecision(16) << frame_in_tmp_in[linha][coluna] << "\t";
						//	 }
						//	frame_entrada_tempo << "\n";
						//}

						for (int k=0; k<BINS;k++){  //zera vetores acumuladores de power
							pvL_k[k]=1e-12;
							pvR_k[k]=1e-12;
							pvLR_k[k]=1e-12;
							pxL_k[k]=1e-12;
							pxR_k[k]=1e-12;
							pxLR_k[k]=1e-12;
							pzvL_k[k]=1e-12;
						    pzvR_k[k]=1e-12;
							pzvLR_k[k]=1e-12;
							pzxL_k[k]=1e-12;
						    pzxR_k[k]=1e-12;
							pzxLR_k[k]=1e-12;
						}

						for (int k=0; k<BINS;k++){
							for (int i=0; i<NS;i++){
								//Noise power
								pvL_k[k]=pvL_k[k]+frames_noise_in_L[k][i]; //ajustar conta de média
								pvR_k[k]=pvR_k[k]+frames_noise_in_R[k][i];
								pvLR_k[k]=pvLR_k[k]+(frames_noise_in_L[k][i]*std::conj(frames_noise_in_R[k][i]));
								//Speech Power
								pxL_k[k]=pxL_k[k]+frames_speech_in_L[k][i];
								pxR_k[k]=pxR_k[k]+frames_speech_in_R[k][i];
								pxLR_k[k]=pxLR_k[k]+(frames_speech_in_L[k][i]* std::conj(frames_speech_in_R[k][i]));
								//Noise Filtered Power
								pzvL_k[k]=pzvL_k[k]+frames_noise_out_L[k][i]; //ajustar conta de média
							    pzvR_k[k]=pzvR_k[k]+frames_noise_out_R[k][i];
							    pzvLR_k[k]=pzvLR_k[k]+(frames_noise_out_L[k][i]* std::conj(frames_noise_out_R[k][i])); //Cross-Power
							    //Speech Filtered Power
							    pzxL_k[k]=pzxL_k[k]+frames_speech_out_L[k][i]; //ajustar conta de média
							    pzxR_k[k]=pzxR_k[k]+frames_speech_out_R[k][i];
							    pzxLR_k[k]=pzxLR_k[k]+(frames_speech_out_L[k][i]* std::conj(frames_speech_out_R[k][i])); //Cross-Power
							}
						}
						/////// Médias e cálculos de power

						 for (int l=0; l<BINS; l++){
							 pvL_k[l]=pow((std::abs(pvL_k[l])/(std::complex<float>)NS),2);
							 pvR_k[l]=pow((std::abs(pvR_k[l])/(std::complex<float>)NS),2);
							 pvLR_k[l]=pvLR_k[l]/(std::complex<float>)NS;

							 pxL_k[l]=pow((std::abs(pxL_k[l])/(std::complex<float>)NS),2);
							 pxR_k[l]=pow((std::abs(pxR_k[l])/(std::complex<float>)NS),2);
							 pxLR_k[l]=pxLR_k[l]/(std::complex<float>)NS;

							 pzvL_k[l]=pow((std::abs(pzvL_k[l])/(std::complex<float>)NS),2);
							 pzvR_k[l]=pow((std::abs(pzvR_k[l])/(std::complex<float>)NS),2);
							 pzvLR_k[l]=pzvLR_k[l]/(std::complex<float>)NS;

							 pzxL_k[l]=pow((std::abs(pzxL_k[l])/(std::complex<float>)NS),2);
							 pzxR_k[l]=pow((std::abs(pzxR_k[l])/(std::complex<float>)NS),2);
							 pzxLR_k[l]=pzxLR_k[l]/(std::complex<float>)NS;

							 power_vectors << std::fixed << std::setprecision(16) << pvLR_k[l] << "\t" << pxLR_k[l] << "\t" << pzvLR_k[l]  << "\t" << pzxLR_k[l] << "\n";
						 }



						 /////// Métricas (ERRO!!! vetores power são números complexos)

						 for (int p=0;p<BINS;p++){

							 Delta_SNRg[p] =std::abs(10*log10((std::abs(pzxL_k[p]) + std::abs(pzxR_k[p]))/std::abs((pzvL_k[p]) + std::abs(pzvR_k[p])))-10*log10((std::abs(pxL_k[p])+std::abs(pxR_k[p]))/(std::abs(pvL_k[p])+std::abs(pvR_k[p]))));
							 //Delta_SNRg[p] =std::abs(10*log10((std::abs(pzxL_k[p] + pzxR_k[p]))/std::abs((pzvL_k[p] + pzvR_k[p])))-10*log10((std::abs(pxL_k[p]+pxR_k[p]))/(std::abs(pvL_k[p]+pvR_k[p]))));
							 Delta_ILDv[p]=std::abs(10*log10(std::abs(pzvL_k[p])/std::abs(pzvR_k[p]))-10*log10(std::abs(pvL_k[p])/std::abs(pvR_k[p])));
							 Delta_ILDx[p]=std::abs(10*log10(std::abs(pzxL_k[p])/std::abs(pzxR_k[p]))-10*log10(std::abs(pxL_k[p])/std::abs(pxR_k[p])));
							 Delta_IPDv[p]=std::abs(std::arg(pzvLR_k[p])-std::arg(pvLR_k[p]));
							 Delta_IPDx[p]=std::abs(std::arg(pzxLR_k[p])-std::arg(pxLR_k[p]));
						 }

						 DSNR=0;  //Variáveis que acumulam devem ser zeradas e serão gravadas no arquivo texto
					     DIPDv=0;
					     DILDv=0;
					     DIPDx=0;
					     DILDx=0;

						 for (int p=0;p<BINS;p++){  //Média dos bins
						      DSNR=DSNR+(Delta_SNRg[p]);
						      DILDv=DILDv+(Delta_ILDv[p]);
						      DILDx=DILDx+(Delta_ILDx[p]);
						      DIPDv=DIPDv+(Delta_IPDv[p]);
							  DIPDx=DIPDx+(Delta_IPDx[p]);
						 }
						 DSNR=DSNR/BINS;
						 DILDx=DILDx/BINS;
						 DILDv=DILDv/BINS;
						 DIPDx=DIPDx/BINS;
						 DIPDv=DIPDv/BINS;

						 GLOBAL_MEASURES << std::fixed << std::setprecision(16) << DSNR << "\t" << (DILDx) << "\t" << (DILDv) << "\t" << (DIPDx) << "\t" << (DIPDv) << "\n";
					}

                    #endif

					for (int i = 0; i < BINS; i++)
						for (int k = 0; k < 2; k++)
						{
							frame_out_tmp_in[i][k] =frame_out_tmp_in[i][k] + frame_out_tmp_in_aux[i][k]; //Frame de saída no tempo apenas para INPUT
						}

					int x = 0;
					for (int outCount = 0; outCount < BLOCK; outCount++)
					{
						for (int j = 0; j < 2; j++)
						{
							output_signal << std::fixed << std::setprecision(16) << frame_out_tmp_in[outCount][j].real();
							if (j==0)
								output_signal << "\t";
						}
						output_signal << "\n";
						x++;
					}

				}
			}
			//wSaida.close();
			//Rv_Matrix.close();
		    output_signal.close();
		    //frame_entrada_tempo.close();
			inputAudio->close();
			inputAudioVAD->close();
			std::cout << "GAMA =" << gama[countGama] << "\n" << "BETA =" << beta << "\n" << "LAMBDA =" << lambda << "\n";
            #ifdef METRICAS
			//deltaSNR.close();
			//deltaILD.close();
			debugFFT.close();
			GLOBAL_MEASURES.close();
			power_vectors.close();
            #endif
		}

	}

	/*************************************************************************************************************************/
	/**************************************************  Metrics Calculation *************************************************/
	/*************************************************************************************************************************/
/*
#ifdef METRICAS

	// Time Variables (Colocar essas variáveis no .hpp)
	float input[6];                                                 //Colocar essas variáveis no .hpp
	//float frame_in_tmp_inFala[256][6];                              //Colocar essas variáveis no .hpp
	//float frame_in_tmp_inRuido[256][6];                             //Colocar essas variáveis no .hpp

	// Frequency Variables
	std::complex<float> frame_out_frq_out_puroFala[256][2];         //Colocar essas variáveis no .hpp   - frame_in_frq_sp
	std::complex<float> frame_out_frq_out_filtradoFala[256][2];     //Colocar essas variáveis no .hpp   - frame_out_frq_sp
	std::complex<float> frame_out_frq_out_puroRuido[256][2];        //Colocar essas variáveis no .hpp   - frame_in_frq_no
	std::complex<float> frame_out_frq_out_filtradoRuido[256][2];    //Colocar essas variáveis no .hpp   - frame_out_frq_no

	// Coefficients
	std::complex<float> coefsSaida[12][256];                        //Colocar essas variáveis no .hpp

	//Original Input Speech and Noise audio files
	std::ifstream inputSpeechAudio(Input_Speech_audioFile);
    std::ifstream inputNoiseAudio(Input_Noise_audioFile);

    //SNR Data
    float snrPuroL = 0;
    float snrPuroR = 0;
    float snrFiltradoL = 0;
    float snrFiltradoR = 0;

    //ILD
    float ild_sp = 0;
    float ild_no = 0;
    float Delta_ild_sp = 0;
    float Delta_ild_no = 0;




    int duplexFreq = floor(1500*BINS/Experiment.FAM)+1;


///// Code starts here!

    for (int numTeste = 0; numTeste < number_of_experiments; numTeste++){
		for (int countGama = 0; countGama < sizeGama; countGama++){

			 //Coefficients
			 std::ifstream coefs(Output_Coefficients[countGama]);
             //Output Files
			 std::ofstream inputSNR(Input_SNR[countGama]);
			 std::ofstream outputSNR(Output_SNR[countGama]);
			 std::ofstream deltaSNR(Delta_SNR[countGama]);
			 std::ofstream deltaILD(Delta_ILD[countGama]);
			 std::ofstream deltaIPD(Delta_IPD[countGama]);

//		if (numTeste == 0)
//				inRuido.open("ruidoF+60.txt");
//			else
//				inRuido.open("ruidoF-60.txt");


			// IMPORTANTES

			for (int i=0; i<BINS;i++)
				for (int j=0; j<NUM_MICS;j++)
				{
					frame_in_tmp_sp[i][j] = 0;
					frame_in_tmp_no[i][j] = 0;
				}

			//int ct_frame = 0;
			cont_frames=0;

			float aux_sp = 0;
			float aux_no = 0;

			float vec_in = 0;
			float vec_out = 0;

			float ILD_in_LF = 0;
			float ILD_out_LF = 0;



			float itd_sp[256];
			float itd_no[256];
			std::complex<float> itd_sp_o[256];
			std::complex<float> itd_no_o[256];
			std::complex<float> itd_sp_i[256];
			std::complex<float> itd_no_i[256];
			float ild_sp_vec[256];
			float ild_no_vec[256];

			float temp_itd_sp[WINDOW];
			float temp_itd_no[WINDOW];
			float Delta_itd_sp = 0;
			float Delta_itd_no = 0;

			float powerFreqPuroFala[2] ={0.0,0.0};
			float powerFreqFiltradoFala[2] = {0.0,0.0};
			float powerFreqPuroRuido[2] ={0.0,0.0};
			float powerFreqFiltradoRuido[2] = {0.0,0.0};



/////////////////////////////////   Load the coefficients   ////////////////////////////////////////////////////////


			   //Começa processamento

			for (cont_samples = 0; cont_samples < tamSign; cont_samples++)
			{

///////////////////////////////////   Speech       ///////////////////////////////////////////////////

				for (int k =0; k < BINS-1; k++)  // montando frame de Bins amostras
					for (int g =0; g < NUM_MICS; g++) // NUM_MICS
						frame_in_tmp_sp[k][g] = frame_in_tmp_sp[k+1][g];   // Monta frame de speech

				std::getline(inputSpeechAudio, line); //Load data from Speech File

				std::istringstream iss(line);
				for (int i = 0; i < NUM_MICS; i++) {
					std::getline(iss, token, '\t');
					float numberDub = sciToDub(token);
					input[i] = numberDub;
				}

				for (int g = 0; g < NUM_MICS; g++)
					frame_in_tmp_sp[BINS-1][g] = input[g];


/////////////////////////////////   Noise   ////////////////////////////////////////////////////////

				for (int k =0; k < BINS-1; k++)
					for (int g =0; g < NUM_MICS; g++)
						frame_in_tmp_no[k][g] = frame_in_tmp_no[k+1][g];

				std::getline(inputNoiseAudio, line);  //Load data from Noise File

				std::istringstream iss2(line);
				for (int i = 0; i < NUM_MICS; i++) {
					std::getline(iss2, token, '\t');
					float numberDub = sciToDub(token);
					input[i] = numberDub;
				}

				for (int g = 0; g < NUM_MICS; g++)
					frame_in_tmp_no[BINS-1][g] = input[g];


/////////////////////////////////   Montagem do frame	////////////////////////////////////////////////

				if (((cont_samples+1) % BLOCK) == 0)
				{
					cont_frames++;  //Increment Frame counter


					float realPart, imagPart; //Já foi declaradolá em cima

					//Montagem dos coeficientes complexos
					for (int j = 0; j < BINS; j++){

							std::getline(coefs, line);
							std::istringstream iss(line);
							for (int i = 0; i < 2*NUM_MICS; i++){
								std::getline(iss, token, ' ');
								realPart = sciToDub(token);
								std::getline(iss, token, '\t');
								imagPart = sciToDub(token);
								std::complex<float> entry(realPart, imagPart);
								coefsSaida[i][j] = entry;
							}
					}

					float somaFala = 0;
					float somaFalaRuido = 0;
					for (int k = 0; k < BINS; k++)
					{
						for (int y = 0; y < NUM_MICS; y++)
						{
							somaFala = somaFala + std::abs(frame_in_tmp_sp[k][y]);  //potência da fala
							somaFalaRuido = somaFalaRuido + std::abs(frame_in_tmp_no[k][y]); //potência ruído
						}
					}

					if (somaFala > 0.0)
					{
						mainFilterForResults(cont_frames,frameEndPriorEst, wm, weighting, frame_in_tmp_sp, qL, qR, frame_out_frq_out_puroFala, frame_out_frq_out_filtradoFala, coefsSaida);

						if (somaFalaRuido > 0.0)
							mainFilterForResults(cont_frames,frameEndPriorEst, wm, weighting, frame_in_tmp_no, qL, qR, frame_out_frq_out_puroRuido, frame_out_frq_out_filtradoRuido, coefsSaida);


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
						// Delta SNR_L - Delta SNR_R    0 -> Left Side  1 - Right Side

						for (int lado = 0; lado < 2; lado++)
						{
							float somaBinPuroFala = 0;
							float somaBinFiltradoFala = 0;
							float somaBinPuroRuido = 0;
							float somaBinFiltradoRuido = 0;

							for (int bin = 0; bin < 256; bin++)
							{

								somaBinPuroFala = somaBinPuroFala + pow(std::abs(frame_out_frq_out_puroFala[bin][lado]),2);
								somaBinFiltradoFala = somaBinFiltradoFala + pow(std::abs(frame_out_frq_out_filtradoFala[bin][lado]),2);
								somaBinPuroRuido = somaBinPuroRuido + pow(std::abs(frame_out_frq_out_puroRuido[bin][lado]),2);
								somaBinFiltradoRuido = somaBinFiltradoRuido + pow(std::abs(frame_out_frq_out_filtradoRuido[bin][lado]),2);


							}
							powerFreqPuroFala[lado] = somaBinPuroFala;
							powerFreqFiltradoFala[lado] = somaBinFiltradoFala;
							powerFreqPuroRuido[lado] = somaBinPuroRuido;
							powerFreqFiltradoRuido[lado] = somaBinFiltradoRuido;
						}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					// IPD

						for (int bin = 0; bin < BINS; bin++)
						{
							itd_sp[bin] =  pow(std::arg(frame_out_frq_out_filtradoFala[bin][0] * std::conj(frame_out_frq_out_filtradoFala[bin][1]) * frame_out_frq_out_puroFala[bin][1] * std::conj(frame_out_frq_out_puroFala[bin][0])),2);
							itd_no[bin] =  pow(std::arg(frame_out_frq_out_filtradoRuido[bin][0] * std::conj(frame_out_frq_out_filtradoRuido[bin][1]) * frame_out_frq_out_puroRuido[bin][1] * std::conj(frame_out_frq_out_puroRuido[bin][0])),2);

							//Aqui correto
							itd_sp_o[bin] = frame_out_frq_out_filtradoFala[bin][0] * std::conj(frame_out_frq_out_filtradoFala[bin][1]);
							itd_sp_i[bin] = frame_out_frq_out_puroFala[bin][0] * std::conj(frame_out_frq_out_puroFala[bin][1]);
							itd_no_i[bin] = frame_out_frq_out_puroRuido[bin][0] * std::conj(frame_out_frq_out_puroRuido[bin][1]);
							itd_no_o[bin] = frame_out_frq_out_filtradoRuido[bin][0] * std::conj(frame_out_frq_out_filtradoRuido[bin][1]);

							//ild_sp_vec[bin] = pow( (log10((abs(frame_out_frq_out_filtradoFala[bin][0] * std::conj(frame_out_frq_out_puroFala[bin][1])))  / (abs( frame_out_frq_out_filtradoFala[bin][1] * std::conj(frame_out_frq_out_puroFala[bin][0]))) ) ) , 2);
							//ild_no_vec[bin] = pow( (log10((abs(frame_out_frq_out_filtradoRuido[bin][0] * std::conj(frame_out_frq_out_puroRuido[bin][1])))  / (abs( frame_out_frq_out_filtradoRuido[bin][1] * std::conj(frame_out_frq_out_puroRuido[bin][0]))) ) ) , 2);

	//						ild_sp(bin) = ild_sp(bin) + ( log10(   abs(    out_frq(bin,3) * conj(out_frq(bin,5))    ) / abs(   out_frq(bin,7) * conj(out_frq(bin,1))  )    ) )^2;

	//						itd_sp(bin) = itd_sp(bin) + ( phase(    out_frq(bin,3) * conj(out_frq(bin,7)) * out_frq(bin,5) * conj(out_frq(bin,1))        ) )^2;

						}



						float soma_itd_no = 0;
						float soma_itd_sp = 0;
						for (int bin = 0; bin < WINDOW; bin++)
						{
							//ex. temp_itd_sp[bin] = abs(std::arg(Rx[1][4][bins])-std::arg(Rx[1][4][bins]));
							temp_itd_sp[bin] = abs(std::arg(itd_sp_o[bin]) - std::arg(itd_sp_i[bin]));
							temp_itd_no[bin] = abs(std::arg(itd_no_o[bin]) - std::arg(itd_no_i[bin]));

							if (temp_itd_no[bin] > M_PI)
							{
								temp_itd_no[bin] = abs((2*M_PI)-temp_itd_no[bin]);
							}
							if (temp_itd_sp[bin] > M_PI)
							{
								temp_itd_sp[bin] = abs((2*M_PI)-temp_itd_sp[bin]);
							}

							soma_itd_no = soma_itd_no + temp_itd_no[bin];
							soma_itd_sp = soma_itd_sp + temp_itd_sp[bin];
						}


						Delta_itd_sp= (soma_itd_sp/WINDOW)/M_PI;
						Delta_itd_no= (soma_itd_no/WINDOW)/M_PI;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

						if (powerFreqPuroFala[0] > 0.00001)
						{
							snrPuroL = 10 * log10(powerFreqPuroFala[0]/powerFreqPuroRuido[0]);
							snrPuroR = 10 * log10(powerFreqPuroFala[1]/powerFreqPuroRuido[1]);
							snrFiltradoL = 10 * log10(powerFreqFiltradoFala[0]/powerFreqFiltradoRuido[0]);
							snrFiltradoR = 10 * log10(powerFreqFiltradoFala[1]/powerFreqFiltradoRuido[1]);

							aux_sp = (powerFreqFiltradoFala[0]*powerFreqPuroFala[1]) / (powerFreqPuroFala[0] * powerFreqFiltradoFala[1]);
							aux_no = (powerFreqFiltradoRuido[0]*powerFreqPuroRuido[1]) / (powerFreqPuroRuido[0] * powerFreqFiltradoRuido[1]);;

						//	vec_in = powerFreqPuroRuido[0] / powerFreqPuroRuido[1];
						//	vec_out = powerFreqFiltradoRuido[0] / powerFreqFiltradoRuido[1];

						//	ILD_in_LF = 10*log10(vec_in); // ILD de entrada, ruido esquerdo / ruido direito
						//	ILD_out_LF = 10*log10(vec_out); // ILD de saida, ruido filtrado esquerdo / ruido filtrado direito

							//Esse aqui é o delta ILD
							ild_sp = abs(10*log10(aux_sp)); // ILD fala, ?
							ild_no = abs(10*log10(aux_no)); // ILD ruido, ?

							//Delta_ild_sp = 0.33333333 * ild_sp; // Delta ILD fala, ?
							//Delta_ild_no = 0.33333333 * ild_no; // Delta ILD ruido, ?
						}
						else
						{
							snrPuroL = -10;
							snrPuroR = -10;
							snrFiltradoL = -10;
							snrFiltradoR = -10;
						}

						//std::cout << "PURO= SNR L: " << snrPuroL << " SNR R: " << snrPuroR <<  "FILTRADO= SNR L: " << snrFiltradoL << " SNR R: " << snrFiltradoR << "\n";
					}
					else
					{
						snrPuroL = -10;
						snrPuroR = -10;
						snrFiltradoL = -10;
						snrFiltradoR = -10;

						//std::cout << "PURO= SNR L: " << snrPuroL << " SNR R: " << snrPuroR <<  "FILTRADO= SNR L: " << snrFiltradoL << " SNR R: " << snrFiltradoR << "\n";
					}

					inputSNR << std::fixed << std::setprecision(16) << snrPuroL << "\t" << snrPuroR << "\n";
					outputSNR << std::fixed << std::setprecision(16) << snrFiltradoL << "\t" << snrFiltradoR << "\n";
					//Fábio
					deltaSNR << std::fixed << std::setprecision(16) << (snrFiltradoL-snrPuroL) << "\t" << (snrFiltradoR-snrPuroR) << "\n";
					//ild << std::fixed << std::setprecision(16) << ild_sp << "\t" << ild_no << "\n";
					//deltaILD << std::fixed << std::setprecision(16) << Delta_ild_sp << "\t" << Delta_ild_no << "\n";
					deltaILD << std::fixed << std::setprecision(16) << ild_sp << "\t" << ild_no << "\n";
					deltaIPD << std::fixed << std::setprecision(16) << Delta_itd_sp << "\t" << Delta_itd_no << "\n";

					//Novas variáveis de métricas
					//DILDx << std::fixed << std::setprecision(16) << ild_sp << "\n";
					//DILDv << std::fixed << std::setprecision(16) << ild_no << "\n";
				}
			}


	// Fechamento dos arquivos
			inputSNR.close();
			outputSNR.close();
			deltaSNR.close();
			deltaILD.close();
			deltaIPD.close();

		}

    }
#endif
*/
 return 0;
}

