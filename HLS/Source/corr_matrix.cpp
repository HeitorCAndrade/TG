#include "mainProject.hpp"

void corr_matrix_estimation(int estimation_mode, int length_fft, int channel, std::complex<float> frame_y[256][NUM_CHAN], std::complex<float> (&rv)[NUM_CHAN][NUM_CHAN][256], std::complex<float> (&ry)[NUM_CHAN][NUM_CHAN][256], std::complex<float> (&rx)[NUM_CHAN][NUM_CHAN][256], int vad, int & cont_y, int & cont_v, float smooth_coef)
{
	int channelNum = 6;
	std::complex<float> frame_y_inv[6];
	std::complex<float> frame_y_conj[6];
	std::complex<float> frame_y_new[6][6];
	

	if (estimation_mode == 0)    //Before frameEndPriorEst
	{
		if(vad == 0) {  //noise only frame (vad=sumVad)
			for(int bin=0; bin<length_fft; bin++)
			{
				for(int j=0; j<channelNum; j++)
				{
					frame_y_inv[j] = frame_y[bin][j];
					frame_y_conj[j] = std::conj(frame_y[bin][j]);
				}				
				
				for(int i = 0; i < channelNum; ++i)
					for(int j = 0; j < channelNum; ++j)
					{
						std::complex<float> number = frame_y_inv[i] * frame_y_conj[j];
						frame_y_new[i][j] = number;
					}
					
				for(int i = 0; i < channelNum; ++i)
					for(int j = 0; j < channelNum; ++j)
						rv[i][j][bin] = rv[i][j][bin] + frame_y_new[i][j];
			}
			cont_v++;
		}
		else if(vad > 0) {   //speech + Noise
			for(int bin=0; bin<length_fft; bin++)
			{
				for(int j=0; j<channelNum; j++)
				{
					frame_y_inv[j] = frame_y[bin][j];
					frame_y_conj[j] = std::conj(frame_y[bin][j]);
				}				
				
				for(int i = 0; i < channelNum; ++i)
					for(int j = 0; j < channelNum; ++j)
					{
						std::complex<float> number = frame_y_inv[i] * frame_y_conj[j];
						frame_y_new[i][j] = number;
					}
					
				for(int i = 0; i < channelNum; ++i)
					for(int j = 0; j < channelNum; ++j)
						ry[i][j][bin] = ry[i][j][bin] + frame_y_new[i][j];
			}
			cont_y++;
		}
		else
			return; // error
		//std::cout << "DEU RUIM! =" << "\n";
	}
	else if (estimation_mode == 1) {  //After frameEndPriorEst
		if(vad == 0) {
			for(int bin=0; bin<length_fft; bin++)
			{
				for(int j=0; j<channelNum; j++)
				{
					frame_y_inv[j] = frame_y[bin][j];
					frame_y_conj[j] = std::conj(frame_y[bin][j]);
				}				
				
				for(int i = 0; i < channelNum; ++i)
					for(int j = 0; j < channelNum; ++j)
					{
						std::complex<float> number = (1-smooth_coef) * frame_y_inv[i] * frame_y_conj[j];
						frame_y_new[i][j] = number;
					}
					
				for(int i = 0; i < channelNum; ++i)
					for(int j = 0; j < channelNum; ++j)
						rv[i][j][bin] = (rv[i][j][bin] * smooth_coef) + frame_y_new[i][j];
			}
		}
		else if(vad > 0) {
			for(int bin=0; bin<length_fft; bin++)
			{
				for(int j=0; j<channelNum; j++)
				{
					frame_y_inv[j] = frame_y[bin][j];
					frame_y_conj[j] = std::conj(frame_y[bin][j]);
				}				
				
				for(int i = 0; i < channelNum; ++i)
					for(int j = 0; j < channelNum; ++j)
					{
						std::complex<float> number = (1-smooth_coef) * frame_y_inv[i] * frame_y_conj[j];
						frame_y_new[i][j] = number;
					}
					
				for(int i = 0; i < channelNum; ++i)
					for(int j = 0; j < channelNum; ++j)
						ry[i][j][bin] = (ry[i][j][bin] * smooth_coef) + frame_y_new[i][j];
			}
		}
		else
			return; // error

		for(int bin = 0; bin < channelNum; bin++)
			for(int i = 0; i < channelNum; ++i)
				for(int j = 0; j < 256; ++j)
					rx[bin][i][j] = ry[bin][i][j] - rv[bin][i][j];
	}
	else
		return; //error

	return;
}
