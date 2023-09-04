//----- Defines ------
#define NUM_MICS 6
#define BINS 256
#define WINDOW 128
#define BLOCK 64
#define METRICAS
#define NS 76
//----- Metrics Variables ------
float DILDv=0;
float DIPDv=0;
float DILDx=0;
float DIPDx=0;
float DSNR=0;

//----- Variables ------
int tamSign = 160193;
int frameEndPriorEst =1001;
int total_number_of_frames =2504;
int frames_metricas =1503;
float gama[1] = {0};
int sizeGama=1;
float beta = 0.001;
float lambda = 0.98684;
int number_of_experiments= 1;
float weighting[256];
float qL[6]={1,0,0,0,0,0};
float qR[6]={0,0,0,1,0,0};
std::complex<float> Rv[6][6][256],Ry[6][6][256],Rx[6][6][256];
//----- Counters ------
int cont_samples = 0;
int cont_frames = 0;
int WM = 100;
int cont_v = 0;
int cont_y = 0;
//----- Buffers ------
int frame_vad[256];
float frame_in_tmp_in[256][6];
float frame_in_tmp_sp[256][6];
float frame_in_tmp_no[256][6];
std::complex<float> frame_in_frq_in[256][6];
std::complex<float> frame_in_frq_sp[256][6];
std::complex<float> frame_in_frq_no[256][6];
std::complex<float> w[12][256];
std::complex<float> frame_out_tmp_in[256][2];
std::complex<float> frame_out_frq_in[256][2];
std::complex<float> frame_out_tmp_in_aux[256][2];
std::complex<float> frame_out_tmp_sp[256][2];
std::complex<float> frame_out_frq_sp[256][2];
std::complex<float> frame_out_tmp_no[256][2];
std::complex<float> frame_out_frq_no[256][2];
//----- Struct ------
struct Data_Experiments{
std::string SPEECH = {"FA03_09.wav"};
std::string NOISE = {"ICRA_No01_16kHz_12s.wav"};
std::string ENVIRONMENT = {"Anechoic"};
int AZS = 0;
int AZN = -60;
int SNR = 0;
int FAM = 16000;
int REPETITION = 2;
int ESTIMATION_TIME = 4;};
//----- Input Files ------
std::string Input_Speech_audioFile ={"Input_SPEECH_FA03_09.wav_AZS_0_NOISE_ICRA_No01_16kHz_12s.wav_AZN_-60_SNR_0_dB_MICS_6_PRIOR_4s_REPETITION_2.txt"};
std::string Input_Noise_audioFile ={"Input_NOISE_FA03_09.wav_AZS_0_NOISE_ICRA_No01_16kHz_12s.wav_AZN_-60_SNR_0_dB_MICS_6_PRIOR_4s_REPETITION_2.txt"};
//----- TestBench Files ------
std::string testBench_audioFile ={"TestBench_SPEECH_FA03_09.wav_AZS_0_NOISE_ICRA_No01_16kHz_12s.wav_AZN_-60_SNR_0_dB_MICS_6_PRIOR_4s_REPETITION_2.txt"};
std::string testBench_VAD ={"VAD_TestBench_FA03_09.wav_AZS_0_NOISE_ICRA_No01_16kHz_12s.wav_AZN_-60_SNR_0_dB_MICS_6_PRIOR_4s_REPETITION_2.txt"};
//----- Output Files ------

std::string Output_File_Name[1]={"Output_File_FA03_09.wav_AZS_0_NOISE_ICRA_No01_16kHz_12s.wav_AZN_-60_BETA_0.001_GAMA_0_SNR_0_dB_MICS_6_PRIOR_4s_REPETITION_2.txt"};
std::string Output_Coefficients[1]={"Output_Coefficients_FA03_09.wav_AZS_0_NOISE_ICRA_No01_16kHz_12s.wav_AZN_-60_BETA_0.001_GAMA_0_SNR_0_dB_MICS_6_PRIOR_4s_REPETITION_2.txt"};
//----- Metric Files ------
std::string MEASURES[1]={"MEASURES_FA03_09.wav_AZS_0_NOISE_ICRA_No01_16kHz_12s.wav_AZN_-60_BETA_0.001_GAMA_0_SNR_0_dB_MICS_6_PRIOR_4s_REPETITION_2.txt"};
std::string debug_FFT[1]={"debug_FFT_FA03_09.wav_AZS_0_NOISE_ICRA_No01_16kHz_12s.wav_AZN_-60_BETA_0.001_GAMA_0_SNR_0_dB_MICS_6_PRIOR_4s_REPETITION_2.txt"};
