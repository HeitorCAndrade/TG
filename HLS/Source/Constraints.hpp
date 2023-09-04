int tamSign = 188879; // tamanho sinal

std::string line, token, data;

float gama[2] = {0.0001, 0.0006}; //Simula com 2 valores de Gama  \\valores de gama devem vir do gerador de cenários

std::complex<float> rv[6][6][256], ry[6][6][256], rx[6][6][256]; // desamarrar