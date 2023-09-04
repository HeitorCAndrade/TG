from fxpmath import Fxp

n_word = 64
n_channel = 6
nfft = 256

zero_fp = Fxp(0, n_word=n_word*n_channel)
line = zero_fp.hex()
line = line[2:]

with open('zero_init_mem.txt', 'w') as f:
    for i in range(nfft):
        f.write(line)
        if (i<nfft-1):
          f.write('\n')