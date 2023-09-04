from fxpmath import Fxp

n_words = 32
n_frac = 20
first_time = 1

fixed_point_0 = []
fixed_point_1 = []
fixed_point_2 = []
fixed_point_3 = []
fixed_point_4 = []
fixed_point_5 = []
n = 0
zero_word = Fxp(0, n_word=n_words)
zero_word_str = zero_word.hex()
zero_word_str = zero_word_str[2:]
if (first_time == 1):
    with open('input_signal.txt') as f:
        for line in f:
            n = n+1
            float_str = line.split(',')
            float1 = float(float_str[0])
            float2 = float(float_str[1])
            float3 = float(float_str[2])
            float4 = float(float_str[3])
            float5 = float(float_str[4])
            float6 = float(float_str[5])
            fixed_point_0.append(Fxp(float1, n_word=n_words, n_frac=n_frac))
            fixed_point_1.append(Fxp(float2, n_word=n_words, n_frac=n_frac))
            fixed_point_2.append(Fxp(float3, n_word=n_words, n_frac=n_frac))
            fixed_point_3.append(Fxp(float4, n_word=n_words, n_frac=n_frac))
            fixed_point_4.append(Fxp(float5, n_word=n_words, n_frac=n_frac))
            fixed_point_5.append(Fxp(float6, n_word=n_words, n_frac=n_frac))



    with open('audio_input_fp_tb.txt', 'w') as f:
        for i in range(n):
            fp0_str = fixed_point_0[i].hex()
            fp0_str = fp0_str[2:]
            fp1_str = fixed_point_1[i].hex()
            fp1_str = fp1_str[2:]
            
            fp2_str = fixed_point_2[i].hex()
            fp2_str = fp2_str[2:]
            fp3_str = fixed_point_3[i].hex()
            fp3_str = fp3_str[2:]

            fp4_str = fixed_point_4[i].hex()
            fp4_str = fp4_str[2:]
            fp5_str = fixed_point_5[i].hex()
            fp5_str = fp5_str[2:]
            line = zero_word_str+fp5_str +  zero_word_str+fp4_str+  zero_word_str+fp3_str+ zero_word_str+ fp2_str+ zero_word_str+ fp1_str +zero_word_str+ fp0_str
            f.write(line)
            f.write('\n')
            f.flush()

""" lines = []
channels = []
with open('audio_input_fp_tb.txt') as f:
    for line in f:
        channels = line.split(' ')
        new_str = ''
        for str_i in channels:
            new_str = new_str + str_i[2:]
        lines.append(new_str)

with open('audio_input_fp_tb_fixed.txt', 'w') as f:
    for line in lines:
        f.write(line)
        f.write('\n') """
