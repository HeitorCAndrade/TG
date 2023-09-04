from fxpmath import Fxp


pattern1_re = []
pattern1_im = []

pattern2_re = []
pattern2_im = []

pattern3_re = []
pattern3_im = []

pattern4_re = []
pattern4_im = []

fp_pt_1_re = []
fp_pt_1_im = []

fp_pt_2_re = []
fp_pt_2_im = []

fp_pt_3_re = []
fp_pt_3_im = []

fp_pt_4_re = []
fp_pt_4_im = []

for i in range(256):
    pattern2_re.append(1)
    pattern2_im.append(0)

    pt2 = i % 4
    pt3 = i % 4
    pt4 = i % 2

    if (pt2 == 0):
        pattern1_re.append(1)
        pattern1_im.append(0)
    
    if (pt2 == 1):
        pattern1_re.append(0)
        pattern1_im.append(-1)
    
    if (pt2 == 2):
        pattern1_re.append(-1)
        pattern1_im.append(0)
    
    if (pt2 == 3):
        pattern1_re.append(0)
        pattern1_im.append(1)

    if (pt3 == 0):
        pattern3_re.append(1)
        pattern3_im.append(0)
    
    if (pt3 == 1):
        pattern3_re.append(0)
        pattern3_im.append(1)

    if (pt3 == 2):
        pattern3_re.append(-1)
        pattern3_im.append(0)

    if (pt3 == 3):
        pattern3_re.append(0)
        pattern3_im.append(-1)
    
    if (pt4 == 0):
        pattern4_re.append(1)
        pattern4_im.append(0)
    
    if (pt4 == 1):
        pattern4_re.append(-1)
        pattern4_im.append(0)


for i in range(256):
    # pattern1_re[i] = pattern1_re[i] * int((-1) ** i)
    # pattern1_im[i] = pattern1_im[i] * int((-1) ** i)

    # pattern2_re[i] = pattern2_re[i] * int((-1) ** i)
    # pattern2_im[i] = pattern2_im[i] * int((-1) ** i)

    # pattern3_re[i] = pattern3_re[i] * int((-1) ** i)
    # pattern3_im[i] = pattern3_im[i] * int((-1) ** i)

    # pattern4_re[i] = pattern4_re[i] * int((-1) ** i)
    # pattern4_im[i] = pattern4_im[i] * int((-1) ** i)

    fp_pt_1_re.append(Fxp(pattern1_re[i], n_word=32, n_frac=20)) 
    fp_pt_1_im.append(Fxp(pattern1_im[i], n_word=32, n_frac=20)) 

    fp_pt_2_re.append(Fxp(pattern2_re[i], n_word=32, n_frac=20)) 
    fp_pt_2_im.append(Fxp(pattern2_im[i], n_word=32, n_frac=20)) 

    fp_pt_3_re.append(Fxp(pattern3_re[i], n_word=32, n_frac=20)) 
    fp_pt_3_im.append(Fxp(pattern3_im[i], n_word=32, n_frac=20)) 

    fp_pt_4_re.append(Fxp(pattern4_re[i], n_word=32, n_frac=20)) 
    fp_pt_4_im.append(Fxp(pattern4_im[i], n_word=32, n_frac=20)) 

    
with open('phase_0.txt', 'w') as f1:
    for i in range(256):
        line_im = fp_pt_1_im[i].hex()
        line_re = fp_pt_1_re[i].hex()
        line_im = line_im[2:]
        line_re = line_re[2:]
        line = line_im + line_re
        f1.write(line)
        f1.write('\n')

with open('phase_1.txt', 'w') as f2:
    for i in range(256):
      line_im = fp_pt_2_im[i].hex()
      line_re = fp_pt_2_re[i].hex()
      line_im = line_im[2:]
      line_re = line_re[2:]
      line = line_im + line_re
      f2.write(line)
      f2.write('\n')

with open('phase_2.txt', 'w') as f3:
    for i in range(256):
      line_im = fp_pt_3_im[i].hex()
      line_re = fp_pt_3_re[i].hex()
      line_im = line_im[2:]
      line_re = line_re[2:]
      line = line_im + line_re
      f3.write(line)
      f3.write('\n')

with open('phase_3.txt', 'w') as f4:
    for i in range(256):
      line_im = fp_pt_4_im[i].hex()
      line_re = fp_pt_4_re[i].hex()
      line_im = line_im[2:]
      line_re = line_re[2:]
      line = line_im + line_re
      f4.write(line)
      f4.write('\n')