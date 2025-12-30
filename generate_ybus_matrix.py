import numpy as np

j = complex(0, 1)
real_min, real_max = -1.9422, 3.3074
imag_min, imag_max = -39.3089, 17.3611

real_random = np.random.uniform(real_min, real_max, (1000, 1000))
imag_random = np.random.uniform(imag_min, imag_max, (1000, 1000))
Ybus = real_random + 1j * imag_random

with open("Ybus_1000x1000_matrix.txt", "w") as f:
    f.write("j = sqrt(-1);\n")
    f.write("Ybus = sparse(1000);\n")
    for i in range(1000):
        for k in range(1000):
            val = Ybus[i, k]
            if abs(val) > 1e-10:
                real_part = f"{val.real:.4f}"
                imag_part = f"{val.imag:.4f}"
                f.write(f"Ybus({i+1:5d},{k+1:5d}) = {real_part}+ j*({imag_part});\n")
