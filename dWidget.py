import numpy as np
import matplotlib.pyplot as plt

# --re^ik+e^-ik------1-j1-A-
#                    |   /
#                    j2 j3
#                    | /
#                    B

# j2 B + j1 A + re^ik = e^ik
# j1^* + j3^* B = Ae^-ik
# j2^* + j3 A = Be^-ik

p = 1/2
lambd = 1 - (2 * p) / 3
j1 = lambd - np.sqrt(1 - lambd**2) * (1 + 1j) / np.sqrt(2)
j2 = lambd * (1 + 1j) / np.sqrt(2) + np.sqrt(1 - lambd**2)
j3 = 1j


