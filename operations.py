import numpy as np
from stabilizers import *

HI = 1 / np.sqrt(2) * np.array([[1, 0, -1, 0],
                            [0, 1, 0, -1],
                            [-1, 0, 1, 0],
                            [0, -1, 0, 1]])

CX = np.array([[1, 0, 0, 0],
               [0, 1, 0, 0],
               [0, 0, 0, 1],
               [0, 0, 1, 0]])

print(HI @ CX @ HI)


