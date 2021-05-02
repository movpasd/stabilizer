import numpy as np

CX = np.array([[1, 0, 0, 0],
               [0, 1, 0, 0],
               [0, 0, 0, 1],
               [0, 0, 1, 0]], dtype=complex)

CZ = np.array([[1, 0, 0, 0],
               [0, 1, 0, 0],
               [0, 0, 1, 0],
               [0, 0, 0, -1]], dtype=complex)

X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)

X_X = np.kron(X, X)
X_Z = np.kron(X, Z)
Z_X = np.kron(Z, X)
Z_Z = np.kron(Z, Z)

ket_plus = np.array([1, 1], dtype=complex)


x = CX @ (np.kron(ket_plus, ket_plus))


print(np.real(x))
print(np.imag(x))