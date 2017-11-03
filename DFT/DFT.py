# For this part of the assignment, please implement your own code for all computations,
# Do not use inbuilt functions like fft from either numpy, opencv or other libraries
import numpy as np
import math

class DFT:

    def forward_transform(self, matrix):
        """Computes the forward Fourier transform of the input matrix
        takes as input:
        matrix: a 2d matrix
        returns a complex matrix representing fourier transform"""

        N = matrix.shape[0]
        result_mat = np.zeros((N,N),dtype=np.complex_)

        #print(result_mat)
        #u,v = result_mat.shape
        #print(u,v)

        for u in range(N):
            for v in range(N):
                for i in range(N):
                    for j in range(N):
                        result_mat[u][v] += matrix[i][j]*((math.cos((2*math.pi/N)*(u*i+v*j)))-1j*(math.sin((2*math.pi/N)*(u*i+v*j))))
        #print(np.fft.fft2(matrix))
        return result_mat

    def inverse_transform(self, matrix):
        """Computes the inverse Fourier transform of the input matrix
        matrix: a 2d matrix (DFT) usually complex
        takes as input:
        returns a complex matrix representing the inverse fourier transform"""
        N = matrix.shape[0]
        result_mat = np.zeros((N, N), dtype=np.complex_)
        print("This is input of Inverse")
        print(matrix)
        for i in range(result_mat.shape[0]):
            for j in range(result_mat.shape[1]):
                for u in range(N):
                    for v in range(N):
                        result_mat[i][j] += matrix[u][v]*((math.cos((2*math.pi/N)*(u*i+v*j)))+1j*(math.sin((2*math.pi/N)*(u*i+v*j))))
        #print(np.fft.ifft2(matrix))
        return result_mat


    def discrete_cosine_tranform(self, matrix):
        """Computes the discrete cosine transform of the input matrix
        takes as input:
        matrix: a 2d matrix
        returns a matrix representing discrete cosine transform"""


        N = matrix.shape[0]
        result_mat = np.zeros((N, N))

        # print(result_mat)
        # u,v = result_mat.shape
        # print(u,v)

        for u in range(result_mat.shape[0]):
            for v in range(result_mat.shape[1]):
                for i in range(N):
                    for j in range(N):
                        result_mat[u][v] += matrix[i][j] * (math.cos((2 * math.pi / N) * (u * i + v * j)))
        return result_mat




    def magnitude(self, matrix):
        """Computes the magnitude of the DFT
        takes as input:
        matrix: a 2d matrix
        returns a matrix representing magnitude of the dft"""
        N = matrix.shape[0]
        result_mat = np.zeros((N, N))

        for i in range(0, np.shape(matrix)[0]):
            for j in range(0, np.shape(matrix)[1]):
                result_mat[i][j] = math.sqrt(matrix[i][j].real * matrix[i][j].real + matrix[i][j].imag * matrix[i][j].imag)

        return result_mat