import ctypes
from ctypes import POINTER, c_double 

clibrary = ctypes.CDLL("src/cpython.so")

clibrary.calculate_transformation_matrix.argtypes = [POINTER(c_double), POINTER(c_double)]
clibrary.calculate_transformation_matrix.restype = POINTER(c_double)

state_1 = (c_double * 6)(1.0, 0.0, 0.0, 4.0, 5.0, 6.0)
state_2 = (c_double * 6)(1.0, 2.0, 3.0, 4.0, 5.0, 6.0)
result_ptr = clibrary.calculate_transformation_matrix(state_1, state_2)
result = [result_ptr[i] for i in range(16)]

print(result)
