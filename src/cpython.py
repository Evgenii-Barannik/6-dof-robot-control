import ctypes
from ctypes import POINTER, c_double, c_size_t

clibrary = ctypes.CDLL("cpython.so")

clibrary.set_frames_py.argtypes = [
    POINTER(c_double),  # table_x
    POINTER(c_double),  # table_y
    POINTER(c_double),  # table_z
    POINTER(c_double),  # slider_x
    POINTER(c_double),  # slider_y
    POINTER(c_double),  # slider_z
    POINTER(c_double),  # state_1
    POINTER(c_double),  # state_2
    c_size_t            # num_of_frames
]
clibrary.set_frames_py.restype = ctypes.c_bool

state_1 = (c_double * 6)(0.0, 0.0, 0.0, -3.1415/10, 3.1415/10, 3.1415/10)
state_2 = (c_double * 6)(0.5, 0.5, 3.0, -3.1415/10, 3.1415/10, 3.1415/10)

num_of_frames = 4

table_x = (c_double * (num_of_frames * 8))()  
table_y = (c_double * (num_of_frames * 8))()
table_z = (c_double * (num_of_frames * 8))()

slider_x = (c_double * (num_of_frames * 6))()  
slider_y = (c_double * (num_of_frames * 6))()
slider_z = (c_double * (num_of_frames * 6))()

success = clibrary.set_frames_py(
    table_x, table_y, table_z,
    slider_x, slider_y, slider_z,
    state_1, state_2,
    num_of_frames
)

print("==== Python output")
if success:
    print("Table Coordinates:")
    for i in range(num_of_frames):
        print(f"Frame {i}:")
        for j in range(8):  # 8 vectors per frame
            print(f"  Point {j}: ({table_x[i * 8 + j]}, {table_y[i * 8 + j]}, {table_z[i * 8 + j]})")

    print("\nSlider Coordinates:")
    for i in range(num_of_frames):
        print(f"Frame {i}:")
        for j in range(6):  # 6 sliders per frame
            x = slider_x[i * 6 + j]
            y = slider_y[i * 6 + j]
            z = slider_z[i * 6 + j]
            if not (x == y == z == float('nan')):
                print(f"  Slider {j}: ({x}, {y}, {z})")
            else:
                print(f"  Slider {j}: NAN")
else:
    print("Failed to calculate frames.")

