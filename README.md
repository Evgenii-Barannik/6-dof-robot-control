This project contains code for the robotic 3D printer with 6 degrees of freedom (6 DOF). The project has C and Python components:
* C component is where the primary calculations are (see `main.c`). Code in C is written starting from linear algebra; specifically, GSL-BLAS library is used. Main function shows how generated `Trajectory` can be used to interact with hardware. To build C code, Conan and Meson are used (see `build.sh`).
* Python component is used for visualization (see `visualisation.ipynb`). It contains Jupyter notebook built around C functions called using C FFI. This Jupyter notebook requires `cpython.so` to run, which will be produced upon succesful execution of the `build.sh` routine.

Visualization showing transition between two arbitrary robotic printer states (denoted A and B):  
![](robot-6dof-control.gif)

Robotic printer has 6 sliders (6 actuators) which are used to control table position (table is shown as a hexagon in the visualization). Each slider is independent and rides along its own vertical rail; all 6 vertical rails are parallel to each other. The geometry described here is different from "Stewart platform" which also has 6 DOF and is different from more commonly used 3 DOF robotic printer. Roughly speaking this project aims to answer the following question:

Given we know initial (X, Y, Z, roll, yaw, pitch) and desired final state (X, Y, Z, roll, yaw, pitch) of the printer table, how can we move sliders to bring printer table from initial to final state?

Suppose we know the geometry and dimensions of our robotic 3D printer. Suppose we also know that any correct state of the robotic printer can be described using 6 coordinates: X, Y, Z, roll, yaw, pitch of the printer table. To be exact those are the coordinates of the end of the needle growing from the center of the printer table. Futhermore we know that any correct printer table state must have corresponding slider positions. In mechanical sense we will need to change slider positions if we desire to control the table position. Answer to the question above is `Trajectory` struct, which is the main computational result. `Trajectory` struct describes how and when sliders should move to make the transition between states. 
