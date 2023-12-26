# 6-dof-robot-control
This code shows basic math to control a robot with six vertical parallel rails (6 DOF). The robot here differs from a common Stewart platform and has a rather exotic geometry. The Python and C components in this project are independent, they are realizations of the same general idea. Python code has better visualization (see gif below). To build C code, Conan and Meson are used (see `build.sh`).

![Python visualization][6dof]

[6dof]: https://github.com/Evgenii-Barannik/6-dof-robot-control/blob/main/6dof.gif
