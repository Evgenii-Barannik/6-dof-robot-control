conan install . --output-folder=build --build=missing
cd build
meson setup --native-file conan_meson_native.ini .. meson-src
meson compile -C meson-src
./meson-src/6dof