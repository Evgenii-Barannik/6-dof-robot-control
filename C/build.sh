conan install . --output-folder=build --build=missing
cd build
meson setup --native-file conan_meson_native.ini .. meson-src

meson compile -C meson-src
if [ $? -ne 0 ]; then
    echo "Meson compilation failed."
    exit 1
fi

gcc -fPIC -shared -o ../src/cpython.so ../src/main.c -lgsl -lgslcblas
if [ $? -ne 0 ]; then
    echo "C-Python compilation failed."
    exit 1
fi

./meson-src/6dof
