#!/bin/bash
# Build Vina-GPU 2.1 for AMD RDNA4 (gfx1201) with MinGW-w64 (run from MSYS2).
# Produces a standalone Vina-GPU-AMD.exe (DLLs copied alongside).
# Prereqs: pacman -S mingw-w64-x86_64-{gcc,boost,opencl-headers,opencl-icd}
# Apply the patches in patches.md to a fresh clone first.
set -e
export PATH=/mingw64/bin:/usr/bin:$PATH
DIR="${1:-$HOME/gpu_docking/Vina-GPU-2.1/AutoDock-Vina-GPU-2.1}"
cd "$DIR"
g++ -o Vina-GPU-AMD.exe \
  -I./lib -I./OpenCL/inc -I/mingw64/include -I/mingw64/include/boost \
  ./main/main.cpp -O3 \
  ./lib/*.cpp ./OpenCL/src/wrapcl.cpp \
  -L/mingw64/lib \
  -lboost_program_options-mt -lboost_filesystem-mt -lboost_thread-mt -lOpenCL \
  -lstdc++ -lstdc++fs -lm -lpthread \
  -DOPENCL_2_0 -DAMD_PLATFORM -DSMALL_BOX -DBOOST_TIMER_ENABLE_DEPRECATED \
  -DCL_TARGET_OPENCL_VERSION=200 -DNDEBUG -DBUILD_KERNEL_FROM_SOURCE
# copy non-system runtime DLLs next to the exe -> standalone (callable from WSL/cmd)
ldd Vina-GPU-AMD.exe | grep -i mingw64 | awk '{print $3}' | while read d; do
  [ -f "$d" ] && cp -u "$d" .
done
echo "built $(stat -c%s Vina-GPU-AMD.exe) bytes + DLLs"
