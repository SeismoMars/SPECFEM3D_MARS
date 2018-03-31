#!/bin/tcsh
cd /tigress/bozdag/MARS/SSR_R1/specfem3d_marsglobe
make clean
make -j 8 xmeshfem3D
make -j 8 xcreate_header_file
./bin/xcreate_header_file
make -j 8 xspecfem3D
make -j 8 xcombine_vol_data
