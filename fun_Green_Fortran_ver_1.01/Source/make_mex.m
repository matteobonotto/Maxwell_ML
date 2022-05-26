clc; close all; clearvars;
clear functions;

mex -V ...
    COMPFLAGS="$COMPFLAGS -openmp"...
    -largeArrayDims ...
    -output ...
    fun_Green_filament_Aphi_SP_f90 ...
    module_Green_axi_MB.f90 ...
    fun_Green_filament_Aphi_SP_4mex.f90

mex -V ...
    COMPFLAGS="$COMPFLAGS -openmp"...
    -largeArrayDims ...
    -output ...
    fun_Green_filament_flux_SP_f90 ...
    module_Green_axi_MB.f90 ...
    fun_Green_filament_flux_SP_4mex.f90

mex -V ...
    COMPFLAGS="$COMPFLAGS -openmp"...
    -largeArrayDims ...
    -output ...
    fun_Green_filament_BrBz_SP_f90 ...
    module_Green_axi_MB.f90 ...
    fun_Green_filament_BrBz_SP_4mex.f90




















