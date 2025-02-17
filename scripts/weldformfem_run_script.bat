#@echo off
call setvars.bat
echo ----------------------------------------------------------------
set WELDFORMFEM_PATH = D:\Luciano\Numerico\WeldFormFEM_bin_cxx\
set PATH=%PATH%;%WELDFORMFEM_PATH%
set file=%1
set solver = \src\WeldFormFEM.exe
echo %PATH%
D:\Luciano\Numerico\WeldFormFEM_bin_cxx\src\WeldFormFEM.exe  %1