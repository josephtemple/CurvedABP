@echo off
REM -----------------------------------------------------------
REM ABP C++ build script for Conda Python 3.12 on Windows
REM -----------------------------------------------------------

REM Activate your Conda environment
call conda activate CurvedABP

REM Go to the C++ project folder
cd /d C:\Users\Joseph\Desktop\Assignments\8th Sem\CurvedABP\repo\cpp

REM Remove old build folder and recreate
rmdir /s /q build
mkdir build
cd build

REM Configure CMake with MSVC and Conda Python
cmake -S .. -G "Visual Studio 17 2022" -A x64 -DPython_EXECUTABLE=%CONDA_PREFIX%\python.exe

REM Build the Release version of the module
cmake --build . --config Release

REM Optional: copy the .pyd to your python folder for easy imports
set PYTHON_DIR=C:\Users\Joseph\Desktop\Assignments\8th Sem\CurvedABP\repo\python
copy "Release\abp_core.cp312-win_amd64.pyd" "%PYTHON_DIR%"

echo.
echo Build complete! The module should now be importable in Python:
echo     import abp_core
pause
