@echo off
setlocal EnableExtensions
title CellMorphR
cd /d "%~dp0"

rem ---------------------------------------------------------------------------
rem CellMorphR launcher (Windows).
rem Finds R in this order (first hit wins):
rem   1. R_RSCRIPT env var  (point it at your Rscript.exe to force a choice)
rem   2. Registry: HKLM / HKCU  SOFTWARE\R-core\R  InstallPath
rem   3. Newest install under "C:\Program Files\R\R-*"
rem   4. Rscript on PATH  (system R, or an active conda env with r-base)
rem   5. R inside a conda install (Miniconda3/Anaconda3/miniforge3/mambaforge)
rem Then runs launcher.R, which installs any missing packages on first run and
rem opens the app in your default browser. Keep this window open while working.
rem ---------------------------------------------------------------------------

set "RHOME="
set "RSCRIPT="

if defined R_RSCRIPT if exist "%R_RSCRIPT%" set "RSCRIPT=%R_RSCRIPT%"

if not defined RSCRIPT (
  for /f "skip=2 tokens=2,*" %%A in ('reg query "HKLM\SOFTWARE\R-core\R" /v InstallPath 2^>nul') do set "RHOME=%%B"
)
if not defined RSCRIPT if not defined RHOME (
  for /f "skip=2 tokens=2,*" %%A in ('reg query "HKCU\SOFTWARE\R-core\R" /v InstallPath 2^>nul') do set "RHOME=%%B"
)
if not defined RSCRIPT if not defined RHOME (
  for /d %%D in ("C:\Program Files\R\R-*") do set "RHOME=%%D"
)
if not defined RSCRIPT if defined RHOME (
  if exist "%RHOME%\bin\x64\Rscript.exe" set "RSCRIPT=%RHOME%\bin\x64\Rscript.exe"
  if not defined RSCRIPT if exist "%RHOME%\bin\Rscript.exe" set "RSCRIPT=%RHOME%\bin\Rscript.exe"
)

if not defined RSCRIPT ( where Rscript >nul 2>nul && set "RSCRIPT=Rscript" )

if not defined RSCRIPT (
  for %%C in (
    "%USERPROFILE%\Miniconda3" "%USERPROFILE%\Anaconda3" "%USERPROFILE%\miniforge3" "%USERPROFILE%\mambaforge"
    "%ProgramData%\Miniconda3" "%ProgramData%\Anaconda3"
  ) do (
    if not defined RSCRIPT if exist "%%~C\Scripts\Rscript.exe"          set "RSCRIPT=%%~C\Scripts\Rscript.exe"
    if not defined RSCRIPT if exist "%%~C\Library\bin\Rscript.exe"      set "RSCRIPT=%%~C\Library\bin\Rscript.exe"
    if not defined RSCRIPT if exist "%%~C\envs\r\Scripts\Rscript.exe"   set "RSCRIPT=%%~C\envs\r\Scripts\Rscript.exe"
    if not defined RSCRIPT if exist "%%~C\envs\r\Library\bin\Rscript.exe" set "RSCRIPT=%%~C\envs\r\Library\bin\Rscript.exe"
  )
)

if not defined RSCRIPT (
  echo.
  echo   Could not find an R installation.
  echo     - Install R from https://cran.r-project.org/  ^(or put it on PATH^), OR
  echo     - set R_RSCRIPT to your Rscript.exe and run this again.
  echo.
  pause
  exit /b 1
)

echo.
echo   Starting CellMorphR using:
echo     %RSCRIPT%
echo.
echo   The app will open in your default browser.
echo   KEEP THIS WINDOW OPEN while you work. Closing it stops the app.
echo.

"%RSCRIPT%" --no-init-file "%~dp0launcher.R"

echo.
echo   CellMorphR has stopped.
pause
