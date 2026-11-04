@echo off
setlocal
cd /d "%~dp0"
set "VENV=%~dp0.venv"
set "PYW="
if exist "%VENV%\Scripts\pythonw.exe" set "PYW=%VENV%\Scripts\pythonw.exe"
if not defined PYW call :makevenv
if defined PYW goto :launch

rem ---- fallback: any Python on PATH (only if the private .venv could not be created) ----
if not defined PYW (
  where pythonw >nul 2>nul && set "PYW=pythonw"
)
if not defined PYW (
  echo Python was not found. Install Python 3.11 - 3.13 from python.org and run this file again.
  pause
  exit /b 1
)
echo WARNING: using the fallback Python "%PYW%" because .venv could not be created.
set "PYC=%PYW:pythonw=python%"
"%PYC%" -c "import gemmi,statsmodels,scipy,pandas,matplotlib,requests,tkinter" >nul 2>nul
if errorlevel 1 (
  echo Installing required packages once - this can take a few minutes ...
  "%PYC%" -m pip install -r "%~dp0requirements.txt"
)

:launch
start "" "%PYW%" "%~dp0gate_gui.py"
endlocal
exit /b 0

:makevenv
echo Creating the private Python environment in "%VENV%" - first run only, needs internet, a few minutes ...
set "BASEPY="
for %%P in ("%LOCALAPPDATA%\Programs\Python\Python313\python.exe" "%LOCALAPPDATA%\Programs\Python\Python312\python.exe" "%LOCALAPPDATA%\Programs\Python\Python311\python.exe") do (
  if not defined BASEPY if exist %%P set BASEPY="%%~P"
)
if not defined BASEPY (
  where py >nul 2>nul && set "BASEPY=py -3"
)
if not defined BASEPY (
  where python >nul 2>nul && set "BASEPY=python"
)
if not defined BASEPY goto :eof
%BASEPY% -m venv "%VENV%"
if errorlevel 1 goto :venvfail
"%VENV%\Scripts\python.exe" -m pip install --disable-pip-version-check -r "%~dp0requirements.txt"
if errorlevel 1 goto :venvfail
"%VENV%\Scripts\python.exe" -c "import numpy,scipy,pandas,statsmodels,matplotlib,gemmi,requests,tkinter"
if errorlevel 1 goto :venvfail
set "PYW=%VENV%\Scripts\pythonw.exe"
goto :eof

:venvfail
echo Could not build .venv - see the messages above.
if exist "%VENV%" rmdir /s /q "%VENV%"
set "PYW="
goto :eof
