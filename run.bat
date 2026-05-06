@echo off
echo Installing / updating dependencies...
python -m pip install -q PyQt6 numpy matplotlib scipy
echo.
echo Launching RIS CubeSat Simulator...
python ris_simulator.py
pause
