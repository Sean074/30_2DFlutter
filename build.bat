@echo off
REM Build the Windows single-file executable.
REM Run this on Windows with: build.bat
REM Output: dist\flutter2d.exe

python -c "import tkinter" 2>nul
if errorlevel 1 (
    echo WARNING: this Python has no tkinter. The binary will not open a plot window;
    echo          it will only save the PNG. Reinstall Python with the tcl/tk option.
)

pip install pyinstaller --quiet
pyinstaller build.spec --clean --noconfirm

echo.
echo Build complete: dist\flutter2d.exe
echo Test with: dist\flutter2d.exe
