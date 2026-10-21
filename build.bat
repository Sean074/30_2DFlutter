@echo off
REM Build the Windows single-file executable.
REM Run this on Windows with: build.bat
REM Output: dist\flutter2d.exe

pip install pyinstaller --quiet
pyinstaller build.spec --clean --noconfirm

echo.
echo Build complete: dist\flutter2d.exe
echo Test with: dist\flutter2d.exe
