#!/bin/bash
# Build the macOS single-file executable.
# Run this on a Mac with: bash build.sh
# Output: dist/flutter2d

set -e

if ! python -c "import tkinter" 2>/dev/null; then
    echo "WARNING: this Python has no tkinter. The binary will plot with the native macOS"
    echo "         window backend instead. Install python-tk for the Tk window."
fi

pip install pyinstaller --quiet
pyinstaller build.spec --clean --noconfirm

echo ""
echo "Build complete: dist/flutter2d"
echo "Test with: ./dist/flutter2d"
