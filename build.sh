#!/bin/bash
# Build the macOS single-file executable.
# Run this on a Mac with: bash build.sh
# Output: dist/flutter2d

set -e

pip install pyinstaller --quiet
pyinstaller build.spec --clean --noconfirm

echo ""
echo "Build complete: dist/flutter2d"
echo "Test with: ./dist/flutter2d"
