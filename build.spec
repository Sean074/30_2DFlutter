# PyInstaller spec for flutter2d — single-file console executable.
# Build on Mac to get the Mac binary; build on Windows to get the Windows binary.
# Usage: pyinstaller build.spec --clean --noconfirm

block_cipher = None

a = Analysis(
    ["main.py"],
    pathex=["."],
    binaries=[],
    datas=[],
    hiddenimports=[
        "numpy",
        "matplotlib",
        "matplotlib.backends.backend_tkagg",
        "matplotlib.backends.backend_agg",
        "tkinter",
        "_tkinter",
        "analysis.flutter",
        "file_io.loader",
        "file_io.exporter",
        "plotting.plots",
        "ui.tui",
    ],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=["pytest", "black", "ruff", "rich", "questionary", "flask"],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.zipfiles,
    a.datas,
    [],
    name="flutter2d",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=True,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
