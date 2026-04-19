# TUI Style Guide — 2D Flutter Analysis

## Design Philosophy

This is a 1990s-style command-line engineering program. The interface is text-only,
functional, and direct. No decoration for its own sake. Every element on screen
earns its place by communicating information.

---

## Visual Language

### Colour Palette

Restrict to colours available on a standard 16-colour ANSI terminal.

| Role              | Colour       | ANSI Code |
|-------------------|--------------|-----------|
| Normal text       | White/Grey   | Default   |
| Section headers   | Cyan (bold)  | `\033[1;36m` |
| Prompts           | White (bold) | `\033[1;37m` |
| Errors            | Red          | `\033[0;31m` |
| Warnings          | Yellow       | `\033[0;33m` |
| Confirmed/OK      | Green        | `\033[0;32m` |
| Calculated output | Cyan         | `\033[0;36m` |

No background colours. No blinking. No 256-colour extensions.

### Typography

- All caps for section titles and field labels.
- Monospaced alignment: pad labels to a fixed width (e.g. 28 chars) so values line up.
- Use `=` for dividers between major sections. Use `-` for minor sub-dividers.
- Width: target 72 characters. Never exceed 80.

```
================================================================================
  2D AEROELASTIC FLUTTER ANALYSIS                            p-k METHOD
================================================================================
```

---

## Workflow Order

The screen sequence mirrors the engineering workflow:

```
1. INPUT          — user provides geometry, structural, and aero parameters
2. CHECK INPUT    — program echoes values back; user confirms before proceeding
3. OUTPUT         — program runs analysis and prints results
4. CHECK OUTPUT   — summary flags, warnings, and flutter speed result
```

Do not skip stages. Do not combine stages on the same screen.

---

## Stage 1 — INPUT

Present one logical group of fields at a time. Group order:

1. **Geometry** — semi-chord, elastic axis, CG offset
2. **Mass properties** — mass, moment of inertia
3. **Structural** — natural frequencies, damping
4. **Aerodynamic** — air density, velocity range

Each field:

```
  SEMI-CHORD b [m]             : _
```

Rules:
- Label is left-aligned, padded to 28 chars.
- Unit is shown in brackets.
- Colon and cursor on same line.
- One field per line.
- Blank line between groups.
- Group heading in caps on its own line, no bracket or colon.

```
  GEOMETRY
  --------
  SEMI-CHORD b [m]             : _
  ELASTIC AXIS a [-1 to 1]     : _
  CG OFFSET x_a [-1 to 1]      : _
```

---

## Stage 2 — CHECK INPUT

After all input is collected, print a full echo of every value before running.

```
================================================================================
  INPUT SUMMARY — CONFIRM BEFORE RUNNING
================================================================================

  GEOMETRY
  --------
  SEMI-CHORD b [m]             : 0.5000
  ELASTIC AXIS a [-1 to 1]     : -0.2000
  CG OFFSET x_a [-1 to 1]      :  0.1000

  ...

  Proceed? [Y/n] : _
```

Rules:
- Values are right-aligned to match their column.
- Values printed to 4 decimal places (or engineering notation for very small/large).
- Any value outside expected bounds is highlighted in yellow with a warning tag.
- User must press Y (or Enter) to continue, or N to re-enter.

---

## Stage 3 — OUTPUT

Print progress and results as the solver runs. Do not buffer silently.

```
================================================================================
  RUNNING ANALYSIS
================================================================================

  Velocity sweep : 0.0 to 100.0 m/s  (50 steps)
  Mode 1 of 2 ....... done
  Mode 2 of 2 ....... done

  Elapsed : 0.34 s

================================================================================
  RESULTS
================================================================================

  V [m/s]       g (MODE 1)    f [rad/s]     g (MODE 2)    f [rad/s]
  ----------    ----------    ----------    ----------    ----------
     10.000        -0.0312       24.512        -0.0801       61.340
     20.000        -0.0274       24.609        -0.0712       61.280
     ...
```

Rules:
- Column headers are short and include units.
- Columns are separated by fixed-width gaps (4 spaces minimum).
- Data is right-aligned within each column.
- Highlight (cyan) any row where `g >= 0` — flutter has occurred.

---

## Stage 4 — CHECK OUTPUT

After results table, print a summary block.

```
================================================================================
  ANALYSIS SUMMARY
================================================================================

  Flutter detected     : YES
  Flutter speed (est.) :  42.60 m/s      (MODE 1, g=0 crossing)
  g limit (0.03)       :  38.10 m/s      (MODE 1)

  Plot saved           : flutter_output.png

  [W] WARNING: reduced frequency k < 0.01 at V > 80 m/s — quasi-steady limit
================================================================================
```

Rules:
- Binary facts (`YES/NO`, `DETECTED/NOT DETECTED`) before numeric results.
- Each warning on its own line, prefixed `[W]`.
- Each error prefixed `[E]` in red.
- End with the closing `===` bar.

---

## Error and Warning Messages

Always inline; never pop-up, never modal.

```
  [E] INVALID INPUT: b must be > 0. Re-enter.
  [W] x_a outside typical range [-0.5, 0.5]. Confirm value is correct.
```

- `[E]` = fatal for that field; re-prompt immediately.
- `[W]` = advisory; shown at check-input and check-output stages.
- No stack traces visible to the user. Log them to file if needed.

---

## Interaction Rules

- Default values shown in brackets: `[0.5]`
- Press Enter to accept default.
- `q` or `Q` at any prompt exits cleanly with a one-line goodbye.
- No mouse input.
- No ncurses / curses screen-clearing mid-workflow; scroll downward only.

---

## Do Not Use

- Spinners, progress bars with ASCII animation
- Box-drawing characters (`┌`, `│`, `─`, etc.)
- Emoji or unicode symbols beyond standard ASCII
- Colour gradients or 256-colour codes
- Interactive menus or arrow-key navigation
- Any library that requires a full terminal (curses, rich, textual)
