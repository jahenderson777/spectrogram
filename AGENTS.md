# AGENTS.md — Spectrogram CLAP Plugin

Guidance for AI agents (and humans) working on this codebase.

---

## What This Project Is

A real-time **spectrogram visualiser** implemented as a [CLAP](https://cleveraudio.org/) audio plugin. It is a pass-through effect: audio flows in and out unmodified, while the plugin GUI displays a scrolling spectrogram of the incoming signal.

- **Plugin ID**: `org.earthics.Spectrogram`
- **Author**: Alexander Henderson / earthics.org
- **Standard**: CLAP 1.x
- **Language**: C++20 + Objective-C (macOS GUI)
- **Platform**: macOS (primary). Win32 and X11 stubs exist but do **not** compile cleanly.

---

## Repository Layout

```
ClapPluginCppTemplate/
├── CMakeLists.txt              # Root build file
├── cmake/
│   ├── Dependencies.cmake      # FetchContent for clap + clap-helpers
│   └── ClapPluginCppTemplate.plist.in   # macOS bundle plist template
└── src/
    ├── Factory.cpp             # CLAP entry point; registers the plugin
    ├── Plugin.h                # Plugin class declaration
    ├── Plugin.cpp              # Plugin logic: audio processing, FFT, painting
    ├── Utils.h                 # dB <-> gain helpers (header-only)
    ├── gui_mac.cpp             # macOS GUI glue (C++, included into Plugin.cpp)
    ├── gui_mac.m               # Cocoa NSView implementation (Objective-C)
    ├── gui_w32.cpp             # Win32 GUI stub (NOT currently buildable)
    ├── gui_x11.cpp             # X11 GUI stub (NOT currently buildable)
    └── tinycolormap.hpp        # Header-only colormap library (MIT, vendored)
```

### Key third-party dependencies (fetched by CMake)

| Library | Source | Role |
|---------|--------|------|
| `clap` | github.com/free-audio/clap | CLAP C API headers |
| `clap-helpers` | github.com/free-audio/clap-helpers | C++ plugin/host wrappers |
| Apple **Accelerate** | system framework | vDSP FFT (macOS only) |
| **tinycolormap** | vendored in `src/` | Turbo/Magma/Jet colormaps |

---

## Building

### Prerequisites

- macOS with Xcode command-line tools
- CMake ≥ 3.25
- Internet access on first build (FetchContent clones clap and clap-helpers)

### Step 1 — Compile the Objective-C file

`CMakeLists.txt` currently adds `src/gui_mac.o` as a **precompiled object**
rather than letting CMake compile `gui_mac.m` directly (known issue #2 below).
You must regenerate this object whenever `gui_mac.m` changes:

```bash
# From the repo root
gcc -c src/gui_mac.m -g -o src/gui_mac.o
```

### Step 2 — Configure and build

```bash
cmake -S . -B build
cmake --build build
```

Output: `build/ClapPluginCppTemplate.clap/` (macOS bundle)

### Step 3 — Install for testing

```bash
cp -r build/ClapPluginCppTemplate.clap ~/Library/Audio/Plug-Ins/CLAP/
```

---

## Architecture

### Audio thread (`Plugin::process`)

Called by the host on the real-time audio thread.  It:

1. Passes audio through unchanged (input → output, stereo).
2. Accumulates incoming samples (left channel only) into `fftBuffer`.
3. When enough samples arrive for a full FFT frame:
   - Applies a **Hann window** via `vDSP_hann_window`.
   - Runs a **real-to-complex FFT** (`vDSP_fft_zrip`, size 2048).
   - Converts to power in dB (`vDSP_zvmags` + `vDSP_vdbcon`).
   - Calls `paintInterpolatedVerticalLines` — paints directly into the GUI
     pixel buffer (3 columns per frame, interpolated from the previous frame).
   - Slides `fftBuffer` forward by `fftOverlap` (512 samples = 75% overlap).
4. On MIDI note-on, the display resets (black fill, `fftX = 0`).

**FFT parameters:**

| Parameter | Value |
|-----------|-------|
| FFT size | 2048 samples |
| Overlap | 512 samples (75%) |
| Window | Hann (normalised) |
| Columns per FFT | 3 (interpolated) |

### GUI / paint thread

The plugin spawns a background thread (`animationThread`) at 10 fps that calls
`GUIPaint`, which forwards to `plugin->paint()` and then schedules a Cocoa
`setNeedsDisplay`.  `Plugin::paint` currently only acquires `bufferMutex`; the
actual spectrogram data is written by the audio thread directly into the pixel
buffer inside `process()`.

### Pixel buffer

`gui->bits` is a flat `uint32_t[GUI_WIDTH * GUI_HEIGHT]` (512×512 = 262 144
pixels) in **BGRA** byte order (matching Cocoa's `NSDrawBitmap` convention).
Pixel addressing: `bits[(GUI_HEIGHT - y - 1) * GUI_WIDTH + x]` (y=0 is bottom).

### Frequency scale

The y-axis uses a **mel-like scale** defined by:

```
mel  = 901 * log10(1 + hz)
hz   = 10^(mel/901) - 1
```

Range: 25 Hz (bottom) → 20 000 Hz (top).

### Colormap

`tinycolormap::ColormapType::Turbo` is used (Google's Turbo colormap).
`normalizedMag` is computed as:

```cpp
(1.2f + y * 0.4f / GUI_HEIGHT) * mag / 140.0f + 0.4f
```

This adds a frequency-dependent brightness boost and a constant floor offset.

---

## Status of Known Issues

### ✅ Fixed

1. **Stray `BB` prefix on `Plugin.cpp` line 1** — removed, committed
   (`279dcdc`).

### 🔴 Open — Build / Infrastructure

2. **Precompiled `gui_mac.o` in CMake**: `CMakeLists.txt` references
   `src/gui_mac.o` as a pre-built object.  Changes to `gui_mac.m` are silently
   ignored until the object is manually rebuilt.
   **Fix**: add `enable_language(OBJC)` to CMakeLists.txt and list
   `src/gui_mac.m` directly in `add_library(...)`, removing the
   `set(GUI_OBJ ...)` workaround.

### 🔴 Open — Correctness

3. **Thread safety of the pixel buffer**: The audio thread writes to
   `gui->bits` (holding `bufferMutex`), and the Cocoa `drawRect` reads from it
   on the main thread without any lock.
   **Fix**: use a double-buffer (ping-pong), or extend `bufferMutex` coverage
   to include `drawRect` reads.

4. **`fftBuffer` only resized in `guiCreate`**: If the host calls `process()`
   before `guiCreate()` (legal in CLAP), `fftBuffer` is empty.  The
   `if (gui && gui->bits)` guard prevents a crash but silently drops audio data.
   **Fix**: resize `fftBuffer` in `activate()` and reset `bufferIndex = 0` there.

5. **`std::async` misuse in animation loop**: `startAnimationLoop` fires
   `GUIPaint` via `std::async(std::launch::async, ...)` but discards the
   `std::future`.  The discarded future's destructor blocks until completion,
   which can cause overlapping paint calls if the frame takes longer than the
   interval.
   **Fix**: call `GUIPaint` directly in the loop body (it only calls
   `setNeedsDisplayInRect`, which is cheap and thread-safe on macOS).

### 🟡 Open — Dead Code / Cleanliness

6. **`audioBuffer` is allocated but never used**: A 32 768-sample circular
   buffer and its mutex exist but audio samples go directly into `fftBuffer`.
   **Fix**: remove `audioBuffer`, `bufferSize`, and the related `bufferMutex`
   lock in `paint()` (replace with a dedicated `fftMutex` if needed).

7. **`filter` / `createLowPassFIRFilter` is unused**: A 101-tap FIR filter is
   computed in the constructor but the `vDSP_desamp` call is commented out.
   **Fix**: remove both if there's no near-term plan to use them.

8. **Large pile of commented-out code**: Several alternative
   `paintVerticalLine` implementations, a waveform renderer, and a multi-band
   FFT sketch remain in `Plugin.cpp`.
   **Fix**: delete them (they're in git history if ever needed).

9. **Unused variables generating compiler warnings**: `note`, `nextEvent`,
   `prevHighMag`, `currHighMag`, `verticalScale` — all flagged by `-Wunused`.
   **Fix**: remove or use each one.

### 🟡 Open — Missing Features

10. **Gain parameter is a no-op**: The plugin exposes a "Gain" parameter with
    proper dB display but never applies it to the audio signal.
    **Fix**: after the pass-through copy in `process()`, multiply output samples
    by `gain_` (use `vDSP_vsmul`).

11. **Win32 / X11 GUI stubs are broken**: `gui_w32.cpp` and `gui_x11.cpp`
    reference `MyPlugin`, `PluginPaint`, `PluginProcessMouseDrag`, etc. —
    names that no longer exist.  They are copy-pasted from an older template.

---

## Recommended Next Steps (Priority Order)

1. **Fix CMake to compile `gui_mac.m` directly** (issue #2) — eliminates the
   fragile manual compile step and makes the build fully reproducible.

2. **Clean up dead code** (issues #6, #7, #8) — makes the file readable.

3. **Fix compiler warnings** (issue #9) — remove unused variables.

4. **Move `fftBuffer` resize to `activate()`** (issue #4).

5. **Fix the animation loop** (issue #5) — drop the `std::async` wrapper.

6. **Add double-buffering** (issue #3) — the most involved fix; needed for
   correctness under a strict host.

---

## Code Style

- C++20; 4-space indentation.
- Private member variables use a trailing underscore (`gain_`); `fftX` is an
  exception that should eventually become `fftX_`.
- No exceptions in audio-thread code (`noexcept` everywhere in CLAP callbacks).
- FFT/vDSP buffers use raw `malloc`/`free` for alignment compatibility with
  vDSP — keep that pattern for vDSP data, but prefer `std::vector` everywhere
  else.
- Prefer `static_cast` over C-style casts.

---

## Testing

No automated test suite.  Manual procedure:

1. Build and install (see "Building" above).
2. Open a CLAP host (e.g. [Bitwig Studio](https://www.bitwig.com/),
   [REAPER](https://www.reaper.fm/) with the CLAP extension, or
   [clap-info](https://github.com/free-audio/clap-info) for quick validation).
3. Insert the plugin on an audio track and verify:
   - The spectrogram scrolls left-to-right as audio plays.
   - Frequencies appear on a mel-like scale (bass at bottom, highs at top).
   - A MIDI note-on event clears and resets the display.
   - Audio passes through without audible artefacts.

---

## File-by-File Reference

| File | Purpose | Thread context |
|------|---------|----------------|
| `Factory.cpp` | CLAP DLL entry; creates `Plugin` instances | init/host thread |
| `Plugin.h` | Class declaration, FFT state, constants | — |
| `Plugin.cpp` | All plugin logic: audio, FFT, painting | audio + GUI threads |
| `Utils.h` | dB/gain conversions; pure functions, no state | any |
| `gui_mac.cpp` | Bridges `Plugin` ↔ Cocoa; `#include`d into `Plugin.cpp` | GUI thread |
| `gui_mac.m` | `NSView` subclass: pixel-buffer rendering + mouse input | main/GUI thread |
| `tinycolormap.hpp` | Lookup-table colormaps (Turbo, Magma, Jet…); header-only | any |
| `cmake/Dependencies.cmake` | FetchContent for clap + clap-helpers | build time |
| `cmake/ClapPluginCppTemplate.plist.in` | macOS bundle metadata template | build time |
