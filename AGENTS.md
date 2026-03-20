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

### Compile the Objective-C file first (one-time step)

The macOS GUI is written in Objective-C (`gui_mac.m`).  CMakeLists.txt currently
adds `gui_mac.o` as a **precompiled object** rather than letting CMake compile
the `.m` file directly.  You must regenerate this object whenever `gui_mac.m`
changes:

```bash
# From the repo root
clang -x objective-c -fobjc-arc -c src/gui_mac.m -o src/gui_mac.o
```

### Configure and build

```bash
mkdir -p build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Debug
make -j$(sysctl -n hw.ncpu)
```

Output: `build/ClapPluginCppTemplate.clap/` (macOS bundle)

### Install for testing

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
actual spectrogram data is written by the audio thread.

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

This adds a frequency-dependent brightness boost and a floor offset.

---

## Known Issues & Technical Debt

These are things that should be fixed before treating the plugin as production-ready:

1. **Typo in Plugin.cpp line 1**: `BB#include "Plugin.h"` — the leading `BB`
   is stray text and must be removed.  The file currently compiles only because
   the precompiled `gui_mac.o` (which includes this file via CMake) was built
   before the typo was introduced, or the compiler silently skips it somehow.
   **Fix**: delete the `BB` prefix.

2. **Precompiled `gui_mac.o` in CMake**: `CMakeLists.txt` references
   `src/gui_mac.o` as a precompiled object.  This means changes to `gui_mac.m`
   are silently ignored until the object is manually rebuilt.  **Fix**: teach
   CMake to compile the `.m` file directly (see "Recommended Fixes" below).

3. **Win32 / X11 stubs are broken**: `gui_w32.cpp` and `gui_x11.cpp` reference
   `MyPlugin`, `PluginPaint`, `PluginProcessMouseDrag`, etc. — names that no
   longer exist.  They are copy-pasted from a template and have not been adapted
   to the `Plugin` class.

4. **Gain parameter is a no-op**: The plugin exposes a "Gain" parameter to the
   host (automatable, with dB display) but never applies it to the audio.

5. **`audioBuffer` is unused**: A 32 768-sample circular buffer is allocated and
   its mutex is locked during FFT, but audio samples are written to `fftBuffer`
   instead.  `audioBuffer` can be removed.

6. **`fftBuffer` only resized in `guiCreate`**: If the host calls `process()`
   before `guiCreate()` (legal in CLAP), `fftBuffer` is empty and `gui` is
   null — the `if (gui && gui->bits)` guard prevents a crash, but FFT data is
   lost until the GUI opens.  **Fix**: resize `fftBuffer` in `activate()`.

7. **`filter` / `createLowPassFIRFilter` is unused**: A 101-tap windowed-sinc
   FIR filter is allocated and computed but the `vDSP_desamp` call is commented
   out.

8. **Thread safety of the pixel buffer**: The audio thread writes to
   `gui->bits` (holding `bufferMutex`), and the GUI thread reads from it
   (without holding any lock) inside `MacPaint → drawRect`.  A double-buffer
   or read lock is needed for correctness.

9. **`std::async` inside animation thread**: `startAnimationLoop` launches
   `GUIPaint` via `std::async(std::launch::async, ...)` but discards the
   returned `std::future`.  The future destructor blocks until the async task
   finishes, which can cause the paint call to overlap with the next tick if
   painting is slow.  **Fix**: call `GUIPaint` directly, or use a `dispatch`
   queue on macOS.

10. **Large pile of commented-out code**: Several alternative `paintVerticalLine`
    implementations, a waveform renderer, and a multi-band FFT sketch are left
    commented out in `Plugin.cpp`.  These should be either deleted or moved to
    a separate branch.

---

## Recommended Fixes (Priority Order)

### 1. Fix the `BB` typo in Plugin.cpp

```cpp
// Line 1 of Plugin.cpp — change:
BB#include "Plugin.h"
// to:
#include "Plugin.h"
```

### 2. Compile gui_mac.m via CMake

Add to `CMakeLists.txt`:

```cmake
# Enable Objective-C
enable_language(OBJC)

add_library(${PROJECT_NAME} MODULE
    src/Utils.h
    src/Plugin.h
    src/Plugin.cpp
    src/Factory.cpp
    src/gui_mac.m          # ← compile directly, not as precompiled object
)
```

Remove the `set(GUI_OBJ src/gui_mac.o)` line and the `${GUI_OBJ}` reference.

### 3. Move fftBuffer resize to activate()

```cpp
bool Plugin::activate(double sampleRate, uint32_t, uint32_t maxFrameCount) noexcept {
    fftBuffer.resize(fftSize, 0.0f);   // ← add this
    bufferIndex = 0;
    // ... existing fade-length calculation ...
    return true;
}
```

### 4. Apply gain in process()

```cpp
// After the pass-through copy loop in process():
if (gain_ != 1.0) {
    for (uint32_t ch = 0; ch < outputChannelsCount; ++ch)
        vDSP_vsmul(output[ch], 1, &gain_, output[ch], 1, process->frames_count);
}
```

---

## Code Style

- C++20; 4-space indentation.
- Member variables: trailing underscore for private (`gain_`, `fftX` is an
  exception — should be `fftX_`).
- No exceptions in audio-thread code (`noexcept` everywhere in CLAP callbacks).
- Prefer `std::` containers over raw `malloc`/`free` in new code; the existing
  FFT buffers use `malloc` for vDSP alignment compatibility — keep that pattern
  for vDSP data, but document it.
- Prefer `static_cast` over C-style casts.

---

## Testing

There is no automated test suite.  Manual testing procedure:

1. Build the plugin (see "Building" above).
2. Install to `~/Library/Audio/Plug-Ins/CLAP/`.
3. Open a CLAP host (e.g. [Clap-Info](https://github.com/free-audio/clap-info),
   [Bitwig Studio](https://www.bitwig.com/), or
   [REAPER](https://www.reaper.fm/) with the CLAP extension).
4. Insert the plugin on an audio track and verify:
   - The spectrogram scrolls left-to-right as audio plays.
   - Frequencies are displayed on a mel-like scale (low frequencies at bottom,
     high at top).
   - A MIDI note-on event resets the display.
   - Audio passes through without audible artifacts.

---

## File-by-File Reference

| File | Purpose | Thread context |
|------|---------|----------------|
| `Factory.cpp` | CLAP DLL entry; creates `Plugin` instances | init/host thread |
| `Plugin.h` | Class declaration, FFT state, constants | — |
| `Plugin.cpp` | All plugin logic including audio processing, FFT, painting | audio + GUI threads |
| `Utils.h` | dB/gain conversions; pure functions, no state | any |
| `gui_mac.cpp` | Bridges Plugin ↔ Cocoa; `#include`d into Plugin.cpp | GUI thread |
| `gui_mac.m` | `NSView` subclass for pixel-buffer rendering and mouse input | main/GUI thread |
| `tinycolormap.hpp` | Lookup-table colormaps; header-only, no state | any |
| `cmake/Dependencies.cmake` | FetchContent declarations for clap + clap-helpers | build time |
| `cmake/ClapPluginCppTemplate.plist.in` | macOS bundle metadata template | build time |
