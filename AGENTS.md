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
spectrogram_plugin/
├── CMakeLists.txt              # Root build file
├── AGENTS.md                   # This file
├── README.md                   # User-facing documentation
├── cmake/
│   ├── Dependencies.cmake      # FetchContent for clap + clap-helpers
│   └── Spectrogram.plist.in    # macOS bundle plist template
├── build_number.txt            # Auto-incrementing build number (persisted across builds)
└── src/
    ├── build_number.h          # Generated — #define BUILD_NUMBER N (gitignored)
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

| Library              | Source                             | Role                      |
| -------------------- | ---------------------------------- | ------------------------- |
| `clap`               | github.com/free-audio/clap         | CLAP C API headers        |
| `clap-helpers`       | github.com/free-audio/clap-helpers | C++ plugin/host wrappers  |
| Apple **Accelerate** | system framework                   | vDSP FFT (macOS only)     |
| **tinycolormap**     | vendored in `src/`                 | Turbo/Magma/Jet colormaps |

---

## Building

### Prerequisites

- macOS with Xcode command-line tools
- CMake ≥ 3.25
- Internet access on first build (FetchContent clones clap and clap-helpers)

### First-time configure

Only needed once, or after changing `CMakeLists.txt`:

```bash
cmake -S . -B build
```

### Build and install

The normal day-to-day command — run this after every code change:

```bash
cmake --build build && cmake --install build
```

This compiles the plugin, **auto-increments the build number** (written to
`src/build_number.h` and displayed in the GUI overlay), copies
`Spectrogram.clap` to `~/Library/Audio/Plug-Ins/CLAP/`, and triggers a
plugin rescan in your DAW.  No manual steps are needed for `gui_mac.m` —
CMake compiles Objective-C directly (`LANGUAGES C CXX OBJC`).

Output bundle: `build/Spectrogram.clap/`

**Note on the legacy bundle**: A second copy is kept at
`~/Library/Audio/Plug-Ins/CLAP/ClapPluginCppTemplate.clap` for compatibility
with any Bitwig projects that were saved when the plugin had the old template
name.  It is identical to `Spectrogram.clap` and must be updated manually
alongside it (`cp -r ~/Library/Audio/Plug-Ins/CLAP/Spectrogram.clap ~/Library/Audio/Plug-Ins/CLAP/ClapPluginCppTemplate.clap`).

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
   - Pushes an `FFTFrame` (magnitude vector) onto `fftQueue_` (protected by
     `fftQueueMutex_`) for the paint thread to consume.
   - Slides `fftBuffer` forward by `fftOverlap` (512 samples = 75% overlap).
4. On MIDI note-on, pushes a reset `FFTFrame` (empty magnitudes, `reset=true`)
   onto `fftQueue_` and resets the FFT accumulator.

The audio thread does **no pixel writes** and holds no GUI state.

**FFT parameters:**

| Parameter       | Value             |
| --------------- | ----------------- |
| FFT size        | 2048 samples      |
| Overlap         | 512 samples (75%) |
| Window          | Hann (normalised) |
| Columns per FFT | 3 (interpolated)  |

### Paint thread / `Plugin::paint()`

`Plugin::paint()` is the sole writer of `gui->bits`.  It drains `fftQueue_`,
and for each frame either clears the buffer (reset frame) or calls
`paintInterpolatedVerticalLines` (3 columns, interpolated from the previous
frame) and advances `fftX`.

`paint()` is currently called from two places:
- **`MacInputEvent`** — on mouse events (main thread).
- **`GUIPaint`** — called by the animation thread (see below).

### Animation thread (`startAnimationLoop`)

The plugin spawns a background `std::thread` (`animationThread`) that wakes at
10 fps, calls `GUIPaint` (which calls `paint()` then `MacPaint`), and goes back
to sleep.

> ⚠️ **Currently disabled for crash diagnostics.**  `startAnimationLoop` /
> `stopAnimationLoop` are commented out in `guiCreate` / `guiDestroy`.
> Redraws therefore only happen on mouse events.  See open issue #6.

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

- **Precompiled `gui_mac.o` in CMake**: `CMakeLists.txt` previously referenced
  `src/gui_mac.o` as a pre-built object; changes to `gui_mac.m` were silently
  ignored.  Fixed by adding `LANGUAGES OBJC` to the project and listing
  `src/gui_mac.m` directly in `add_library(...)`.

- **`std::async` misuse in animation loop**: `startAnimationLoop` previously
  fired `GUIPaint` via `std::async(std::launch::async, ...)` and discarded the
  `std::future`, causing the destructor to block and overlapping paint calls.
  Fixed by calling `GUIPaint` directly in the loop body.

- **Audio thread doing pixel writes**: `process()` previously called
  `paintInterpolatedVerticalLines` and managed `fftX` directly.  Replaced with
  a producer/consumer model: the audio thread pushes `FFTFrame` structs onto
  `fftQueue_`; all pixel writes happen in `Plugin::paint()` on the paint thread.

### 🔴 Open — Correctness

6. **Animation thread disabled (crash diagnostics)**: `startAnimationLoop` /
   `stopAnimationLoop` are currently commented out in `guiCreate` /
   `guiDestroy` while investigating a `BitwigAudioEngine` `EXC_BAD_ACCESS`
   crash that appeared after the producer/consumer refactor.  Redraws only
   happen on mouse events until this is resolved and the thread is re-enabled.

7. **`MacPaint` / `MacSetDebugText` called from a background thread**: When the
   animation thread is re-enabled, `GUIPaint` calls both `MacSetDebugText` and
   `MacPaint` (`[NSView setNeedsDisplayInRect:]`) from a non-main thread.
   AppKit requires all view operations on the main thread.
   **Fix**: wrap both calls in `dispatch_async(dispatch_get_main_queue(), ^{ … })`
   inside `gui_mac.m`.

8. **Thread safety of the pixel buffer**: The animation thread writes to
   `gui->bits` inside `Plugin::paint()`, while Cocoa's `drawRect:` reads the
   same buffer on the main thread without any lock.
   **Fix**: use a double-buffer (ping-pong), or guard `drawRect:` reads with
   the same mutex used in `paint()`.

### 🟡 Open — Dead Code / Cleanliness

9. **`audioBuffer` is allocated but never used**: A 32 768-sample circular
   buffer and its mutex exist but audio samples go directly into `fftBuffer`.
   **Fix**: remove `audioBuffer`, `bufferSize`, and the related `bufferMutex`
   lock in `paint()` (replace with a dedicated `fftMutex` if needed).

10. **`filter` / `createLowPassFIRFilter` is unused**: A 101-tap FIR filter is
    computed in `activate()` but the `vDSP_desamp` call is commented out.
    **Fix**: remove both if there's no near-term plan to use them.

11. **Large pile of commented-out code**: Several alternative
    `paintVerticalLine` implementations, a waveform renderer, and a multi-band
    FFT sketch remain in `Plugin.cpp`.
    **Fix**: delete them (they're in git history if ever needed).

12. **Unused variables generating compiler warnings**: `nextEvent`,
    `prevHighMag`, `currHighMag` — flagged by `-Wunused`.
    **Fix**: remove or use each one.

### 🟡 Open — Missing Features

13. **Gain parameter is a no-op**: The plugin exposes a "Gain" parameter with
    proper dB display but never applies it to the audio signal.
    **Fix**: after the pass-through copy in `process()`, multiply output samples
    by `gain_` (use `vDSP_vsmul`).

14. **Win32 / X11 GUI stubs are broken**: `gui_w32.cpp` and `gui_x11.cpp`
    reference `MyPlugin`, `PluginPaint`, `PluginProcessMouseDrag`, etc. —
    names that no longer exist.  They are copy-pasted from an older template.

---

## Recommended Next Steps (Priority Order)

1. **Diagnose and fix the animation-thread crash** (issue #6) — re-enable the
   animation loop once the root cause of the `BitwigAudioEngine` crash is
   confirmed.

2. **Dispatch AppKit calls to the main thread** (issue #7) — wrap
   `setNeedsDisplayInRect:` and `MacSetDebugText` in `dispatch_async` so the
   animation thread is safe to use.

3. **Add double-buffering** (issue #8) — needed for pixel-buffer correctness
   once the animation thread is active again.

4. **Clean up dead code** (issues #9, #10, #11) — makes the file readable.

5. **Fix compiler warnings** (issue #12) — remove unused variables.

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

### Automated — clap-validator

[clap-validator](https://github.com/free-audio/clap-validator) exercises the
plugin's CLAP API conformance without needing a DAW.  A copy lives at
`~/Library/Audio/Plug-Ins/CLAP/clap-validator`.

```bash
~/Library/Audio/Plug-Ins/CLAP/clap-validator validate build/Spectrogram.clap
```

Expected result: **0 failures** (some tests are skipped because the plugin does
not implement optional extensions like `clap.state`).

Run this after any change to `Factory.cpp`, `Plugin.h`, or the CLAP entry point
before loading the plugin in a DAW.

### Manual — Bitwig Studio

1. Build and install (see "Building" above).
2. Open [Bitwig Studio](https://www.bitwig.com/) and rescan plugins
   (**Settings → Plug-ins → Rescan**) if the plugin does not appear.
3. Insert **Spectrogram** (or **ClapPluginCppTemplate** for legacy projects) on
   an audio track and verify:
   - The spectrogram window opens and scrolls left-to-right as audio plays.
   - Frequencies appear on a mel-like scale (bass at bottom, highs at top).
   - A MIDI note-on event clears and resets the display.
   - Audio passes through without audible artefacts.

---

## File-by-File Reference

| File                         | Purpose                                                  | Thread context      |
| ---------------------------- | -------------------------------------------------------- | ------------------- |
| `Factory.cpp`                | CLAP DLL entry; creates `Plugin` instances               | init/host thread    |
| `Plugin.h`                   | Class declaration, FFT state, constants                  | —                   |
| `Plugin.cpp`                 | All plugin logic: audio, FFT, painting                   | audio + GUI threads |
| `Utils.h`                    | dB/gain conversions; pure functions, no state            | any                 |
| `gui_mac.cpp`                | Bridges `Plugin` ↔ Cocoa; `#include`d into `Plugin.cpp`  | GUI thread          |
| `gui_mac.m`                  | `NSView` subclass: pixel-buffer rendering + mouse input  | main/GUI thread     |
| `tinycolormap.hpp`           | Lookup-table colormaps (Turbo, Magma, Jet…); header-only | any                 |
| `cmake/Dependencies.cmake`   | FetchContent for clap + clap-helpers                     | build time          |
| `cmake/Spectrogram.plist.in` | macOS bundle metadata template                           | build time          |
