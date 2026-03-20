# Spectrogram

A real-time audio spectrogram visualiser CLAP plugin, written in C++.

Built by [Alexander Henderson](https://earthics.org).

## Features

- Real-time FFT-based spectrogram display
- Stereo audio input support
- Embedded GUI (macOS native)

## Building

```sh
git clone git@github.com:jahenderson777/spectrogram.git
cd spectrogram
```

The Objective-C GUI file must be compiled separately before the main CMake build
(see [known issue #2](AGENTS.md)):

```sh
gcc -c src/gui_mac.m -g -o src/gui_mac.o
```

Then configure and build:

```sh
cmake -S . -B build
cmake --build build
```

The built plugin will be at `./build/Spectrogram.clap`.

## Installing

Copy the plugin to your system CLAP folder so your DAW can find it:

```sh
cmake --install build
```

This installs `Spectrogram.clap` to `~/Library/Audio/Plug-Ins/CLAP/`.
After installing, trigger a plugin rescan in your DAW (e.g. Bitwig Studio:
**Settings → Plug-ins → Rescan**).

## Validation

[clap-validator](https://github.com/free-audio/clap-validator) can be used to
check the plugin for correctness before loading it in a DAW:

```sh
clap-validator validate build/Spectrogram.clap
```

A copy of `clap-validator` is kept at `~/Library/Audio/Plug-Ins/CLAP/clap-validator`
for convenience.

## Dependencies

- [CLAP SDK](https://github.com/free-audio/clap)
- [CLAP C++ Helpers](https://github.com/free-audio/clap-helpers)

These are downloaded and configured automatically at build time via CMake's
FetchContent if they are not already available on your system.

## Platform Support

- macOS (primary)

## License

MIT — Copyright © 2025 Alexander Henderson
