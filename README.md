# Litmus

Litmus is an automated molecular integrals generator.

## Building

Litmus uses CMake (>= 3.16) and a C++17 compiler. To build:

```bash
cmake -S . -B build              # configure (defaults to a Release build)
cmake --build build -j           # compile, produces build/litmus.x
```

To make a debug build instead, configure with
`cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug`.

## Running

The run is configured in `src/litmus.cpp` (integral type, maximum angular
momentum, geometrical-derivative orders). Adjust it, rebuild, then run:

```bash
./build/litmus.x
```

Generated integral source files are written to the current working directory.

## Testing

Unit tests use [GoogleTest](https://github.com/google/googletest), fetched
automatically at configure time (requires network access on the first
configure). After building, run the suite with CTest:

```bash
ctest --test-dir build --output-on-failure
```

Tests live under `tests/`. Pass `-DLITMUS_BUILD_TESTS=OFF` at configure time to
skip building them (e.g. for offline or packaging builds).
