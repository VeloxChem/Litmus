# Litmus

Litmus is an automated molecular integrals generator.

New to the codebase? Start with [ONBOARDING.md](ONBOARDING.md) for an
architecture map, conventions, and pitfalls.

## Building

Litmus uses CMake (>= 3.16) and a C++17 compiler. To build:

```bash
cmake -S . -B build              # configure (defaults to a Release build)
cmake --build build -j           # compile, produces build/litmus.x
```

To make a debug build instead, configure with
`cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug`.

## Running

A run is described by a config file (a minimal TOML subset) and launched with
the `run` subcommand:

```bash
./build/litmus.x run examples/t2c_geom.toml
```

The config selects the integral family (`type`), the maximum angular momentum
(`lmax`), the operator label (`integral`), the geometric-derivative orders
(`geom`), and a few family-specific keys. See `litmus.x --help` for the full
schema and the ready-to-run samples under [`examples/`](examples). A minimal
config looks like:

```toml
type     = "t2c_geom_cpu"   # run-type family
lmax     = 2                # maximum angular momentum
integral = "none"           # operator label
geom     = [0, 2, 0]        # geometric-derivative orders
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
