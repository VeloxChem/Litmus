# Litmus — Onboarding & Handoff

This document orients a new developer to the Litmus codebase: what it does, how
it is structured, how to build/run/test it, and the conventions and pitfalls to
keep in mind. For the short build/run/test commands see [README.md](README.md);
this document is the deeper map.

## What Litmus is

Litmus is an **automated molecular integrals generator**. It does not compute
integrals at runtime — it is a *code generator* that emits C++ source files
implementing the recurrence relations for molecular integrals (overlap, kinetic,
nuclear attraction, electron repulsion, multipoles, ECPs, geometric
derivatives, …). Those generated files are consumed by VeloxChem.

The core idea: integral recurrences (Obara–Saika and Hermite/HRR/VRR schemes)
are represented **symbolically** as data structures, manipulated into a
fully-reduced recursion graph, and then printed as C++ code.

```
litmus.cpp ──▶ generators/ ──▶ recursions/ ──▶ algebra/ ──▶ general/
 (config)      (emit C++)      (apply OS rules) (symbolic types) (I/O, strings)
```

## Repository layout

```
src/
  general/      I/O and string helpers (file_stream, string_formater)
  algebra/      symbolic value types and templates (the vocabulary)
  recursions/   Obara–Saika "drivers": apply recurrence rules to build graphs
  generators/   emit C++ source from recursion graphs
  litmus.cpp    the CLI entry point; parses a config file and dispatches
tests/
  algebra/      unit tests for the value types/templates
  general/      unit tests for the utilities
  recursions/   one test file per driver
  generators/   unit tests for the code generators/emitters
CMakeLists.txt  root build; src/* and tests/* have their own
.github/workflows/ci.yml   CI: build + test on ubuntu-latest and macos-latest
```

Rough sizes: `algebra` ~27 files, `recursions` ~111 files (54 drivers + headers
+ defs), `generators` ~138 files, plus a test suite of **472 tests**.

## The algebra module — the vocabulary

Everything is built from small immutable symbolic value types. Read these first
(`src/algebra/`):

- **`TensorComponent`** — a Cartesian component `(ax, ay, az)`; `order()` is the
  angular momentum, `shift(axis, value)` raises/lowers it.
- **`Tensor`** — a shell of a given order; `label()` gives `S/P/D/F/...`.
- **`OperatorComponent` / `Operator`** — an integral operator (name + tensor
  shape + target/center), e.g. `"1"` (overlap), `"T"` (kinetic), `"A"` (nuclear),
  `"1/|r-r'|"` (ERI), `"r"` (multipole), `"p"` (linear momentum), `"U_L"`/`"U_l"`
  (ECP), `"G(r)"`/`"GX(r)"`/`"GR2(r)"` (three-center families).
- **`OneCenter(Component)`**, **`TwoCenterPair(Component)`** — the bra/ket
  expansion centers (one Gaussian, or a pair).
- **`Integral<Bra,Ket>` / `IntegralComponent<Bra,Ket>`** — an integral: bra, ket,
  integrand operator, Boys/recursion `order`, and a vector of **prefix**
  operators (used for geometric derivatives). `Component` = a single Cartesian
  component; the non-component form is a whole shell.
- **`Factor` / `Fraction`** — symbolic coefficients attached to recursion terms
  (e.g. `PA`, `1/eta`, `zeta/b_e^2`) and rational prefactors.
- **`RecursionTerm` / `RecursionExpansion` (a.k.a. *dist*) / `RecursionGroup`** —
  a term is `prefactor · factors · integral`; an expansion is a root term plus
  the terms it expands into; a group is a set of expansions.

These types live in `std::set`/`std::map`, so every `operator<` **must be a
strict weak ordering** — this is load-bearing, not cosmetic.

## The recursions module — the drivers

A *driver* knows one recurrence. Naming encodes the integral arity and scheme:

- `t2c_*` / `t3c_*` / `t4c_*` — **symbolic** term drivers (2/3/4-center). They
  build `RecursionExpansion`s with explicit `Factor`s.
- `v2i_*` / `v3i_*` / `v4i_*` — **enumerators**. They produce the *set* of unique
  integrals a kernel depends on (no factors), used to drive code generation.
- `*_geom_*` — geometric-derivative variants; they operate on integrals that
  carry **prefix** operators, and their `is_*` predicate checks
  `prefixes_order()` (e.g. `{1,0,0,0}`).
- `*_center_*`, `*_hrr_*`, `*_vrr_*`, `*_proj_ecp`, `*_loc_ecp`, `*_trans_*` —
  the remaining schemes.

The two-center types are reused by some "three-center" drivers (the third center
rides on the operator, e.g. `t3c_ovl` uses `R2CTerm` with operator `"G(r)"`).
Type aliases per arity live in `t2c_defs.hpp` / `t3c_defs.hpp` / `t4c_defs.hpp`
(`R2CTerm`, `R2CDist`, `R2Group`, `T2CIntegral`, `I2CIntegral`, `VT2CIntegrals`,
`SI2CIntegrals`, …). **Always check a driver's own `.hpp` for which `*_defs` it
includes** before assuming its types.

Typical driver shape:

```cpp
bool is_overlap(const R2CTerm&) const;                  // gate by operator/prefix
std::optional<R2CDist> bra_vrr(const R2CTerm&, char ax) const;  // one step
R2CDist apply_bra_vrr(const R2CTerm&) const;            // recurse to base
R2Group create_recursion(const VT2CIntegrals&) const;   // whole reduction
```

## Build, run, test

```bash
cmake -S . -B build               # configure (Release by default)
cmake --build build -j            # build → build/litmus.x
ctest --test-dir build --output-on-failure   # run all 472 tests
```

- C++17, CMake ≥ 3.16. GoogleTest is fetched via FetchContent on first configure
  (needs network once). `-DLITMUS_BUILD_TESTS=OFF` skips tests.
- Modules are CMake **OBJECT libraries** (`ltm_general`, `ltm_algebra`,
  `ltm_recursions`, `ltm_generators`); this layout breaks a recursions↔generators
  include cycle. `litmus_headers` is an INTERFACE target exposing the include
  dirs.
- The recursion and generator test targets glob `tests/recursions/*.cpp` and
  `tests/generators/*.cpp` (`CONFIGURE_DEPENDS`), so a new test file is picked up
  on the next configure. The algebra/general targets list files explicitly.

To change what gets generated, write a config file and run
`./build/litmus.x run <config>` — no recompile needed. `src/litmus.cpp` parses
the config (via `cfg::` in `src/general/config.{hpp,cpp}`) and dispatches on the
`type` key to the matching generator; an unknown `type` is a hard error listing
the valid families. See `litmus.x --help` and the samples in `examples/`. Output
files land in the current working directory.

The config format is a minimal TOML subset parsed by hand (zero dependencies):
`key = value` lines, `#` comments, and values that are quoted/bare strings,
integers, booleans, or `[1, 2, 3]` integer arrays. The raw parser is generic
(`cfg::Config`, typed getters with fallbacks); two schemas sit on top of it.

**Legacy schema** (the original 13 generator families). Keys: `type` (required),
`lmax`, `integral`, `geom` (arity 3/4/5 per family), `aux_lmax` (t3c),
`proj_lmax` (proj-ecp), `rec_form` (t2c, `[1, 0]`), `use_rs` (t2c/g2c).

**New-style schema** (`cfg::RunConfiguration` in
`src/general/run_configuration.{hpp,cpp}`) decomposes the monolithic `type` into
orthogonal, typed dimensions for the next generation of generators. Keys:
`integral_type` (required: `two_center`/`three_center`/`four_center`, also
`2c`/`3c`/`4c`), `max_ang_mom` (required), `min_ang_mom` (default 0),
`operator_type` (default `overlap`; the integrand — `overlap`, `kinetic_energy`,
`nuclear_potential`, `electron_repulsion`, `dipole_momentum`, `linear_momentum`,
`local_ecp`, `projected_ecp`, `three_center_overlap`, `three_center_r2`,
`three_center_r_dot_r2`, with `to_string` round-tripping to the generator label),
`hardware` (default `cpu`), `language` (default `C++`), `storage_form` (default
`VeloxChemSparse`), `signature` (default `VeloxChemScreened`). Each enumerated
field is validated against its allowed spellings (case/`_`/`-` insensitive) and
the angular-momentum range is checked. `litmus run` recognizes a config as
new-style when it carries an `integral_type` key; the first generator that
consumes `RunConfiguration` now exists — see *The new-style two-center
generator* below. See also `examples/four_center.toml`.

## The new-style two-center generator (work in progress)

`TwoCenterGenerator` (`src/generators/two_center_generators.{hpp,cpp}`) is the
first generator built on `cfg::RunConfiguration`. `litmus run` dispatches a
new-style config with `integral_type = two_center` to it (three/four-center
still print a "not wired in yet" notice). It currently supports the operators
`overlap`, `kinetic_energy`, and `electron_repulsion`; every other
`OperatorType` raises a `cfg::ConfigError`. The operator `switch`es deliberately
carry no `default`, so adding an `OperatorType` enumerator trips `-Wswitch` at
each dispatch site that needs updating.

`generate()` loops the bra (A) and ket (B) over `[min_ang_mom, max_ang_mom]` and,
for each target `(a|o|b)`, splits the work into three integral sets:

- **HRR transfer group** (`_generate_hrr_integral_group`) — the momentum-transfer
  triangle that grows the *smaller* side up to the target. For `a <= b` it grows
  the bra (bra `i = 1..a`, ket `j = b .. a+b-i`); for `a > b` it grows the ket
  symmetrically.
- **VRR base group** (`_generate_vrr_base_integral_group`) — the seed ladder HRR
  consumes, keeping the smaller side at 0: `(0|o|b)...(0|o|a+b)` when `a <= b`,
  else `(a|o|0)...(a+b|o|0)`. These are "the integrals output from HRR", fed into
  the VRR recursion.
- **Full VRR group** (`_generate_vrr_integral_group`) — the base fed into the real
  vertical-recursion drivers (`V2IOverlapDriver`; `V2IKineticEnergyDriver` then a
  `V2IOverlapDriver` closure; `V2IElectronRepulsionDriver`). This is the
  recursion closure of the base, so the base is a subset of it; the *remaining*
  VRR integrals (full minus base) are what VRR generates to evaluate the seeds.

The drivers and their chaining mirror `T2CCPUGenerator::_generate_integral_group`
(check that method when adding an operator). Angular momentum is read off an
integral with `integral[0]` (bra) and `integral[1]` (ket); the integrand is
`integral.integrand()`.

For each target `generate()` hands the three groups to a **two-center emitter**
(`src/generators/two_center_emitters.{hpp,cpp}`) which writes a `.hpp`/`.cpp`
pair per integral. The emitter is selected by `make_two_center_emitter()`, a
factory that dispatches on the config's `hardware` × `language` (nested `switch`es
with no `default`, so an unsupported target trips `-Wswitch` and is otherwise
rejected with a `cfg::ConfigError`). The only emitter today is
`CppCpuTwoCenterEmitter` (C++ on CPU). `generate()` still prints a one-line
summary per target to stdout (`Generated SP kernel (0 HRR, 1 VRR base, 3 VRR
rest)`). Try:

```bash
printf 'integral_type="two_center"\noperator_type="kinetic_energy"\nmax_ang_mom=2\n' > /tmp/t2c.toml
cd /tmp && /path/to/litmus.x run /tmp/t2c.toml   # writes ObaraSaikaTwoCenter*.{hpp,cpp} into the cwd
```

**What the emitter produces.** Files are named
`ObaraSaikaTwoCenter<Operator><BraKet>.{hpp,cpp}` (e.g.
`ObaraSaikaTwoCenterKineticEnergyPP`). Each `.hpp` is a guarded declaration in the operator's
namespace — `os2c::ovl` (overlap), `os2c::kin` (kinetic energy), `os2c::eri`
(electron repulsion) — of a `compute_<la>_<lb>` kernel (lowercase shell labels,
e.g. `compute_p_p`); the `(bra|op|ket)^order` caption is used for the `@brief`
doc line. The `.cpp` carries the matching definition whose body now emits the
integral computation **workflow** as real C++ calls (not comment placeholders),
in order:

- `npairs = pair.number_of_pairs()` and a zeroed result buffer
  `CDenseMatrix((2*la+1)*(2*lb+1), npairs)`.
- the Obara–Saika distances for the side the VRR builds —
  `osfunc::compute_pb` when the seeds grow the ket (`(0|o|b)`), `compute_pa` when
  they grow the bra (`(b|o|0)`); a pure `(0|o|0)` target needs neither. The build
  side is read off the integral set by `vrr_build_side`.
- the primitive-overlap seed `povl_ss = osfunc::compute_overlap(pair)` (overlap
  and kinetic only — see `seeds_with_primitive_overlap`).
- one VRR step per non-seed VRR integral, lowest total momentum first:
  `p<abbr>_<shells> = os2c::vrr::<abbr>::compute_<L>(pair, pa|pb, <two lower
  integrals>)`, where `<abbr>` keys off each integral's *own* integrand
  (`integrand_abbrev`: `ovl`/`kin`/`eri`) and `L` is the built-side momentum.
- contraction of the VRR **base** integrals HRR consumes:
  `c<abbr>_<shells> = CDenseMatrix(...); osfunc::contract(c..., p...)`.
- the HRR transfer to the target (only when both sides carry momentum):
  `c<abbr>_<shells> = os2c::hrr::compute_<la>_<lb>(pair, <contracted base>)`. HRR
  is operator-independent, so its namespace carries no operator abbreviation.
- the Cartesian→spherical store: `osfunc::transform<la, lb>(buffer, <target>)`,
  then `return buffer`.

The `(s|s)` target is special-cased: the contracted primitive overlaps are
already the spherical result, so it emits `osfunc::contract(buffer, povl_ss)` and
returns, skipping the transform. `ObaraSaikaFunc.hpp` is always included;
`CartesianToSphericalFunc.hpp` only when a transform is actually emitted.

The `osfunc::*`, `os2c::vrr::*`, and `os2c::hrr::*` functions the body calls are
**not themselves generated yet** — emitting those primitive recurrence kernels is
the next plug-in point. Note also that the VRR/HRR/transform workflow currently
only fires for primitive-overlap-seeded operators; `electron_repulsion` (Boys-
seeded) falls through to an allocate-and-return stub for now. The
namespace/file/operator tags are derived freshly in the emitter (not via the
legacy `t2c::` helpers).

**Signature/storage awareness.** The parameter list is the cross-product of two
dimensions, each behind a no-`default` `switch`: `signature`
(`VeloxChemScreened`) selects the **inputs** — `const CScreenedBasisFunctionPair&
pair` — and `storage_form` (`VeloxChemSparse`) selects the **return type** —
`CDenseMatrix`. The matching VeloxChem headers (`ScreenedBasisFunctionPair.hpp`,
`DenseMatrix.hpp`) are included from the chosen types. The emitter is covered by
`tests/generators/test_two_center_emitters.cpp` (a new `generator_tests` target;
generator tests are globbed like the recursion tests). The tests chdir into a temp
directory because the emitter writes relative to the working directory.

## Conventions & pitfalls (read before editing)

- **`operator<` is a strict weak ordering.** A historical bug returned the wrong
  field on a tie; types in `std::set`/`std::map` silently misbehave if this is
  broken.
- **Guard your `std::optional`s.** `shift`, `shift_order`, `shift_prefix`,
  `shift_operator`, `bra_vrr`, etc. return `std::optional`; dereference only
  inside an `if (const auto x = ...)`. Several real bugs were unguarded `*opt`.
- **`replace()` is `const` and returns a new value.** `integral.replace(op);` on
  its own is a no-op — a live bug once shipped this way. Use the return value.
- **`for (const auto axis : "xyz")` also iterates the trailing `'\0'`.** It is
  safe today only because `shift('\0', …)` returns `nullopt` *before* any
  `_rxyz[axes::to_index('\0')]` indexing (which would be `_rxyz[-1]`). If you
  move an `_rxyz[...]` lookup above that guard you get an out-of-bounds read.
  Prefer `{'x','y','z'}` in new code.
- **libstdc++ (GCC/Linux) is stricter than libc++ (macOS).** Headers like
  `<cmath>` (for `std::floor`) are not always transitively included on Linux.
  If you call `std::floor/max/sort/find`, include the header explicitly — macOS
  builds can hide the omission; CI's ubuntu job will catch it.
- **Declared-but-undefined methods are legal until called.** A few drivers have
  historically carried dead declarations (link error only if invoked). Don't add
  a declaration without a definition.

## Testing approach

Tests are the safety net for refactors and for the recurrence math. Patterns
used throughout `tests/recursions/`:

- Build a term/integral with small angular momenta, call the driver method,
  assert the **term count** and the **set of `Factor` names**, plus the
  predicate's true/false branches and a rejection (wrong operator) case.
- A small helper collects factor names:
  ```cpp
  std::set<std::string> factor_names(const R2CDist& d) {
      std::set<std::string> n;
      for (const auto& f : d.unique_factors()) n.insert(f.name());
      return n;
  }
  ```
- Derive expected counts by **tracing the driver code** (each reached
  `t.add(...)` is one term; `nullopt` shifts are skipped). When a count is hard
  to pin down, prefer a robust assertion (`EXPECT_GE`, or `contains()` on a
  reduced base) over a wrong exact number.
- Enumerator (`v*i_*`) tests assert produced `std::set` sizes / membership; note
  that reductions often bottom out at a plain-overlap `"1"` auxiliary, so "every
  integral has operator X" is usually wrong.

Adding a test: drop `tests/recursions/test_<driver>.cpp`, reconfigure (the glob
picks it up), build, run. The fastest loop when an exact count is uncertain is to
write the test, run it, and read the actual value from the failure output.

## CI

`.github/workflows/ci.yml` builds and runs the full suite on **ubuntu-latest**
(GCC/libstdc++) and **macos-latest** (Apple Clang/libc++) on every push to `main`
and on PRs. The ubuntu job is the one that catches missing standard includes.

## Where to start reading

1. `src/algebra/tensor_component.{hpp,cpp}` and `integral_component.hpp` — the
   vocabulary.
2. `src/recursions/t2c_ovl_driver.cpp` — the simplest complete driver.
3. `tests/recursions/test_t2c_ovl_driver.cpp` — how a driver is exercised.
4. `src/recursions/t2c_defs.hpp` — the type aliases everything else uses.
5. `src/litmus.cpp` — the top-level wiring into the generators.
6. `src/generators/two_center_generators.cpp` — the new-style generator and the
   HRR/VRR integral-group split.
7. `src/generators/two_center_emitters.cpp` — the C++/CPU emitter that turns those
   groups into a `.hpp`/`.cpp` kernel pair, wiring the VRR/HRR/transform workflow
   as `osfunc::`/`os2c::vrr`/`os2c::hrr` calls (the active work; generating those
   primitive recurrence kernels themselves is the next step).
