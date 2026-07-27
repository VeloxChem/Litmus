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

Rough sizes: `algebra` ~29 files, `recursions` ~111 files (54 drivers + headers
+ defs), `generators` ~152 files, plus a test suite of **518 tests**.

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
- **`sphar::`** (`spherical_harmonics.{hpp,cpp}`) — exact Cartesian-to-spherical
  transformation coefficients. A `SphericalFactor` is a `Fraction` times
  `sqrt(square-free radicand)`, the closed form every real solid-harmonic
  coefficient takes in the Helgaker–Jørgensen–Olsen convention; `operator*`
  reduces the combined radical (`sqrt(3) * sqrt(3) -> 3`). `spherical_factors(l,
  m)` / `spherical_component_factors(l, component)` give one shell's expansion,
  and `two_center_spherical_factors(la, lb, ...)` their tensor product — the
  two-center transform, which is **operator-independent**. This is what lets the
  emitters fuse the transform into the recurrence kernels instead of emitting a
  separate transform pass.

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
ctest --test-dir build --output-on-failure   # run all 518 tests
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

**Legacy schema** (14 generator families). Keys: `type` (required),
`lmax`, `integral`, `geom` (arity 3/4/5 per family), `aux_lmax` (t3c),
`proj_lmax` (proj-ecp), `rec_form` (t2c, `[1, 0]`), `use_rs` (t2c/g2c). The
`spherical_momentum` family is the odd one out: it emits VeloxChem's
`SphericalMomentum.hpp` from `sphar::` and reads `lmax` as the highest shell to
tabulate.

**New-style schema** (`cfg::RunConfiguration` in
`src/general/run_configuration.{hpp,cpp}`) decomposes the monolithic `type` into
orthogonal, typed dimensions for the next generation of generators. Keys:
`integral_type` (`two_center`/`three_center`/`four_center`, also `2c`/`3c`/`4c`)
**or** `recursion_type` (`hrr_bra_ket`/`hrr_bra`/`hrr_ket` for the horizontal-
recurrence transfer, `vrr_cartesian`/`vrr_spherical` for the vertical one — i.e.
which *recurrence kernel* to emit rather than a whole integral) — **exactly one
of the two is required**, and
supplying both is a `ConfigError`; they are stored as `std::optional` fields so a
consumer tests which one is set. `max_ang_mom` (required), `min_ang_mom`
(default 0),
`operator_type` (default `overlap`; the integrand — `overlap`, `kinetic_energy`,
`nuclear_potential`, `electron_repulsion`, `dipole_momentum`, `linear_momentum`,
`local_ecp`, `projected_ecp`, `three_center_overlap`, `three_center_r2`,
`three_center_r_dot_r2`, with `to_string` round-tripping to the generator label),
`hardware` (default `cpu`), `language` (default `C++`), `storage_form` (default
`VeloxChemSparse`), `signature` (default `VeloxChemScreened`). Each enumerated
field is validated against its allowed spellings (case/`_`/`-` insensitive) and
the angular-momentum range is checked. `litmus run` recognizes a config as
new-style when it carries an `integral_type` *or* `recursion_type` key, and
dispatches it to one of the three generators described below. See also
`examples/four_center.toml`.

## The new-style two-center generators (work in progress)

Three generators are built on `cfg::RunConfiguration`. Each writes its
`.hpp`/`.cpp` output into the current working directory and prints one line per
kernel to stdout:

| config key | generator | emits |
| --- | --- | --- |
| `integral_type = two_center` | `TwoCenterGenerator` | orchestration kernels — `os2c::<op>::compute_<la>_<lb>` |
| `recursion_type = hrr_bra_ket`/`hrr_bra`/`hrr_ket` | `TwoCenterHrrGenerator` | `os2c::hrr::compute_<la>_<lb>` |
| `recursion_type = vrr_cartesian` | `TwoCenterVrrGenerator` | `os2c::vrr::ovl::compute_<l>` |
| `recursion_type = vrr_spherical` | `TwoCenterVrrGenerator` | `os2c::ovl::compute_<l>_sph` |

`integral_type = three_center`/`four_center` still print a "not wired in yet"
notice. The enumerated `switch`es throughout deliberately carry no `default`, so
adding an `OperatorType` / `RecursionType` / `Hardware` / … enumerator trips
`-Wswitch` at every dispatch site that needs updating.

### The osfunc ABI

Everything these three emit speaks VeloxChem's `osfunc::` data structures. This is
the single most important thing to know before editing an emitter:

- **All data is `osfunc::CArray<double>`** — an `(nrows × npairs)` buffer with one
  screened atom pair per column. `.row(i)` yields the pointer the SIMD loops walk;
  `.nrows()`/`.ncols()` give the extents. Contracted Cartesian blocks pack their
  components into the rows; primitive blocks pack `nprims * ncomponents` rows.
- **The pair is `osfunc::CBasisFunctionPair`** — `.bra()`, `.ket()`,
  `.number_of_pairs()`, `.number_of_primitive_pairs()`.
- **Recurrence kernels return `void`** and write through a trailing out-parameter
  (`osfunc::CArray<double>& out`). Only the orchestration kernel returns a value,
  `-> osfunc::CArray<double>`.
- **There is no standalone Cartesian→spherical transform call.** The transform is
  *fused* into the HRR kernel and into the spherical VRR kernel using the exact
  `sphar::` coefficients; nothing emits an `osfunc::transform` any more.
- Generated headers include `Array.hpp`, plus `BasisFunctionPair.hpp` when the
  kernel takes a pair. The orchestration `.cpp` also includes `ObaraSaikaFunc.hpp`
  and the header of every `os2c` kernel it calls.

### The orchestration generator (`integral_type = two_center`)

`TwoCenterGenerator` (`src/generators/two_center_generators.{hpp,cpp}`) supports
the operators `overlap`, `kinetic_energy`, and `electron_repulsion`; every other
`OperatorType` raises a `cfg::ConfigError`.

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
`CppCpuTwoCenterEmitter` (C++ on CPU). `generate()` prints a one-line summary per
target to stdout (`Generated PP kernel (1 HRR, 2 VRR base, 1 VRR rest)`). Try:

```bash
printf 'integral_type="two_center"\nmax_ang_mom=2\n' > /tmp/t2c.toml
cd /tmp && /path/to/litmus.x run /tmp/t2c.toml   # writes ObaraSaikaTwoCenter*.{hpp,cpp} into the cwd
```

**What the emitter produces.** Files are named
`ObaraSaikaTwoCenter<Operator><BraKet>.{hpp,cpp}` (e.g.
`ObaraSaikaTwoCenterOverlapPP`). Each `.hpp` is a guarded declaration in the
operator's namespace — `os2c::ovl` (overlap), `os2c::kin` (kinetic energy),
`os2c::eri` (electron repulsion) — of a `compute_<la>_<lb>` kernel (lowercase
shell labels, e.g. `compute_p_p`); the `(bra|op|ket)^order` caption is used for
the `@brief` doc line. The `.cpp` carries the matching definition. Every body
opens with `npairs = pair.number_of_pairs()` and a zeroed spherical result buffer
`osfunc::CArray<double> buffer((2*la+1)*(2*lb+1), npairs)`, and closes with
`return buffer`. In between it takes one of four shapes:

- **a non-overlap operator** — nothing but a `TODO` comment. Only the overlap
  recurrence kernels are generated so far, so `kinetic_energy` and
  `electron_repulsion` return the zeroed buffer.
- **`(s|s)`** — the contracted primitive overlaps *are* the spherical result:
  `osfunc::contract(buffer, osfunc::compute_overlap(pair));`.
- **`(s|lb)` or `(la|s)`** — one fused spherical VRR call. The seed is
  `ss = osfunc::compute_overlap(pair)`; the distances are `osfunc::compute_pb(pair)`
  when the recurrence grows the ket (`la <= lb`) and `compute_pa` when it grows the
  bra. Then `os2c::ovl::compute_<l>_sph(pair, ss, pb, buffer)` does the recurrence,
  the contraction and the transform in a single pass.
- **`(la|lb)` with both sides > 0** — the full ladder. After
  `nprims = pair.number_of_primitive_pairs()`, one
  `os2c::vrr::ovl::compute_<l>(pair, <two lower>, pb, <out>)` per VRR integral
  (lowest total momentum first, into `nprims * ncart` row buffers), then
  `osfunc::contract` of each base integral the HRR consumes, then
  `ab = osfunc::compute_ab(pair)` and
  `os2c::hrr::compute_<la>_<lb>(<contracted base…>, ab, buffer)`, which transfers
  the momentum with the transform fused. HRR is operator-independent, so its
  namespace carries no operator abbreviation.

The tail of a generated `(d|d)` overlap kernel, for orientation:

```cpp
osfunc::CArray<double> sg(nprims * 15, npairs);
os2c::vrr::ovl::compute_g(pair, sf, sd, pb, sg);

osfunc::CArray<double> csg(15, npairs);
osfunc::contract(csg, sg);

const auto ab = osfunc::compute_ab(pair);
os2c::hrr::compute_d_d(csd, csf, csg, ab, buffer);
```

The namespace/file/operator tags are derived freshly in the emitter (not via the
legacy `t2c::` helpers).

**Signature/storage awareness.** The parameter list is the cross-product of two
dimensions, each behind a no-`default` `switch`: `signature`
(`VeloxChemScreened`) selects the **inputs** — `const osfunc::CBasisFunctionPair&
pair` — and `storage_form` (`VeloxChemSparse`) selects the **return type** —
`osfunc::CArray<double>`. The matching VeloxChem headers
(`BasisFunctionPair.hpp`, `Array.hpp`) are included from the chosen types.

### The HRR kernel generator (`recursion_type = hrr_*`)

`TwoCenterHrrGenerator` and `format_hrr_kernel(la, lb)`
(`two_center_hrr_generators.{hpp,cpp}`, `two_center_hrr_emitter.{hpp,cpp}`) emit
`ObaraSaikaTwoCenterHrr<AB>.{hpp,cpp}` in `os2c::hrr`. The emitter runs
`T2CHRRDriver` to reduce the bra (or ket) fully down to `s`, folds in the exact
`sphar::two_center_spherical_factors` coefficients, and emits one SIMD loop per
spherical component over the atom-pair columns. Irrational and non-terminating
factors (`sqrt(3)`, `1/3`) are hoisted into `const double` locals, and terms are
grouped by AB product with the common rational coefficient factored out.

The kernel takes one contracted Cartesian `CArray` per base integral, the AB
distances, and the spherical target as an out-parameter; it loops over
`out.nrows() / ncomps` integral blocks, so one call transfers a whole batch.
Both sides must carry momentum (`la, lb >= 1`): `hrr_bra` keeps the
bra-incremented targets (`la <= lb`), `hrr_ket` the ket-incremented ones
(`la > lb`), and `hrr_bra_ket` all of them.

### The VRR kernel generator (`recursion_type = vrr_*`)

`TwoCenterVrrGenerator` and `two_center_vrr_emitter.{hpp,cpp}` emit
`ObaraSaikaTwoCenterOverlapVrr{Cart,Sph}<L>.{hpp,cpp}`, one pair per ket shell
`lb >= 1`; a `VrrFlavor` table selects the namespace, file tag and emitter
functions per recursion type. Both flavors build `(s|lb)` from the `(s|s)` seed:

- **`vrr_cartesian`** → `os2c::vrr::ovl::compute_<l>(pair, <two lower>, pc, out)`:
  a *single* Obara–Saika ket step, per-primitive, retaining the `fe = 1/(2 eta)`
  term. Its output feeds the HRR and is contracted downstream.
- **`vrr_spherical`** → `os2c::ovl::compute_<l>_sph(pair, ss, pc, out)`: the fully
  reduced recurrence with the contraction *and* the Cartesian→spherical transform
  folded in. Worth knowing why this is cheap: the traceless real solid harmonics
  cancel every `fe` (trace) term, so the spherical kernels come out as pure `pc`
  polynomials — no `fe`, no exponents (verified through `l = 6`). They GCD-factor
  the common rational and radical and wrap long inner sums one term per line.

All three emitters are covered by `tests/generators/` —
`test_two_center_emitters.cpp`, `test_two_center_hrr_emitter.cpp`,
`test_two_center_vrr_emitter.cpp` (target `generator_tests`; generator tests are
globbed like the recursion tests). The tests chdir into a temp directory because
the generators write relative to the working directory.

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
   groups into a `.hpp`/`.cpp` orchestration kernel, wiring the VRR/HRR workflow
   as `osfunc::`/`os2c::vrr`/`os2c::hrr` calls over the `CArray` ABI.
8. `src/generators/two_center_vrr_emitter.cpp` and
   `two_center_hrr_emitter.cpp` — the recurrence kernels themselves, where the
   symbolic recursion and the `sphar::` transform coefficients are turned into
   SIMD loops (the active work).

**Current frontier.** The overlap VRR and HRR kernels are generated and compile
clean against VeloxChem's `osfunc::CArray`. Still open: `kinetic_energy` and
`electron_repulsion` are accepted by the config and by the orchestration
generator's group split, but their recurrence kernels do not exist, so the
orchestration emitter falls through to a zeroed-buffer stub for them — ERI
additionally needs a Boys-function seed rather than a primitive overlap. The
three- and four-center `integral_type` values validate but have no generators.
