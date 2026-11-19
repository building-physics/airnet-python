# airnet contributor guidance

## Mission and scope

`airnet` is an object-oriented Python implementation of the AIRNET bulk
airflow model. The goal is a safe, understandable, testable implementation
that reproduces AIRNET's single-shot calculations and is convenient for
experimentation. Numerical behavior, file-format compatibility, and useful
Python APIs matter more than reproducing the legacy program's structure.

Prioritize, in order:

1. Correctly parse and represent AIRNET network inputs.
2. Produce converged node pressures and link flows that agree with trustworthy
   AIRNET results within documented numerical tolerances.
3. Detect invalid models, numerical failure, and non-convergence explicitly;
   never return plausible-looking partial results as success.
4. Keep the domain model object-oriented and usable independently of the CLI.
5. Keep documentation synchronized with every code change, including public
   behavior, APIs, inputs, outputs, limitations, and tested support.
6. Make focused changes backed by unit and end-to-end parity tests.

AIRNET's interactive model editing, repeated runs driven by stdin, terminal
prompts, and byte-for-byte report formatting are out of scope. Do not add new
physical models merely because they are useful elsewhere: this package should
implement behavior supported by AIRNET. A cleaner API, better diagnostics,
and stronger validation are welcome when they do not silently change the
model.

This is engineering calculation software, but it is not independently
validated building-safety software. State the implemented and tested scope
accurately in documentation and releases.

## Sources of truth

The adjacent `../AIRNET-C` repository is the compatibility reference. Use the
evidence in this order:

1. Historical reports in `../AIRNET-C/c/tests/reference/` and their archived
   originals in `../AIRNET-C/c/tests/legacy/RESULTS.CMB`.
2. Equations, definitions, units, and conventions in the NIST documentation
   under `../AIRNET-C/doc/`.
3. Results from the working C executable and the recovered C implementation
   under `../AIRNET-C/c/src/`.
4. The F77 sources under `../AIRNET-C/f77/` as supporting evidence only.

The C implementation was recovered from incomplete sources, and the F77
version is not known to have worked. Neither source tree is an unquestioned
oracle. When sources disagree, preserve the evidence in a regression test or
an explanatory comment and document the decision. Do not modify historical
reference reports.

The C build already present in `../AIRNET-C/out/cmake/` may be run when a new
comparison is needed. Do not rebuild AIRNET-C as a routine Python development
step, and do not spend time repairing F77 unless a task explicitly requests
it. C cases are stdin-driven and write `afout.*`; follow
`../AIRNET-C/AGENTS.md` and run them only in an isolated scratch directory.

## Repository map

- `src/airnet/reader.py`: legacy network-file parsing and input diagnostics.
- `src/airnet/afedat.py`: airflow-element objects and their constitutive
  relations/Jacobians.
- `src/airnet/model.py`: nodes, links, sparse system assembly, initialization,
  nonlinear solution, result export, and command entry points.
- `tests/`: fast parser, model, solver, and CLI unit tests.
- `regression-tests/`: end-to-end numerical comparisons.
- `examples/`: user-facing model inputs; examples are not substitutes for
  assertions against reference results.
- `docs/`: user and developer documentation.

Generated CSV files, C reports, and temporary case directories do not belong
in source or fixture directories. Put them in a temporary directory or an
ignored build area.

## What parity means

Treat support for a feature or historical case as a progression:

1. **Parse:** all required records are read without losing values.
2. **Construct:** typed Python objects are created and cross-references are
   validated.
3. **Solve:** the element participates correctly in initialization, residual
   and Jacobian assembly, iteration, and convergence checks.
4. **Match:** pressures, pressure drops, both flow components where applicable,
   flow direction, and convergence behavior agree with reference evidence.

Do not call an element or case supported at an earlier step. A parser test is
not a simulation test, and a successful return is not numerical parity.

Compare structured numerical results rather than report text. Use explicit
absolute-plus-relative tolerances appropriate to each quantity. Also check:

- mass balance at every variable-pressure node;
- expected flow direction and one-way/two-way flow count;
- finite values and a genuinely satisfied convergence criterion;
- units, density/viscosity treatment, stack and wind contributions;
- behavior on reversed links and at zero pressure drop.

Formatting differences, signed zero, sparse-solver iteration details, and
interactive transcript text are not parity failures unless they expose a
physical or numerical difference.

## Development workflow

Before changing code:

- Read the smallest historical input, run script, and reference report that
  exercise the behavior.
- Trace the corresponding C equations and data fields; consult the NIST report
  when names or conventions are unclear.
- Record the expected units, signs, flow count, edge cases, and comparison
  tolerances.
- Inspect `git status` and preserve unrelated user changes.

While changing code:

- Treat documentation as part of the implementation, not as release cleanup.
  Every code change requires a documentation-impact review, and all affected
  narrative documentation, docstrings, README material, CLI examples, and
  support-matrix entries must be updated in the same change. A code change is
  not complete while any documentation describes the previous behavior.
- Keep parsing, domain objects, physical equations, solver assembly, and
  presentation separate.
- Do not mutate caller-owned input records while constructing a `Model`.
- Prefer explicit typed state and meaningful exceptions over assertions for
  user/model errors.
- Keep element equations local to element classes. The model should assemble
  their residual/Jacobian contributions without type-specific special cases
  unless AIRNET's equations require one.
- Preserve AIRNET units, sign conventions, link orientation, and two-flow
  semantics. Add comments with source/provenance when a legacy choice is not
  obvious.
- Never hide non-convergence by returning the iteration limit as though it
  were success.
- Do not add a dependency when a small, well-tested implementation using the
  existing stack is sufficient.

## Tests and checks

Install and test through the project environment so SciPy and the package are
available:

```powershell
hatch run pytest tests regression-tests
hatch run ruff check .
hatch run docs:build
```

For a focused change, run the smallest relevant test first, then nearby cases,
then the full suite. Parser changes should exercise malformed input and all
affected historical files. Element/solver changes require direct equation and
Jacobian tests plus at least one end-to-end reference comparison. Check a
Jacobian against finite differences where practical. Build both documentation
formats after any documentation or public-interface change; warning-free HTML
and PDF builds are part of completion.

Do not weaken tolerances, overwrite expected files, or copy current Python
output into a fixture merely to make a test pass. Explain any intentional
deviation from historical results in the test and documentation.

Minimum validation for numerical work:

1. Unit tests for the changed parser/domain/element behavior.
2. The smallest historical single-shot scenario that exercises it.
3. Node mass balance, finite results, and convergence assertions.
4. Neighboring element/case tests to catch sign and orientation errors.
5. The complete Python suite once the focused checks pass.

Report exactly what was run, what passed, and what remains unsupported or
unverified. Do not describe the Python implementation as matching AIRNET based
only on successful parsing, construction, or a build.
