# Backend Overrides

This page describes how backend library paths are selected and how users can override them.

## Resolution Order

At module initialization, backend libraries are resolved in this order:

1. Package artifacts from `Artifacts.toml` (default user path)
2. Environment variables
3. Local generated config file `deps/library_config.jl` (developer fallback)

At runtime, explicit API calls always take precedence over automatic initialization:

- `set_backend_library!("/path/to/lib"; precision=:double)`
- `set_backend_library!("/path/to/lib"; precision=:quad)`

## Environment Variables

Set one or both variables before starting Julia:

- `SPHEROIDALWAVES_LIBRARY_DOUBLE` for double precision backend
- `SPHEROIDALWAVES_LIBRARY_QUAD` for quad precision backend

Windows PowerShell example:

```powershell
$env:SPHEROIDALWAVES_LIBRARY_DOUBLE = "C:\path\to\spheroidal_batch_double.dll"
$env:SPHEROIDALWAVES_LIBRARY_QUAD = "C:\path\to\spheroidal_batch_quad.dll"
julia
```

Linux/macOS shell example:

```bash
export SPHEROIDALWAVES_LIBRARY_DOUBLE="/path/to/libspheroidal_batch_double.so"
export SPHEROIDALWAVES_LIBRARY_QUAD="/path/to/libspheroidal_batch_quad.so"
julia
```

## Programmatic Override

Use explicit runtime configuration when paths are known inside application code:

```julia
using SpheroidalWaves

set_backend_library!("/path/to/libspheroidal_batch_double.so"; precision=:double)
set_backend_library!("/path/to/libspheroidal_batch_quad.so"; precision=:quad)
```

## Inspect Active Paths

```julia
backend_library(precision=:double)
backend_library(precision=:quad)
```

## Notes

- Artifact entries must be populated for true out-of-box backend loading on fresh systems.
- If artifacts are missing for a platform, package import attempts a one-time local build automatically.
- `deps/library_config.jl` is generated locally by build tooling and is intentionally gitignored.
- If a configured path does not exist, it is ignored and a warning is emitted.
- If no backend is configured for a requested precision, calls fail with a clear error message.

