# pymath Agent Guidelines

This project is a strict port of CPython's math and cmath modules to Rust.

## Core Principle

**Every function must match CPython exactly** - same logic, same special case handling, same error conditions.

- `math` module → `Modules/mathmodule.c`
- `cmath` module → `Modules/cmathmodule.c`

## Porting Rules

### 1. Use Existing Helpers

CPython uses helpers like `math_1`, `math_2`, `FUNC1`, `FUNC2`, etc. We have Rust equivalents:

| CPython | Rust |
|---------|------|
| `FUNC1(name, func, can_overflow, ...)` | `math_1(x, func, can_overflow)` |
| `FUNC2(name, func, ...)` | `math_2(x, y, func)` |
| `m_log`, `m_log10`, `m_log2` | Same names in `exponential.rs` |

If a function uses a helper in CPython, use the corresponding helper here.

### 2. Create Missing Helpers

If CPython has a helper we don't have yet, implement it. Examples:
- `math_1_fn` for Rust function pointers (vs C function pointers)
- Special case handlers for specific functions

### 3. Error Handling

CPython sets `errno` and calls `is_error()`. We return `Result` directly:

```rust
// CPython:
errno = EDOM;
if (errno && is_error(r, 1)) return NULL;

// Rust:
return Err(crate::Error::EDOM);
```

**Never use `set_errno(libc::EDOM)` or similar** - just return `Err()` directly.
The only valid `set_errno` call is `set_errno(0)` to clear errno before libm calls.

### 4. Special Cases

CPython has explicit special case handling for IEEE specials (NaN, Inf, etc.). Copy this logic exactly:

```rust
// Example from pow():
if !x.is_finite() || !y.is_finite() {
    if x.is_nan() {
        return Ok(if y == 0.0 { 1.0 } else { x }); // NaN**0 = 1
    }
    // ... more special cases
}
```

### 5. Reference CPython Source

Always check CPython source in the `cpython/` directory:

**For math module** (`Modules/mathmodule.c`):
- Function implementations
- Helper macros (`FUNC1`, `FUNC1D`, `FUNC2`)
- Special case comments
- Error conditions

**For cmath module** (`Modules/cmathmodule.c`):
- `special_type()` enum and function
- 7x7 special value tables (e.g., `tanh_special_values`)
- `SPECIAL_VALUE` macro usage
- Complex-specific error handling

### 6. Fused Multiply-Add (mul_add)

For bit-exact matching with CPython, use `crate::mul_add(a, b, c)` instead of `a * b + c` in specific cases.

**Why this matters**: CPython compiled with clang on macOS may use FMA (fused multiply-add) instructions for expressions like `1.0 + x * x`. FMA computes `a * b + c` in a single operation without intermediate rounding, which can produce results that differ by 1-2 ULP from separate multiply and add operations.

**When to use `mul_add`**:
- Expressions of the form `a * b + c` or `c + a * b` in complex math functions
- Especially in formulas like `1.0 + x * x` → `mul_add(x, x, 1.0)`

**Example** (from `c_tanh`):
```rust
// Wrong - may differ from CPython by 1-2 ULP
let denom = 1.0 + txty * txty;
let r_re = tx * (1.0 + ty * ty) / denom;

// Correct - matches CPython exactly
let denom = mul_add(txty, txty, 1.0);
let r_re = tx * mul_add(ty, ty, 1.0) / denom;
```

**Example** (from `c_asinh`):
```rust
// mul_add for cross-product calculations
let r_re = m::asinh(mul_add(s1.re, s2.im, -(s2.re * s1.im)));
let r_im = m::atan2(z.im, mul_add(s1.re, s2.re, -(s1.im * s2.im)));
```

**Example** (from `c_atanh`):
```rust
// mul_add for squared terms
let one_minus_re = 1.0 - z.re;
let r_re = m::log1p(4.0 * z.re / mul_add(one_minus_re, one_minus_re, ay * ay)) / 4.0;
let r_im = -m::atan2(-2.0 * z.im, mul_add(one_minus_re, 1.0 + z.re, -(ay * ay))) / 2.0;
```

**Feature flag**: The `mul_add` feature controls whether hardware FMA is used:

- `mul_add` enabled: Uses `f64::mul_add()` (hardware FMA instruction)
- `mul_add` disabled (default): Falls back to `a * b + c` (separate operations)

Note: macOS CI always enables `mul_add` because CPython on macOS uses FMA.

**How to identify missing mul_add usage**:
1. If a test fails with 1-2 ULP difference
2. Look for `a * b + c` or `c + a * b` patterns in the failing function
3. Replace with `mul_add(a, b, c)` and re-test

### 7. Platform-specific sincos (macOS, cmath only)

On macOS, Python's cmath module uses Apple's `__sincos_stret` function, which computes sin and cos together with slightly different results than calling them separately (up to 1 ULP difference).

For bit-exact matching on macOS, use `m::sincos(x)` which returns `(sin, cos)` tuple:

```rust
// Instead of:
let sin_x = m::sin(x);
let cos_x = m::cos(x);

// Use:
let (sin_x, cos_x) = m::sincos(x);
```

Required in cmath functions that use both sin and cos of the same angle:

- `cosh`, `sinh` - for the imaginary argument
- `exp` - for the imaginary argument
- `rect` - for the phi angle

On non-macOS platforms, `m::sincos(x)` falls back to calling sin and cos separately.

## Testing

### EDGE_VALUES

All float functions must be tested with `crate::test::EDGE_VALUES` which includes:
- Zeros: `0.0`, `-0.0`
- Infinities: `INFINITY`, `NEG_INFINITY`
- NaNs: `NAN`, `-NAN`, and NaN with different payload
- Subnormals
- Boundary values: `MIN_POSITIVE`, `MAX`, `MIN`
- Large values near infinity
- Trigonometric special values: `PI`, `PI/2`, `PI/4`, `TAU`

### Error Type Verification

Tests must verify both:
1. Correct values for Ok results
2. Correct error types (EDOM vs ERANGE) for Err results

Python `ValueError` → `Error::EDOM`
Python `OverflowError` → `Error::ERANGE`

## File Structure

### Core
- `src/lib.rs` - Root module, `mul_add` function
- `src/err.rs` - Error types (EDOM, ERANGE)
- `src/test.rs` - Test helpers, `EDGE_VALUES`, `EDGE_INTS`

### System libm bindings
- `src/m_sys.rs` - Raw FFI declarations (`extern "C"`)
- `src/m.rs` - Safe wrappers, platform-specific `sincos`

### math module
- `src/math.rs` - Main module, `math_1`, `math_2` helpers, `hypot`, constants
- `src/math/exponential.rs` - exp, log, pow, sqrt, cbrt, etc.
- `src/math/trigonometric.rs` - sin, cos, tan, asin, acos, atan, etc.
- `src/math/misc.rs` - frexp, ldexp, modf, fmod, copysign, isclose, ulp, etc.
- `src/math/gamma.rs` - gamma, lgamma, erf, erfc
- `src/math/aggregate.rs` - fsum, prod, sumprod, dist (vector operations)
- `src/math/integer.rs` - gcd, lcm, isqrt, comb, perm, factorial (requires `_bigint` feature)

### cmath module (requires `complex` feature)
- `src/cmath.rs` - Main module, `special_type`, `special_value!` macro, shared constants
- `src/cmath/exponential.rs` - sqrt, exp, log, log10
- `src/cmath/trigonometric.rs` - sin, cos, tan, sinh, cosh, tanh, asin, acos, atan, asinh, acosh, atanh
- `src/cmath/misc.rs` - phase, polar, rect, abs, isfinite, isnan, isinf, isclose
