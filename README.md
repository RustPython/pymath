pymath
======

**0 ULP (bit-exact) compatibility with CPython's math and cmath modules.**

Every function produces identical results to Python at the binary representation level - not just "close enough", but exactly the same bits.

## Overview

`pymath` is a strict port of CPython's math library to Rust. Each function has been carefully translated from CPython's C implementation, preserving the same algorithms, constants, and corner case handling.

## Compatibility Status

### math (55/55)

- [x] `ceil`
- [x] `copysign`
- [x] `fabs`
- [x] `floor`
- [x] `fmod`
- [x] `frexp`
- [x] `isclose`
- [x] `isfinite`
- [x] `isinf`
- [x] `isnan`
- [x] `ldexp`
- [x] `modf`
- [x] `nextafter`
- [x] `remainder`
- [x] `trunc`
- [x] `ulp`
- [x] `cbrt`
- [x] `exp`
- [x] `exp2`
- [x] `expm1`
- [x] `log`
- [x] `log10`
- [x] `log1p`
- [x] `log2`
- [x] `pow`
- [x] `sqrt`
- [x] `acos`
- [x] `acosh`
- [x] `asin`
- [x] `asinh`
- [x] `atan`
- [x] `atan2`
- [x] `atanh`
- [x] `cos`
- [x] `cosh`
- [x] `sin`
- [x] `sinh`
- [x] `tan`
- [x] `tanh`
- [x] `erf`
- [x] `erfc`
- [x] `gamma`
- [x] `lgamma`
- [x] `dist`
- [x] `fsum`
- [x] `hypot`
- [x] `prod`
- [x] `sumprod`
- [x] `degrees`
- [x] `radians`
- [x] `pi`
- [x] `e`
- [x] `tau`
- [x] `inf`
- [x] `nan`

### math.integer (6/6, requires `num-bigint` or `malachite-bigint` feature)

- [x] `comb`
- [x] `factorial`
- [x] `gcd`
- [x] `isqrt`
- [x] `lcm`
- [x] `perm`

### cmath (24/24, requires `complex` feature)

- [x] `abs`
- [x] `isclose`
- [x] `isfinite`
- [x] `isinf`
- [x] `isnan`
- [x] `phase`
- [x] `polar`
- [x] `rect`
- [x] `exp`
- [x] `log`
- [x] `log10`
- [x] `sqrt`
- [x] `acos`
- [x] `acosh`
- [x] `asin`
- [x] `asinh`
- [x] `atan`
- [x] `atanh`
- [x] `cos`
- [x] `cosh`
- [x] `sin`
- [x] `sinh`
- [x] `tan`
- [x] `tanh`

## Usage

```rust
use pymath::math::{sqrt, gamma, lgamma};

fn main() {
    let sqrt_result = sqrt(2.0).unwrap();
    let gamma_result = gamma(4.5).unwrap();
    let lgamma_result = lgamma(4.5).unwrap();
}
```

### Bit-exact verification

```python
# Python 3.14
>>> import math
>>> math.gamma(0.5).hex()
'0x1.c5bf891b4ef6ap+0'
```

```rust
// Rust - identical bits
assert_eq!(
    pymath::math::gamma(0.5).unwrap().to_bits(),
    0x3ffc5bf891b4ef6a
);
```

*Bit representation may vary across platforms, but CPython and pymath built on the same environment will always produce identical results.*

## Features

- `complex` (default) - Enable cmath module for complex number functions
- `num-bigint` - Enable integer functions using num-bigint
- `malachite-bigint` - Enable integer functions using malachite-bigint
- `mul_add` - Use hardware FMA for bit-exact macOS compatibility

## Module Structure

- `pymath::math` - Real number math functions (Python's `math` module)
- `pymath::cmath` - Complex number functions (Python's `cmath` module)
- `pymath::m` - Direct libm bindings
