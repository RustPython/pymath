pymath
======

A binary representation compatible Rust implementation of Python's math library.

## Overview

`pymath` is a line-by-line port of the CPython math library to Rust. This library aims to provide mathematical functions that produce identical results to Python's standard math library at the bit level.

## Implementation

Each function has been carefully translated from CPython's C implementation to Rust, preserving the same algorithms, constants, and corner case handling. The code maintains the same numerical properties, but in Rust!

## Module Structure

- `pymath::math` - Real number math functions (Python's `math` module)
- `pymath::cmath` - Complex number functions (Python's `cmath` module, requires `complex` feature)
- `pymath::m` - Direct libm bindings

## Usage

```rust
use pymath::math::{gamma, lgamma};

fn main() {
    // Get the same results as Python's math.gamma and math.lgamma
    let gamma_result = gamma(4.5).unwrap();
    let lgamma_result = lgamma(4.5).unwrap();
}
```
