pymath
======

A binary representation compatible Rust implementation of Python's math library.

## Overview

`pymath` is a line-by-line port of the CPython math library to Rust. This library aims to provide mathematical functions that produce identical results to Python's standard math library at the bit level.

## Implementation

Each function has been carefully translated from CPython's C implementation to Rust, preserving the same algorithms, constants, and corner case handling. The code maintains the same numerical properties, but in Rust!

## Usage

```rust
use pymath::{gamma, lgamma};

fn main() {
    // Get the same results as Python's math.gamma and math.lgamma
    let gamma_result = gamma(4.5).unwrap();
    let lgamma_result = lgamma(4.5).unwrap();
}
```
