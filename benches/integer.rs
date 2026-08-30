//! Benchmarks for `math.integer`, the big-integer functions.
//!
//! Requires the `num-bigint` feature.

use divan::{Bencher, black_box};
use num_bigint::BigInt;
use pymath::math::integer;

fn main() {
    divan::main();
}

/// A deterministic big integer with roughly `bits` significant bits.
fn big(bits: u32) -> BigInt {
    let mut n = BigInt::from(1u32);
    let mut i = 0;
    while i < bits {
        n = (n << 13) + BigInt::from(6_700_417u32);
        i += 13;
    }
    n
}

const FACTORIAL_ARGS: [i64; 4] = [20, 200, 2_000, 20_000];

#[divan::bench(args = FACTORIAL_ARGS)]
fn factorial(bencher: Bencher, n: i64) {
    bencher.bench(|| integer::factorial(black_box(n)));
}

const COMB_ARGS: [(i64, i64); 3] = [(30, 15), (500, 250), (5_000, 2_500)];

#[divan::bench(args = COMB_ARGS)]
fn comb(bencher: Bencher, (n, k): (i64, i64)) {
    bencher.bench(|| integer::comb(black_box(n), black_box(k)));
}

#[divan::bench(args = COMB_ARGS)]
fn perm(bencher: Bencher, (n, k): (i64, i64)) {
    bencher.bench(|| integer::perm(black_box(n), Some(black_box(k))));
}

const BITS: [u32; 3] = [64, 1_024, 8_192];

#[divan::bench(consts = BITS)]
fn isqrt<const B: u32>(bencher: Bencher) {
    let n = big(B);
    bencher.bench(|| integer::isqrt(black_box(&n)));
}

#[divan::bench(consts = BITS)]
fn gcd<const B: u32>(bencher: Bencher) {
    let a = big(B);
    let b = big(B) * BigInt::from(7u32) + BigInt::from(1u32);
    let args = [&a, &b];
    bencher.bench(|| integer::gcd(black_box(&args)));
}

#[divan::bench(consts = BITS)]
fn lcm<const B: u32>(bencher: Bencher) {
    let a = big(B);
    let b = big(B) * BigInt::from(7u32) + BigInt::from(1u32);
    let args = [&a, &b];
    bencher.bench(|| integer::lcm(black_box(&args)));
}

#[divan::bench(consts = BITS)]
fn log_bigint<const B: u32>(bencher: Bencher) {
    let n = big(B);
    bencher.bench(|| pymath::math::log_bigint(black_box(&n), None));
}
