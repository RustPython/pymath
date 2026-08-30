//! Benchmarks for the complex-valued `cmath` module.

use divan::{Bencher, black_box};
use num_complex::Complex64;

fn main() {
    divan::main();
}

/// Number of inputs evaluated per iteration.
const N: usize = 128;

/// A spiral of complex values covering all four quadrants.
fn spiral(scale: f64) -> Vec<Complex64> {
    (0..N)
        .map(|i| {
            let t = i as f64 / N as f64 * std::f64::consts::TAU;
            let r = scale * (0.25 + i as f64 / N as f64);
            Complex64::new(r * t.cos(), r * t.sin())
        })
        .collect()
}

/// Generate one benchmark per `Complex64 -> Result<Complex64>` function.
macro_rules! bench_1 {
    ($($name:ident($scale:expr)),* $(,)?) => {
        $(
            #[divan::bench]
            fn $name(bencher: Bencher) {
                let zs = spiral($scale);
                bencher.bench(|| {
                    let mut acc = Complex64::new(0.0, 0.0);
                    for &z in &zs {
                        acc += pymath::cmath::$name(black_box(z))
                            .unwrap_or(Complex64::new(0.0, 0.0));
                    }
                    acc
                });
            }
        )*
    };
}

bench_1! {
    exp(2.0),
    log10(10.0),
    sqrt(10.0),
    sin(2.0),
    cos(2.0),
    tan(2.0),
    sinh(2.0),
    cosh(2.0),
    tanh(2.0),
    asin(0.8),
    acos(0.8),
    atan(2.0),
    asinh(2.0),
    acosh(2.0),
    atanh(0.8),
}

#[divan::bench]
fn log_natural(bencher: Bencher) {
    let zs = spiral(10.0);
    bencher.bench(|| {
        let mut acc = Complex64::new(0.0, 0.0);
        for &z in &zs {
            acc += pymath::cmath::log(black_box(z), None).unwrap_or(Complex64::new(0.0, 0.0));
        }
        acc
    });
}

#[divan::bench]
fn log_base(bencher: Bencher) {
    let zs = spiral(10.0);
    let base = Complex64::new(3.0, 0.0);
    bencher.bench(|| {
        let mut acc = Complex64::new(0.0, 0.0);
        for &z in &zs {
            acc += pymath::cmath::log(black_box(z), Some(black_box(base)))
                .unwrap_or(Complex64::new(0.0, 0.0));
        }
        acc
    });
}

#[divan::bench]
fn abs_phase(bencher: Bencher) {
    let zs = spiral(10.0);
    bencher.bench(|| {
        let mut acc = 0.0;
        for &z in &zs {
            let z = black_box(z);
            acc += pymath::cmath::abs(z) + pymath::cmath::phase(z).unwrap_or_default();
        }
        acc
    });
}

#[divan::bench]
fn polar_rect(bencher: Bencher) {
    let zs = spiral(10.0);
    bencher.bench(|| {
        let mut acc = Complex64::new(0.0, 0.0);
        for &z in &zs {
            let (r, phi) = pymath::cmath::polar(black_box(z)).unwrap_or((0.0, 0.0));
            acc += pymath::cmath::rect(r, phi).unwrap_or(Complex64::new(0.0, 0.0));
        }
        acc
    });
}
