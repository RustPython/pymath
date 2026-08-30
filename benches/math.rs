//! Benchmarks for the real-valued `math` module.
//!
//! Every benchmark evaluates the function over a fixed sweep of representative
//! finite inputs so that the measurement is dominated by the function itself
//! rather than by the harness.

use divan::{Bencher, black_box};

fn main() {
    divan::main();
}

/// Number of inputs evaluated per iteration.
const N: usize = 256;

/// Evenly spaced values over `[start, end)`.
fn sweep(start: f64, end: f64) -> Vec<f64> {
    let step = (end - start) / N as f64;
    (0..N).map(|i| start + step * i as f64).collect()
}

/// Generate one benchmark per `f64 -> Result<f64>` function.
macro_rules! bench_1 {
    ($module:path; $($name:ident($start:expr, $end:expr)),* $(,)?) => {
        $(
            #[divan::bench]
            fn $name(bencher: Bencher) {
                use $module as m;
                let xs = sweep($start, $end);
                bencher.bench(|| {
                    let mut acc = 0.0;
                    for &x in &xs {
                        acc += m::$name(black_box(x)).unwrap_or_default();
                    }
                    acc
                });
            }
        )*
    };
}

mod trigonometric {
    use super::*;

    bench_1! {
        pymath::math;
        sin(-10.0, 10.0),
        cos(-10.0, 10.0),
        tan(-1.5, 1.5),
        asin(-1.0, 1.0),
        acos(-1.0, 1.0),
        atan(-10.0, 10.0),
        sinh(-5.0, 5.0),
        cosh(-5.0, 5.0),
        tanh(-5.0, 5.0),
        asinh(-10.0, 10.0),
        acosh(1.0, 20.0),
        atanh(-0.999, 0.999),
    }

    #[divan::bench]
    fn atan2(bencher: Bencher) {
        let ys = sweep(-10.0, 10.0);
        let xs = sweep(-7.0, 13.0);
        bencher.bench(|| {
            let mut acc = 0.0;
            for (&y, &x) in ys.iter().zip(&xs) {
                acc += pymath::math::atan2(black_box(y), black_box(x)).unwrap_or_default();
            }
            acc
        });
    }
}

mod exponential {
    use super::*;

    bench_1! {
        pymath::math;
        exp(-20.0, 20.0),
        exp2(-20.0, 20.0),
        expm1(-1.0, 1.0),
        log10(1e-8, 1e8),
        log2(1e-8, 1e8),
        log1p(-0.5, 10.0),
        sqrt(0.0, 1e6),
        cbrt(-1e6, 1e6),
    }

    #[divan::bench]
    fn log_natural(bencher: Bencher) {
        let xs = sweep(1e-8, 1e8);
        bencher.bench(|| {
            let mut acc = 0.0;
            for &x in &xs {
                acc += pymath::math::log(black_box(x), None).unwrap_or_default();
            }
            acc
        });
    }

    #[divan::bench]
    fn log_base(bencher: Bencher) {
        let xs = sweep(1e-8, 1e8);
        bencher.bench(|| {
            let mut acc = 0.0;
            for &x in &xs {
                acc += pymath::math::log(black_box(x), Some(black_box(3.0))).unwrap_or_default();
            }
            acc
        });
    }

    #[divan::bench]
    fn pow(bencher: Bencher) {
        let xs = sweep(0.1, 20.0);
        let ys = sweep(-5.0, 5.0);
        bencher.bench(|| {
            let mut acc = 0.0;
            for (&x, &y) in xs.iter().zip(&ys) {
                acc += pymath::math::pow(black_box(x), black_box(y)).unwrap_or_default();
            }
            acc
        });
    }
}

mod gamma {
    use super::*;

    bench_1! {
        pymath::math;
        erf(-4.0, 4.0),
        erfc(-4.0, 4.0),
        gamma(0.5, 20.0),
        lgamma(0.5, 200.0),
    }

    /// `erf`/`erfc` switch algorithms around |x| = 1.5, so cover both branches.
    #[divan::bench]
    fn erf_tail(bencher: Bencher) {
        let xs = sweep(1.5, 6.0);
        bencher.bench(|| {
            let mut acc = 0.0;
            for &x in &xs {
                acc += pymath::math::erf(black_box(x)).unwrap_or_default();
            }
            acc
        });
    }

    #[divan::bench]
    fn gamma_negative(bencher: Bencher) {
        let xs = sweep(-19.5, -0.5);
        bencher.bench(|| {
            let mut acc = 0.0;
            for &x in &xs {
                acc += pymath::math::gamma(black_box(x)).unwrap_or_default();
            }
            acc
        });
    }
}

mod misc {
    use super::*;

    bench_1! {
        pymath::math;
        fabs(-1e6, 1e6),
    }

    #[divan::bench]
    fn floor_ceil_trunc(bencher: Bencher) {
        let xs = sweep(-1e6, 1e6);
        bencher.bench(|| {
            let mut acc = 0.0;
            for &x in &xs {
                let x = black_box(x);
                acc += pymath::math::floor(x) + pymath::math::ceil(x) + pymath::math::trunc(x);
            }
            acc
        });
    }

    #[divan::bench]
    fn frexp_modf(bencher: Bencher) {
        let xs = sweep(-1e6, 1e6);
        bencher.bench(|| {
            let mut acc = 0.0;
            for &x in &xs {
                let x = black_box(x);
                let (m, e) = pymath::math::frexp(x);
                let (fract, int) = pymath::math::modf(x);
                acc += m + e as f64 + fract + int;
            }
            acc
        });
    }

    #[divan::bench]
    fn fmod(bencher: Bencher) {
        let xs = sweep(-1e3, 1e3);
        let ys = sweep(0.5, 17.0);
        bencher.bench(|| {
            let mut acc = 0.0;
            for (&x, &y) in xs.iter().zip(&ys) {
                acc += pymath::math::fmod(black_box(x), black_box(y)).unwrap_or_default();
            }
            acc
        });
    }

    #[divan::bench]
    fn remainder(bencher: Bencher) {
        let xs = sweep(-1e3, 1e3);
        let ys = sweep(0.5, 17.0);
        bencher.bench(|| {
            let mut acc = 0.0;
            for (&x, &y) in xs.iter().zip(&ys) {
                acc += pymath::math::remainder(black_box(x), black_box(y)).unwrap_or_default();
            }
            acc
        });
    }

    #[divan::bench]
    fn fma(bencher: Bencher) {
        let xs = sweep(-1e3, 1e3);
        let ys = sweep(0.5, 17.0);
        let zs = sweep(-5.0, 5.0);
        bencher.bench(|| {
            let mut acc = 0.0;
            for ((&x, &y), &z) in xs.iter().zip(&ys).zip(&zs) {
                acc +=
                    pymath::math::fma(black_box(x), black_box(y), black_box(z)).unwrap_or_default();
            }
            acc
        });
    }

    #[divan::bench]
    fn nextafter_steps(bencher: Bencher) {
        let xs = sweep(-1e3, 1e3);
        bencher.bench(|| {
            let mut acc = 0.0;
            for &x in &xs {
                acc += pymath::math::nextafter(black_box(x), black_box(f64::INFINITY), Some(64));
            }
            acc
        });
    }

    #[divan::bench]
    fn ldexp_ulp(bencher: Bencher) {
        let xs = sweep(-1e3, 1e3);
        bencher.bench(|| {
            let mut acc = 0.0;
            for &x in &xs {
                let x = black_box(x);
                acc += pymath::math::ldexp(x, 12).unwrap_or_default() + pymath::math::ulp(x);
            }
            acc
        });
    }

    #[divan::bench]
    fn isclose(bencher: Bencher) {
        let xs = sweep(-1e3, 1e3);
        bencher.bench(|| {
            let mut count = 0usize;
            for &x in &xs {
                let x = black_box(x);
                if pymath::math::isclose(x, x + 1e-9, None, None).unwrap_or(false) {
                    count += 1;
                }
            }
            count
        });
    }
}

mod aggregate {
    use super::*;

    /// Sizes exercising both the small and the large paths of the reductions.
    const SIZES: [usize; 3] = [8, 128, 4096];

    fn values(len: usize) -> Vec<f64> {
        (0..len)
            .map(|i| (i as f64 + 1.0) * 1.000_000_1_f64.powi(i as i32 % 64))
            .collect()
    }

    #[divan::bench(consts = SIZES)]
    fn fsum<const LEN: usize>(bencher: Bencher) {
        let xs = values(LEN);
        bencher.bench(|| pymath::math::fsum(black_box(&xs).iter().copied()));
    }

    #[divan::bench(consts = SIZES)]
    fn dist<const LEN: usize>(bencher: Bencher) {
        let p = values(LEN);
        let q: Vec<f64> = values(LEN).iter().map(|v| v * 0.5).collect();
        bencher.bench(|| pymath::math::dist(black_box(&p), black_box(&q)));
    }

    #[divan::bench(consts = SIZES)]
    fn sumprod<const LEN: usize>(bencher: Bencher) {
        let p = values(LEN);
        let q: Vec<f64> = values(LEN).iter().map(|v| v * 0.5).collect();
        bencher.bench(|| pymath::math::sumprod(black_box(&p), black_box(&q)));
    }

    #[divan::bench(consts = SIZES)]
    fn sumprod_int<const LEN: usize>(bencher: Bencher) {
        let p: Vec<i64> = (0..LEN as i64).collect();
        let q: Vec<i64> = (0..LEN as i64).map(|i| i * 3 + 1).collect();
        bencher.bench(|| pymath::math::sumprod_int(black_box(&p), black_box(&q)));
    }

    #[divan::bench(consts = SIZES)]
    fn prod<const LEN: usize>(bencher: Bencher) {
        let xs: Vec<f64> = (0..LEN).map(|i| 1.0 + (i % 7) as f64 * 1e-3).collect();
        bencher.bench(|| pymath::math::prod(black_box(&xs).iter().copied(), None));
    }
}
