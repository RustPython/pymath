//! Complex exponential and logarithmic functions.

use super::{
    CM_LARGE_DOUBLE, CM_LOG_LARGE_DOUBLE, INF, M_LN2, N, P, P12, P14, P34, U, c, special_type,
    special_value,
};
use crate::{Error, Result, m};
use num_complex::Complex64;

// Local constants
const M_LN10: f64 = core::f64::consts::LN_10;

/// Scale factors for subnormal handling in sqrt.
const CM_SCALE_UP: i32 = 2 * (f64::MANTISSA_DIGITS as i32 / 2) + 1; // 54 for IEEE 754
const CM_SCALE_DOWN: i32 = -(CM_SCALE_UP + 1) / 2; // -27

// Special value tables

#[rustfmt::skip]
static EXP_SPECIAL_VALUES: [[Complex64; 7]; 7] = [
    [c(0.0, 0.0), c(U, U), c(0.0, -0.0), c(0.0, 0.0),  c(U, U), c(0.0, 0.0), c(0.0, 0.0)],
    [c(N, N),     c(U, U), c(U, U),      c(U, U),      c(U, U), c(N, N),     c(N, N)],
    [c(N, N),     c(U, U), c(1.0, -0.0), c(1.0, 0.0),  c(U, U), c(N, N),     c(N, N)],
    [c(N, N),     c(U, U), c(1.0, -0.0), c(1.0, 0.0),  c(U, U), c(N, N),     c(N, N)],
    [c(N, N),     c(U, U), c(U, U),      c(U, U),      c(U, U), c(N, N),     c(N, N)],
    [c(INF, N),   c(U, U), c(INF, -0.0), c(INF, 0.0),  c(U, U), c(INF, N),   c(INF, N)],
    [c(N, N),     c(N, N), c(N, -0.0),   c(N, 0.0),    c(N, N), c(N, N),     c(N, N)],
];

#[rustfmt::skip]
static LOG_SPECIAL_VALUES: [[Complex64; 7]; 7] = [
    [c(INF, -P34), c(INF, -P),    c(INF, -P),    c(INF, P),    c(INF, P),   c(INF, P34), c(INF, N)],
    [c(INF, -P12), c(U, U),       c(U, U),       c(U, U),      c(U, U),     c(INF, P12), c(N, N)],
    [c(INF, -P12), c(U, U),       c(-INF, -P),   c(-INF, P),   c(U, U),     c(INF, P12), c(N, N)],
    [c(INF, -P12), c(U, U),       c(-INF, -0.0), c(-INF, 0.0), c(U, U),     c(INF, P12), c(N, N)],
    [c(INF, -P12), c(U, U),       c(U, U),       c(U, U),      c(U, U),     c(INF, P12), c(N, N)],
    [c(INF, -P14), c(INF, -0.0),  c(INF, -0.0),  c(INF, 0.0),  c(INF, 0.0), c(INF, P14), c(INF, N)],
    [c(INF, N),    c(N, N),       c(N, N),       c(N, N),      c(N, N),     c(INF, N),   c(N, N)],
];

#[rustfmt::skip]
static SQRT_SPECIAL_VALUES: [[Complex64; 7]; 7] = [
    [c(INF, -INF), c(0.0, -INF), c(0.0, -INF), c(0.0, INF), c(0.0, INF), c(INF, INF), c(N, INF)],
    [c(INF, -INF), c(U, U),      c(U, U),      c(U, U),     c(U, U),     c(INF, INF), c(N, N)],
    [c(INF, -INF), c(U, U),      c(0.0, -0.0), c(0.0, 0.0), c(U, U),     c(INF, INF), c(N, N)],
    [c(INF, -INF), c(U, U),      c(0.0, -0.0), c(0.0, 0.0), c(U, U),     c(INF, INF), c(N, N)],
    [c(INF, -INF), c(U, U),      c(U, U),      c(U, U),     c(U, U),     c(INF, INF), c(N, N)],
    [c(INF, -INF), c(INF, -0.0), c(INF, -0.0), c(INF, 0.0), c(INF, 0.0), c(INF, INF), c(INF, N)],
    [c(INF, -INF), c(N, N),      c(N, N),      c(N, N),     c(N, N),     c(INF, INF), c(N, N)],
];

/// Complex square root.
///
/// Uses symmetries to reduce to the case when x = z.real and y = z.imag
/// are nonnegative, with careful handling of overflow and subnormals.
#[inline]
pub fn sqrt(z: Complex64) -> Result<Complex64> {
    special_value!(z, SQRT_SPECIAL_VALUES);

    if z.re == 0.0 && z.im == 0.0 {
        return Ok(Complex64::new(0.0, z.im));
    }

    let ax = m::fabs(z.re);
    let ay = m::fabs(z.im);

    let s = if ax < f64::MIN_POSITIVE && ay < f64::MIN_POSITIVE {
        // Handle subnormal case
        let ax_scaled = m::ldexp(ax, CM_SCALE_UP);
        m::ldexp(
            m::sqrt(ax_scaled + m::hypot(ax_scaled, m::ldexp(ay, CM_SCALE_UP))),
            CM_SCALE_DOWN,
        )
    } else {
        let ax8 = ax / 8.0;
        2.0 * m::sqrt(ax8 + m::hypot(ax8, ay / 8.0))
    };

    let d = ay / (2.0 * s);

    if z.re >= 0.0 {
        Ok(Complex64::new(s, m::copysign(d, z.im)))
    } else {
        Ok(Complex64::new(d, m::copysign(s, z.im)))
    }
}

/// Complex exponential.
#[inline]
pub fn exp(z: Complex64) -> Result<Complex64> {
    // Handle special values
    if !z.re.is_finite() || !z.im.is_finite() {
        let r = if z.re.is_infinite() && z.im.is_finite() && z.im != 0.0 {
            if z.re > 0.0 {
                Complex64::new(
                    m::copysign(INF, m::cos(z.im)),
                    m::copysign(INF, m::sin(z.im)),
                )
            } else {
                Complex64::new(
                    m::copysign(0.0, m::cos(z.im)),
                    m::copysign(0.0, m::sin(z.im)),
                )
            }
        } else {
            EXP_SPECIAL_VALUES[special_type(z.re) as usize][special_type(z.im) as usize]
        };
        // need to set errno = EDOM if y is +/- infinity and x is not a NaN and not -infinity
        if z.im.is_infinite() && (z.re.is_finite() || (z.re.is_infinite() && z.re > 0.0)) {
            return Err(Error::EDOM);
        }
        return Ok(r);
    }

    let (sin_im, cos_im) = m::sincos(z.im);
    let (r_re, r_im);
    if z.re > CM_LOG_LARGE_DOUBLE {
        let l = m::exp(z.re - 1.0);
        r_re = l * cos_im * core::f64::consts::E;
        r_im = l * sin_im * core::f64::consts::E;
    } else {
        let l = m::exp(z.re);
        r_re = l * cos_im;
        r_im = l * sin_im;
    }

    // detect overflow
    if r_re.is_infinite() || r_im.is_infinite() {
        return Err(Error::ERANGE);
    }
    Ok(Complex64::new(r_re, r_im))
}

/// Complex natural logarithm.
#[inline]
pub fn log(z: Complex64) -> Result<Complex64> {
    special_value!(z, LOG_SPECIAL_VALUES);

    let ax = m::fabs(z.re);
    let ay = m::fabs(z.im);

    let r_re = if ax > CM_LARGE_DOUBLE || ay > CM_LARGE_DOUBLE {
        m::log(m::hypot(ax / 2.0, ay / 2.0)) + M_LN2
    } else if ax < f64::MIN_POSITIVE && ay < f64::MIN_POSITIVE {
        if ax > 0.0 || ay > 0.0 {
            // catch cases where hypot(ax, ay) is subnormal
            m::log(m::hypot(
                m::ldexp(ax, f64::MANTISSA_DIGITS as i32),
                m::ldexp(ay, f64::MANTISSA_DIGITS as i32),
            )) - f64::MANTISSA_DIGITS as f64 * M_LN2
        } else {
            // log(+/-0. +/- 0i)
            return Err(Error::EDOM);
        }
    } else {
        let h = m::hypot(ax, ay);
        if (0.71..=1.73).contains(&h) {
            let am = if ax > ay { ax } else { ay }; // max(ax, ay)
            let an = if ax > ay { ay } else { ax }; // min(ax, ay)
            m::log1p((am - 1.0) * (am + 1.0) + an * an) / 2.0
        } else {
            m::log(h)
        }
    };

    let r_im = m::atan2(z.im, z.re);
    Ok(Complex64::new(r_re, r_im))
}

/// Complex base-10 logarithm.
#[inline]
pub fn log10(z: Complex64) -> Result<Complex64> {
    let r = log(z)?;
    Ok(Complex64::new(r.re / M_LN10, r.im / M_LN10))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_cmath_func<F>(func_name: &str, rs_func: F, re: f64, im: f64)
    where
        F: Fn(Complex64) -> Result<Complex64>,
    {
        crate::cmath::tests::test_cmath_func(func_name, rs_func, re, im);
    }

    fn test_sqrt(re: f64, im: f64) {
        test_cmath_func("sqrt", sqrt, re, im);
    }
    fn test_exp(re: f64, im: f64) {
        test_cmath_func("exp", exp, re, im);
    }
    fn test_log(re: f64, im: f64) {
        test_cmath_func("log", log, re, im);
    }
    fn test_log10(re: f64, im: f64) {
        test_cmath_func("log10", log10, re, im);
    }

    use crate::test::EDGE_VALUES;

    #[test]
    fn edgetest_sqrt() {
        for &re in &EDGE_VALUES {
            for &im in &EDGE_VALUES {
                test_sqrt(re, im);
            }
        }
    }

    #[test]
    fn edgetest_exp() {
        for &re in &EDGE_VALUES {
            for &im in &EDGE_VALUES {
                test_exp(re, im);
            }
        }
    }

    #[test]
    fn edgetest_log() {
        for &re in &EDGE_VALUES {
            for &im in &EDGE_VALUES {
                test_log(re, im);
            }
        }
    }

    #[test]
    fn edgetest_log10() {
        for &re in &EDGE_VALUES {
            for &im in &EDGE_VALUES {
                test_log10(re, im);
            }
        }
    }

    proptest::proptest! {
        #[test]
        fn proptest_sqrt(re: f64, im: f64) {
            test_sqrt(re, im);
        }

        #[test]
        fn proptest_exp(re: f64, im: f64) {
            test_exp(re, im);
        }

        #[test]
        fn proptest_log(re: f64, im: f64) {
            test_log(re, im);
        }

        #[test]
        fn proptest_log10(re: f64, im: f64) {
            test_log10(re, im);
        }
    }
}
