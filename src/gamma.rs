use crate::Error;
// use std::f64::consts::PI;

// Import C library functions properly
unsafe extern "C" {
    fn exp(x: f64) -> f64;
    fn log(x: f64) -> f64;
    fn sin(x: f64) -> f64;
    fn cos(x: f64) -> f64;
    fn pow(x: f64, y: f64) -> f64;
    fn floor(x: f64) -> f64;
    fn fabs(x: f64) -> f64;
    fn round(x: f64) -> f64;
    fn fmod(x: f64, y: f64) -> f64;
}

const PI: f64 = f64::from_bits(0x400921fb54442d18); // = 3.141592653589793238462643383279502884197;
const LOG_PI: f64 = f64::from_bits(0x3ff250d048e7a1bd); // = 1.144729885849400174143427351353058711647;

const LANCZOS_N: usize = 13;
const LANCZOS_G: f64 = f64::from_bits(0x40181945b9800000); // = 6.024680040776729583740234375;
const LANCZOS_G_MINUS_HALF: f64 = f64::from_bits(0x40161945b9800000); // = 5.524680040776729583740234375;
const LANCZOS_NUM_COEFFS: [f64; LANCZOS_N] = [
    f64::from_bits(0x4215ea5143c1a49e), // 23531376880.410759
    f64::from_bits(0x4223fc7075f54c57), // 42919803642.649101
    f64::from_bits(0x4220a132818ab61a), // 35711959237.355667
    f64::from_bits(0x4210b0b522e8261a), // 17921034426.037209
    f64::from_bits(0x41f67fc1b3a5a1e8), // 6039542586.3520279
    f64::from_bits(0x41d57418f5d3f33f), // 1439720407.3117216
    f64::from_bits(0x41adab0c7bb95f2a), // 248874557.86205417
    f64::from_bits(0x417df876f95dcc98), // 31426415.585400194
    f64::from_bits(0x4145f1e95080f44c), // 2876370.6289353725
    f64::from_bits(0x4106b6421f8787eb), // 186056.26539522348
    f64::from_bits(0x40bf87ac0858d804), // 8071.6720023658163
    f64::from_bits(0x406a5a607bbc3b52), // 210.82427775157936
    f64::from_bits(0x40040d931ff62705), // 2.5066282746310002
];
const LANCZOS_DEN_COEFFS: [f64; LANCZOS_N] = [
    f64::from_bits(0x0000000000000000), // 0.0
    f64::from_bits(0x418308a800000000), // 39916800.0
    f64::from_bits(0x419cbd6980000000), // 120543840.0
    f64::from_bits(0x41a1fda6b0000000), // 150917976.0
    f64::from_bits(0x4199187170000000), // 105258076.0
    f64::from_bits(0x4185eeb690000000), // 45995730.0
    f64::from_bits(0x41697171e0000000), // 13339535.0
    f64::from_bits(0x41441f7b00000000), // 2637558.0
    f64::from_bits(0x4115d0bc00000000), // 357423.0
    f64::from_bits(0x40dfe78000000000), // 32670.0
    f64::from_bits(0x409e140000000000), // 1925.0
    f64::from_bits(0x4050800000000000), // 66.0
    f64::from_bits(0x3ff0000000000000), // 1.0
];

const NGAMMA_INTEGRAL: usize = 23;
const GAMMA_INTEGRAL: [f64; NGAMMA_INTEGRAL] = [
    f64::from_bits(0x3ff0000000000000), // 1.0
    f64::from_bits(0x3ff0000000000000), // 1.0
    f64::from_bits(0x4000000000000000), // 2.0
    f64::from_bits(0x4018000000000000), // 6.0
    f64::from_bits(0x4038000000000000), // 24.0
    f64::from_bits(0x405e000000000000), // 120.0
    f64::from_bits(0x4086800000000000), // 720.0
    f64::from_bits(0x40b3b00000000000), // 5040.0
    f64::from_bits(0x40e3b00000000000), // 40320.0
    f64::from_bits(0x4116260000000000), // 362880.0
    f64::from_bits(0x414baf8000000000), // 3628800.0
    f64::from_bits(0x418308a800000000), // 39916800.0
    f64::from_bits(0x41bc8cfc00000000), // 479001600.0
    f64::from_bits(0x41f7328cc0000000), // 6227020800.0
    f64::from_bits(0x42344c3b28000000), // 87178291200.0
    f64::from_bits(0x4273077775800000), // 1307674368000.0
    f64::from_bits(0x42b3077775800000), // 20922789888000.0
    f64::from_bits(0x42f437eeecd80000), // 355687428096000.0
    f64::from_bits(0x4336beecca730000), // 6402373705728000.0
    f64::from_bits(0x437b02b930689000), // 1.21645100408832e+17
    f64::from_bits(0x43c0e1b3be415a00), // 2.43290200817664e+18
    f64::from_bits(0x4406283be9b5c620), // 5.109094217170944e+19
    f64::from_bits(0x444e77526159f06c), // 1.1240007277776077e+21
];

fn lanczos_sum(x: f64) -> f64 {
    let mut num = 0.0;
    let mut den = 0.0;
    // evaluate the rational function lanczos_sum(x).  For large
    // x, the obvious algorithm risks overflow, so we instead
    // rescale the denominator and numerator of the rational
    // function by x**(1-LANCZOS_N) and treat this as a
    // rational function in 1/x.  This also reduces the error for
    // larger x values.  The choice of cutoff point (5.0 below) is
    // somewhat arbitrary; in tests, smaller cutoff values than
    // this resulted in lower accuracy.
    if x < 5.0 {
        for i in (0..LANCZOS_N).rev() {
            num = num * x + LANCZOS_NUM_COEFFS[i];
            den = den * x + LANCZOS_DEN_COEFFS[i];
        }
    } else {
        for i in 0..LANCZOS_N {
            num = num / x + LANCZOS_NUM_COEFFS[i];
            den = den / x + LANCZOS_DEN_COEFFS[i];
        }
    }
    num / den
}

fn m_sinpi(x: f64) -> f64 {
    // this function should only ever be called for finite arguments
    debug_assert!(x.is_finite());
    let y = unsafe { fmod(fabs(x), 2.0) };
    let n = unsafe { round(2.0 * y) } as i32;
    let r = match n {
        0 => unsafe { sin(PI * y) },
        1 => unsafe { cos(PI * (y - 0.5)) },
        2 => {
            // N.B. -sin(pi*(y-1.0)) is *not* equivalent: it would give
            // -0.0 instead of 0.0 when y == 1.0.
            unsafe { sin(PI * (1.0 - y)) }
        }
        3 => unsafe { -cos(PI * (y - 1.5)) },
        4 => unsafe { sin(PI * (y - 2.0)) },
        _ => unreachable!(),
    };
    (1.0f64).copysign(x) * r
}

pub fn tgamma(x: f64) -> Result<f64, Error> {
    // special cases
    if !x.is_finite() {
        if x.is_nan() || x > 0.0 {
            // tgamma(nan) = nan, tgamma(inf) = inf
            return Ok(x);
        } else {
            // tgamma(-inf) = nan, invalid
            return Err((f64::NAN, Error::EDOM).1);
        }
    }
    if x == 0.0 {
        // tgamma(+-0.0) = +-inf, divide-by-zero
        let v = if x.is_sign_positive() {
            f64::INFINITY
        } else {
            f64::NEG_INFINITY
        };
        return Err((v, Error::EDOM).1);
    }
    // integer arguments
    if x == unsafe { floor(x) } {
        if x < 0.0 {
            // tgamma(n) = nan, invalid for
            return Err((f64::NAN, Error::EDOM).1);
        }
        if x < NGAMMA_INTEGRAL as f64 {
            return Ok(GAMMA_INTEGRAL[x as usize - 1]);
        }
    }
    let absx = unsafe { fabs(x) };
    // tiny arguments:  tgamma(x) ~ 1/x for x near 0
    if absx < 1e-20 {
        let r = 1.0 / x;
        if r.is_infinite() {
            return Err((f64::INFINITY, Error::ERANGE).1);
        } else {
            return Ok(r);
        }
    }
    // large arguments: assuming IEEE 754 doubles, tgamma(x) overflows for
    // x > 200, and underflows to +-0.0 for x < -200, not a negative
    // integer.
    if absx > 200.0 {
        if x < 0.0 {
            return Ok(0.0 / m_sinpi(x));
        } else {
            return Err((f64::INFINITY, Error::ERANGE).1);
        }
    }

    let y = absx + LANCZOS_G_MINUS_HALF;
    let z = if absx > LANCZOS_G_MINUS_HALF {
        // note: the correction can be foiled by an optimizing
        // compiler that (incorrectly) thinks that an expression like
        // a + b - a - b can be optimized to 0.0.  This shouldn't
        // happen in a standards-conforming compiler.
        let q = y - absx;
        q - LANCZOS_G_MINUS_HALF
    } else {
        let q = y - LANCZOS_G_MINUS_HALF;
        q - absx
    };
    let z = z * LANCZOS_G / y;
    let r = if x < 0.0 {
        // Using C's math functions through libc to match CPython
        let term1 = -PI / m_sinpi(absx);
        let term2 = term1 / absx;
        let exp_y = unsafe { exp(y) };
        let term3 = term2 * exp_y;
        let lanczos = lanczos_sum(absx);
        let mut r = term3 / lanczos;
        r -= z * r;

        if absx < 140.0 {
            unsafe { r / pow(y, absx - 0.5) }
        } else {
            let sqrtpow = unsafe { pow(y, absx / 2.0 - 0.25) };
            r /= sqrtpow;
            r /= sqrtpow;
            r
        }
    } else {
        let lanczos = lanczos_sum(absx);
        let exp_y = unsafe { exp(y) };
        let mut r = lanczos / exp_y;
        r += z * r;

        if absx < 140.0 {
            unsafe { r * pow(y, absx - 0.5) }
        } else {
            let sqrtpow = unsafe { pow(y, absx / 2.0 - 0.25) };
            r *= sqrtpow;
            r *= sqrtpow;
            r
        }
    };

    if r.is_infinite() {
        return Err((f64::INFINITY, Error::ERANGE).1);
    } else {
        return Ok(r);
    }
}

// natural log of the absolute value of the Gamma function.
// For large arguments, Lanczos' formula works extremely well here.
pub fn lgamma(x: f64) -> Result<f64, Error> {
    // special cases
    if !x.is_finite() {
        if x.is_nan() {
            return Ok(x); // lgamma(nan) = nan
        } else {
            return Ok(f64::INFINITY); // lgamma(+-inf) = +inf
        }
    }

    // integer arguments
    if x == unsafe { floor(x) } && x <= 2.0 {
        if x <= 0.0 {
            // lgamma(n) = inf, divide-by-zero for integers n <= 0
            return Err(Error::EDOM);
        } else {
            // lgamma(1) = lgamma(2) = 0.0
            return Ok(0.0);
        }
    }

    let absx = unsafe { fabs(x) };
    // tiny arguments: lgamma(x) ~ -log(fabs(x)) for small x
    if absx < 1e-20 {
        return Ok(-unsafe { log(absx) });
    }

    // Using C's math functions through libc to match CPython
    let lanczos_sum_val = lanczos_sum(absx);
    let log_lanczos = unsafe { log(lanczos_sum_val) };

    // Subtract lanczos_g as a separate step
    let mut r = log_lanczos - LANCZOS_G;

    // Calculate (absx - 0.5) term
    let factor = absx - 0.5;

    // Calculate log term
    let log_term = unsafe { log(absx + LANCZOS_G - 0.5) };

    // Calculate the multiplication and subtraction
    let step2 = factor * (log_term - 1.0);

    // Combine the results
    r += step2;

    if x < 0.0 {
        // Calculate each component separately as in CPython
        let sinpi_val = m_sinpi(absx);
        let abs_sinpi = unsafe { fabs(sinpi_val) };
        let log_abs_sinpi = unsafe { log(abs_sinpi) };
        let log_absx = unsafe { log(absx) };

        // Combine in exactly the same order as CPython
        r = LOG_PI - log_abs_sinpi - log_absx - r;
    }

    if r.is_infinite() {
        return Err(Error::ERANGE);
    }
    Ok(r)
}

#[cfg(test)]
mod tests {
    use super::*;
    use pyo3::Python;
    use pyo3::prelude::*;

    use proptest::prelude::*;

    fn unwrap<'a, T: 'a>(
        py: Python,
        py_v: PyResult<Bound<'a, PyAny>>,
        v: Result<T, crate::Error>,
    ) -> Option<(T, T)>
    where
        T: PartialEq + std::fmt::Debug + FromPyObject<'a>,
    {
        match py_v {
            Ok(py_v) => {
                let py_v: T = py_v.extract().unwrap();
                Some((py_v, v.unwrap()))
            }
            Err(e) => {
                if e.is_instance_of::<pyo3::exceptions::PyValueError>(py) {
                    assert_eq!(v.err(), Some(Error::EDOM));
                } else if e.is_instance_of::<pyo3::exceptions::PyOverflowError>(py) {
                    assert_eq!(v.err(), Some(Error::ERANGE));
                } else {
                    panic!();
                }
                None
            }
        }
    }

    #[test]
    fn test_literal() {
        // Verify single constants
        assert_eq!(PI, 3.141592653589793238462643383279502884197);
        assert_eq!(LOG_PI, 1.144729885849400174143427351353058711647);
        assert_eq!(LANCZOS_G, 6.024680040776729583740234375);
        assert_eq!(LANCZOS_G_MINUS_HALF, 5.524680040776729583740234375);

        // Verify LANCZOS_NUM_COEFFS
        assert_eq!(LANCZOS_NUM_COEFFS[0], 23531376880.410759);
        assert_eq!(LANCZOS_NUM_COEFFS[1], 42919803642.649101);
        assert_eq!(LANCZOS_NUM_COEFFS[2], 35711959237.355667);
        assert_eq!(LANCZOS_NUM_COEFFS[3], 17921034426.037209);
        assert_eq!(LANCZOS_NUM_COEFFS[4], 6039542586.3520279);
        assert_eq!(LANCZOS_NUM_COEFFS[5], 1439720407.3117216);
        assert_eq!(LANCZOS_NUM_COEFFS[6], 248874557.86205417);
        assert_eq!(LANCZOS_NUM_COEFFS[7], 31426415.585400194);
        assert_eq!(LANCZOS_NUM_COEFFS[8], 2876370.6289353725);
        assert_eq!(LANCZOS_NUM_COEFFS[9], 186056.26539522348);
        assert_eq!(LANCZOS_NUM_COEFFS[10], 8071.6720023658163);
        assert_eq!(LANCZOS_NUM_COEFFS[11], 210.82427775157936);
        assert_eq!(LANCZOS_NUM_COEFFS[12], 2.5066282746310002);

        // Verify LANCZOS_DEN_COEFFS
        assert_eq!(LANCZOS_DEN_COEFFS[0], 0.0);
        assert_eq!(LANCZOS_DEN_COEFFS[1], 39916800.0);
        assert_eq!(LANCZOS_DEN_COEFFS[2], 120543840.0);
        assert_eq!(LANCZOS_DEN_COEFFS[3], 150917976.0);
        assert_eq!(LANCZOS_DEN_COEFFS[4], 105258076.0);
        assert_eq!(LANCZOS_DEN_COEFFS[5], 45995730.0);
        assert_eq!(LANCZOS_DEN_COEFFS[6], 13339535.0);
        assert_eq!(LANCZOS_DEN_COEFFS[7], 2637558.0);
        assert_eq!(LANCZOS_DEN_COEFFS[8], 357423.0);
        assert_eq!(LANCZOS_DEN_COEFFS[9], 32670.0);
        assert_eq!(LANCZOS_DEN_COEFFS[10], 1925.0);
        assert_eq!(LANCZOS_DEN_COEFFS[11], 66.0);
        assert_eq!(LANCZOS_DEN_COEFFS[12], 1.0);

        // Verify GAMMA_INTEGRAL
        assert_eq!(GAMMA_INTEGRAL[0], 1.0);
        assert_eq!(GAMMA_INTEGRAL[1], 1.0);
        assert_eq!(GAMMA_INTEGRAL[2], 2.0);
        assert_eq!(GAMMA_INTEGRAL[3], 6.0);
        assert_eq!(GAMMA_INTEGRAL[4], 24.0);
        assert_eq!(GAMMA_INTEGRAL[5], 120.0);
        assert_eq!(GAMMA_INTEGRAL[6], 720.0);
        assert_eq!(GAMMA_INTEGRAL[7], 5040.0);
        assert_eq!(GAMMA_INTEGRAL[8], 40320.0);
        assert_eq!(GAMMA_INTEGRAL[9], 362880.0);
        assert_eq!(GAMMA_INTEGRAL[10], 3628800.0);
        assert_eq!(GAMMA_INTEGRAL[11], 39916800.0);
        assert_eq!(GAMMA_INTEGRAL[12], 479001600.0);
        assert_eq!(GAMMA_INTEGRAL[13], 6227020800.0);
        assert_eq!(GAMMA_INTEGRAL[14], 87178291200.0);
        assert_eq!(GAMMA_INTEGRAL[15], 1307674368000.0);
        assert_eq!(GAMMA_INTEGRAL[16], 20922789888000.0);
        assert_eq!(GAMMA_INTEGRAL[17], 355687428096000.0);
        assert_eq!(GAMMA_INTEGRAL[18], 6402373705728000.0);
        assert_eq!(GAMMA_INTEGRAL[19], 1.21645100408832e+17);
        assert_eq!(GAMMA_INTEGRAL[20], 2.43290200817664e+18);
        assert_eq!(GAMMA_INTEGRAL[21], 5.109094217170944e+19);
        assert_eq!(GAMMA_INTEGRAL[22], 1.1240007277776077e+21);
    }

    #[test]
    fn test_specific_lgamma_value() {
        let x = 0.003585187864492183;
        let rs_lgamma = lgamma(x).unwrap();

        pyo3::prepare_freethreaded_python();
        Python::with_gil(|py| {
            let math = PyModule::import(py, "math").unwrap();
            let py_lgamma = math
                .getattr("lgamma")
                .unwrap()
                .call1((x,))
                .unwrap()
                .extract::<f64>()
                .unwrap();

            println!("x = {}", x);
            println!("Python lgamma = {} ({:x})", py_lgamma, unsafe {
                std::mem::transmute::<f64, u64>(py_lgamma)
            });
            println!("Rust lgamma = {} ({:x})", rs_lgamma, unsafe {
                std::mem::transmute::<f64, u64>(rs_lgamma)
            });

            // Print intermediate values
            let absx = x.abs();
            let sinpi_val = m_sinpi(absx);

            println!("absx = {}", absx);
            println!("m_sinpi = {}", sinpi_val);

            // Compare with Python's sin(pi * x)
            let py_code = PyModule::from_code(
                py,
                c"import math\ndef sinpi(x): return math.sin(math.pi * x)\n",
                c"",
                c"",
            )
            .unwrap();
            let py_sinpi = py_code
                .getattr("sinpi")
                .unwrap()
                .call1((absx,))
                .unwrap()
                .extract::<f64>()
                .unwrap();
            println!("Python sinpi = {}", py_sinpi);

            let py_lgamma_repr = unsafe { std::mem::transmute::<f64, u64>(py_lgamma) };
            let rs_lgamma_repr = unsafe { std::mem::transmute::<f64, u64>(rs_lgamma) };
            println!("Bit difference: {}", py_lgamma_repr ^ rs_lgamma_repr);

            assert_eq!(py_lgamma_repr, rs_lgamma_repr);
        });
    }

    proptest! {
        #[test]
        fn test_tgamma(x: f64) {
            let rs_gamma = tgamma(x);

            pyo3::prepare_freethreaded_python();
            Python::with_gil(|py| {
                let math = PyModule::import(py, "math").unwrap();
                let py_gamma_func = math
                    .getattr("gamma")
                    .unwrap();
                let r = py_gamma_func.call1((x,));
                let Some((py_gamma, rs_gamma)) = unwrap(py, r, rs_gamma) else {
                    return;
                };
                let py_gamma_repr = py_gamma.to_bits();
                let rs_gamma_repr = rs_gamma.to_bits();
                assert_eq!(py_gamma_repr, rs_gamma_repr, "x = {x}, py_gamma = {py_gamma}, rs_gamma = {rs_gamma}");
            });
        }

        #[test]
        fn test_lgamma(x: f64) {
            let rs_lgamma = lgamma(x);

            pyo3::prepare_freethreaded_python();
            Python::with_gil(|py| {
                let math = PyModule::import(py, "math").unwrap();
                let py_lgamma_func = math
                    .getattr("lgamma")
                    .unwrap();
                let r = py_lgamma_func.call1((x,));
                let Some((py_lgamma, rs_lgamma)) = unwrap(py, r, rs_lgamma) else {
                    return;
                };
                let py_lgamma_repr = py_lgamma.to_bits();
                let rs_lgamma_repr = rs_lgamma.to_bits();
                assert_eq!(py_lgamma_repr, rs_lgamma_repr, "x = {x}, py_lgamma = {py_lgamma}, rs_gamma = {rs_lgamma}");
            });
        }
    }
}
