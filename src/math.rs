//! Real number mathematical functions matching Python's math module.

// Submodules
mod aggregate;
mod exponential;
mod gamma;
#[cfg(feature = "_bigint")]
pub mod integer;
mod misc;
mod trigonometric;

// Re-export from submodules
pub use aggregate::{dist, fsum, prod, sumprod};
pub use exponential::{cbrt, exp, exp2, expm1, log, log1p, log2, log10, pow, sqrt};
pub use gamma::{erf, erfc, gamma, lgamma};
pub use misc::{
    ceil, copysign, fabs, floor, fmod, frexp, isclose, isfinite, isinf, isnan, ldexp, modf,
    nextafter, remainder, trunc, ulp,
};
pub use trigonometric::{
    acos, acosh, asin, asinh, atan, atan2, atanh, cos, cosh, sin, sinh, tan, tanh,
};

/// Simple libm wrapper macro for functions that don't need errno handling.
macro_rules! libm_simple {
    // 1-arg: (f64) -> f64
    (@1 $($name:ident),* $(,)?) => {
        $(
            #[inline]
            pub fn $name(x: f64) -> f64 {
                crate::m::$name(x)
            }
        )*
    };
    // 2-arg: (f64, f64) -> f64
    (@2 $($name:ident),* $(,)?) => {
        $(
            #[inline]
            pub fn $name(x: f64, y: f64) -> f64 {
                crate::m::$name(x, y)
            }
        )*
    };
}

pub(crate) use libm_simple;

/// math_1: wrapper for 1-arg functions
/// - isnan(r) && !isnan(x) -> domain error
/// - isinf(r) && isfinite(x) -> overflow (can_overflow=true) or domain error (can_overflow=false)
/// - isfinite(r) && errno -> check errno (unnecessary on most platforms)
#[inline]
pub(crate) fn math_1(x: f64, func: fn(f64) -> f64, can_overflow: bool) -> crate::Result<f64> {
    crate::err::set_errno(0);
    let r = func(x);
    if r.is_nan() && !x.is_nan() {
        return Err(crate::Error::EDOM);
    }
    if r.is_infinite() && x.is_finite() {
        return Err(if can_overflow {
            crate::Error::ERANGE
        } else {
            crate::Error::EDOM
        });
    }
    // This branch unnecessary on most platforms
    #[cfg(not(any(windows, target_os = "macos")))]
    if r.is_finite() && crate::err::get_errno() != 0 {
        return crate::err::is_error(r);
    }
    Ok(r)
}

/// math_2: wrapper for 2-arg functions
/// - isnan(r) && !isnan(x) && !isnan(y) -> domain error
/// - isinf(r) && isfinite(x) && isfinite(y) -> range error
#[inline]
pub(crate) fn math_2(x: f64, y: f64, func: fn(f64, f64) -> f64) -> crate::Result<f64> {
    let r = func(x, y);
    if r.is_nan() && !x.is_nan() && !y.is_nan() {
        return Err(crate::Error::EDOM);
    }
    if r.is_infinite() && x.is_finite() && y.is_finite() {
        return Err(crate::Error::ERANGE);
    }
    Ok(r)
}

/// math_1a: wrapper for 1-arg functions that set errno properly.
/// Used when the libm function is known to set errno correctly
/// (EDOM for invalid, ERANGE for overflow).
#[inline]
pub(crate) fn math_1a(x: f64, func: fn(f64) -> f64) -> crate::Result<f64> {
    crate::err::set_errno(0);
    let r = func(x);
    crate::err::is_error(r)
}

/// Return the Euclidean distance, sqrt(x*x + y*y).
///
/// Uses high-precision vector_norm algorithm instead of libm hypot()
/// for consistent results across platforms and better handling of overflow/underflow.
#[inline]
pub fn hypot(x: f64, y: f64) -> f64 {
    let ax = x.abs();
    let ay = y.abs();
    let max = if ax > ay { ax } else { ay };
    let found_nan = x.is_nan() || y.is_nan();
    aggregate::vector_norm_2(ax, ay, max, found_nan)
}

// Mathematical constants

/// The mathematical constant π = 3.141592...
pub const PI: f64 = std::f64::consts::PI;

/// The mathematical constant e = 2.718281...
pub const E: f64 = std::f64::consts::E;

/// The mathematical constant τ = 6.283185...
pub const TAU: f64 = std::f64::consts::TAU;

/// Positive infinity.
pub const INF: f64 = f64::INFINITY;

/// A floating point "not a number" (NaN) value.
pub const NAN: f64 = f64::NAN;

// Angle conversion functions

/// Convert angle x from radians to degrees.
#[inline]
pub fn degrees(x: f64) -> f64 {
    x * (180.0 / PI)
}

/// Convert angle x from degrees to radians.
#[inline]
pub fn radians(x: f64) -> f64 {
    x * (PI / 180.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use pyo3::prelude::*;

    // Angle conversion tests
    fn test_degrees(x: f64) {
        let rs_result = Ok(degrees(x));

        pyo3::Python::attach(|py| {
            let math = PyModule::import(py, "math").unwrap();
            let py_func = math.getattr("degrees").unwrap();
            let r = py_func.call1((x,));
            let Some((py_result, rs_result)) = crate::test::unwrap(py, r, rs_result) else {
                return;
            };
            let py_result_repr = py_result.to_bits();
            let rs_result_repr = rs_result.to_bits();
            assert_eq!(
                py_result_repr, rs_result_repr,
                "x = {x}, py_result = {py_result}, rs_result = {rs_result}"
            );
        });
    }

    fn test_radians(x: f64) {
        let rs_result = Ok(radians(x));

        pyo3::Python::attach(|py| {
            let math = PyModule::import(py, "math").unwrap();
            let py_func = math.getattr("radians").unwrap();
            let r = py_func.call1((x,));
            let Some((py_result, rs_result)) = crate::test::unwrap(py, r, rs_result) else {
                return;
            };
            let py_result_repr = py_result.to_bits();
            let rs_result_repr = rs_result.to_bits();
            assert_eq!(
                py_result_repr, rs_result_repr,
                "x = {x}, py_result = {py_result}, rs_result = {rs_result}"
            );
        });
    }

    #[test]
    fn edgetest_degrees() {
        for &x in &crate::test::EDGE_VALUES {
            test_degrees(x);
        }
    }

    #[test]
    fn edgetest_radians() {
        for &x in &crate::test::EDGE_VALUES {
            test_radians(x);
        }
    }

    // Constants test
    #[test]
    fn test_constants() {
        assert!((PI - 3.141592653589793).abs() < 1e-15);
        assert!((E - 2.718281828459045).abs() < 1e-15);
        assert!((TAU - 6.283185307179586).abs() < 1e-15);
        assert!(INF.is_infinite() && INF > 0.0);
        assert!(NAN.is_nan());
    }

    // hypot tests
    fn test_hypot(x: f64, y: f64) {
        crate::test::test_math_2(x, y, "hypot", |x, y| Ok(hypot(x, y)));
    }

    #[test]
    fn edgetest_hypot() {
        for &x in &crate::test::EDGE_VALUES {
            for &y in &crate::test::EDGE_VALUES {
                test_hypot(x, y);
            }
        }
    }

    proptest::proptest! {
        #[test]
        fn proptest_degrees(x: f64) {
            test_degrees(x);
        }

        #[test]
        fn proptest_radians(x: f64) {
            test_radians(x);
        }

        #[test]
        fn proptest_hypot(x: f64, y: f64) {
            test_hypot(x, y);
        }
    }
}
