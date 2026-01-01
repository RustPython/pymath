#[cfg(feature = "complex")]
pub mod cmath;
mod err;
mod gamma;
mod m;
#[cfg(test)]
mod test;

#[cfg(feature = "complex")]
pub use cmath::{c_acos, c_asin, c_asinh, c_sqrt};
pub use err::{Error, Result};
pub use gamma::{gamma, lgamma};

/// Fused multiply-add operation.
/// When `mul_add` feature is enabled, uses hardware FMA instruction.
/// Otherwise, uses separate multiply and add operations.
#[inline(always)]
pub(crate) fn mul_add(a: f64, b: f64, c: f64) -> f64 {
    if cfg!(feature = "mul_add") {
        a.mul_add(b, c)
    } else {
        a * b + c
    }
}

macro_rules! libm {
    // Reset errno and handle errno when return type contains Result
    (fn $name:ident($arg:ident: $ty:ty) -> Result<$ret:ty>) => {
        #[inline(always)]
        pub fn $name($arg: $ty) -> Result<$ret> {
            errno::set_errno(errno::Errno(0));
            let r = unsafe { m::$name($arg) };
            crate::is_error(r)
        }
    };
    // Skip errno checking when return type is not Result
    (fn $name:ident($arg:ident: $ty:ty) -> $ret:ty) => {
        #[inline(always)]
        pub fn $name($arg: $ty) -> $ret {
            unsafe { m::$name($arg) }
        }
    };
}

macro_rules! pyo3_proptest {
    ($fn_name:ident(Result<_>), $test_name:ident, $proptest_name:ident, $edgetest_name:ident) => {
        #[cfg(test)]
        fn $test_name(x: f64) {
            use pyo3::prelude::*;

            let rs_result = $fn_name(x);

            pyo3::Python::attach(|py| {
                let math = PyModule::import(py, "math").unwrap();
                let py_func = math
                    .getattr(stringify!($fn_name))
                    .unwrap();
                let r = py_func.call1((x,));
                let Some((py_result, rs_result)) = crate::test::unwrap(py, r, rs_result) else {
                    return;
                };
                let py_result_repr = py_result.to_bits();
                let rs_result_repr = rs_result.to_bits();
                assert_eq!(py_result_repr, rs_result_repr, "x = {x}, py_result = {py_result}, rs_result = {rs_result}");
            });
        }

        crate::pyo3_proptest!(@proptest, $test_name, $proptest_name);
        crate::pyo3_proptest!(@edgetest, $test_name, $edgetest_name);
    };
    ($fn_name:ident(_), $test_name:ident, $proptest_name:ident, $edgetest_name:ident) => {
        #[cfg(test)]
        fn $test_name(x: f64) {
            use pyo3::prelude::*;

            let rs_result = Ok($fn_name(x));

            pyo3::Python::attach(|py| {
                let math = PyModule::import(py, "math").unwrap();
                let py_func = math
                    .getattr(stringify!($fn_name))
                    .unwrap();
                let r = py_func.call1((x,));
                let Some((py_result, rs_result)) = crate::test::unwrap(py, r, rs_result) else {
                    return;
                };
                let py_result_repr = py_result.to_bits();
                let rs_result_repr = rs_result.to_bits();
                assert_eq!(py_result_repr, rs_result_repr, "x = {x}, py_result = {py_result}, rs_result = {rs_result}");
            });
        }
        crate::pyo3_proptest!(@proptest, $test_name, $proptest_name);
    };
    (@proptest, $test_name:ident, $proptest_name:ident) => {
        #[cfg(test)]
        proptest::proptest! {
            #[test]
            fn $proptest_name(x: f64) {
                $test_name(x);
            }
        }
    };
    (@edgetest, $test_name:ident, $edgetest_name:ident) => {
        #[test]
        fn $edgetest_name() {
            $test_name(f64::MIN);
            $test_name(-f64::MIN);
            $test_name(f64::NAN);
            $test_name(-f64::NAN);
            $test_name(f64::INFINITY);
            $test_name(-f64::NEG_INFINITY);
            $test_name(0.0);
            $test_name(-0.0);
        }
    };
}

libm!(fn erf(n: f64) -> f64);
pyo3_proptest!(erf(_), test_erf, proptest_erf, edgetest_erf);

libm!(fn erfc(n: f64) -> f64);
pyo3_proptest!(erfc(_), test_erfc, proptest_erfc, edgetest_erfc);

/// Call is_error when errno != 0, and where x is the result libm
/// returned.  is_error will usually set up an exception and return
/// true (1), but may return false (0) without setting up an exception.
// fn is_error(x: f64) -> crate::Result<f64> {
//     match errno::errno() {
//         errno::Errno(0) => Ok(x),
//         errno::Errno(libc::ERANGE) if x.abs() < 1.5 => Ok(0f64),
//         errno::Errno(errno) => Err(errno.try_into().unwrap()),
//     }
// }
use pyo3_proptest;
