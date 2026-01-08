use crate::Error;
use pyo3::{Python, prelude::*};

/// Edge values for testing floating-point functions.
/// Includes: zeros, infinities, various NaNs, subnormals, and values at different scales.
pub(crate) const EDGE_VALUES: [f64; 30] = [
    // Zeros
    0.0,
    -0.0,
    // Infinities
    f64::INFINITY,
    f64::NEG_INFINITY,
    // Standard NaNs
    f64::NAN,
    -f64::NAN,
    // Additional NaN with different payload (quiet NaN with payload 1)
    f64::from_bits(0x7FF8_0000_0000_0001_u64),
    // Subnormal (denormalized) values
    f64::MIN_POSITIVE * 0.5, // smallest subnormal
    -f64::MIN_POSITIVE * 0.5,
    // Boundary values
    f64::MIN_POSITIVE, // smallest positive normal
    f64::MAX,          // largest finite
    f64::MIN,          // most negative finite (not smallest!)
    // Near-infinity large values
    f64::MAX * 0.5,
    -f64::MAX * 0.5,
    1e308,
    -1e308,
    // Small scale
    1e-10,
    -1e-10,
    1e-300,
    // Normal scale
    1.0,
    -1.0,
    0.5,
    -0.5,
    2.0,
    // Trigonometric special values (where sin/cos/tan have exact or near-zero results)
    std::f64::consts::PI, // sin(PI) ≈ 0
    -std::f64::consts::PI,
    std::f64::consts::FRAC_PI_2, // cos(PI/2) ≈ 0
    -std::f64::consts::FRAC_PI_2,
    std::f64::consts::FRAC_PI_4, // tan(PI/4) = 1
    std::f64::consts::TAU,       // sin(2*PI) ≈ 0, cos(2*PI) = 1
];

/// Edge integer values for testing functions like ldexp
pub(crate) const EDGE_INTS: [i32; 9] = [0, 1, -1, 100, -100, 1024, -1024, i32::MAX, i32::MIN];

pub(crate) fn unwrap<'py>(
    py: Python<'py>,
    py_v: PyResult<Bound<'py, PyAny>>,
    v: Result<f64, crate::Error>,
) -> Option<(f64, f64)> {
    match py_v {
        Ok(py_v) => {
            let py_v: f64 = py_v.extract().ok().expect("failed to extract");
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

/// Test a 1-argument function that returns Result<f64>
pub(crate) fn test_math_1(x: f64, func_name: &str, rs_func: impl Fn(f64) -> crate::Result<f64>) {
    let rs_result = rs_func(x);

    pyo3::Python::attach(|py| {
        let math = pyo3::types::PyModule::import(py, "math").unwrap();
        let py_func = math.getattr(func_name).unwrap();
        let r = py_func.call1((x,));
        let Some((py_result, rs_result)) = unwrap(py, r, rs_result) else {
            return;
        };
        if py_result.is_nan() && rs_result.is_nan() {
            return;
        }
        assert_eq!(
            py_result.to_bits(),
            rs_result.to_bits(),
            "{func_name}({x}): py={py_result} vs rs={rs_result}"
        );
    });
}

/// Test a 2-argument function that returns Result<f64>
pub(crate) fn test_math_2(
    x: f64,
    y: f64,
    func_name: &str,
    rs_func: impl Fn(f64, f64) -> crate::Result<f64>,
) {
    let rs_result = rs_func(x, y);

    pyo3::Python::attach(|py| {
        let math = pyo3::types::PyModule::import(py, "math").unwrap();
        let py_func = math.getattr(func_name).unwrap();
        let r = py_func.call1((x, y));

        match r {
            Ok(py_val) => {
                let py_f: f64 = py_val.extract().unwrap();
                let rs_val = rs_result.unwrap_or_else(|e| {
                    panic!("{func_name}({x}, {y}): py={py_f} but rs returned error {e:?}")
                });
                if py_f.is_nan() && rs_val.is_nan() {
                    return;
                }
                assert_eq!(
                    py_f.to_bits(),
                    rs_val.to_bits(),
                    "{func_name}({x}, {y}): py={py_f} vs rs={rs_val}"
                );
            }
            Err(e) => {
                // Check error type matches
                let rs_err = rs_result.as_ref().err();
                if e.is_instance_of::<pyo3::exceptions::PyValueError>(py) {
                    assert_eq!(
                        rs_err,
                        Some(&Error::EDOM),
                        "{func_name}({x}, {y}): py raised ValueError but rs={rs_err:?}"
                    );
                } else if e.is_instance_of::<pyo3::exceptions::PyOverflowError>(py) {
                    assert_eq!(
                        rs_err,
                        Some(&Error::ERANGE),
                        "{func_name}({x}, {y}): py raised OverflowError but rs={rs_err:?}"
                    );
                } else {
                    panic!("{func_name}({x}, {y}): py raised unexpected error {e}");
                }
            }
        }
    });
}
