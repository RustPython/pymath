use crate::Error;
use pyo3::{Python, prelude::*};

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
