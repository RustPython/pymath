use crate::Error;
use pyo3::{Python, prelude::*};

pub(crate) fn unwrap<'a, T: 'a>(
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
