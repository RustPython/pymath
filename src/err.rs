use std::ffi::c_int;

// The values are defined in libc
#[derive(Debug, PartialEq, Eq)]
pub enum Error {
    EDOM = 33,
    ERANGE = 34,
}

pub type Result<T> = std::result::Result<T, Error>;

impl TryFrom<c_int> for Error {
    type Error = c_int;

    fn try_from(value: c_int) -> std::result::Result<Self, Self::Error> {
        match value {
            33 => Ok(Error::EDOM),
            34 => Ok(Error::ERANGE),
            _ => Err(value),
        }
    }
}
