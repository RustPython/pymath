// The values are defined in libc
#[derive(Debug, PartialEq, Eq)]
pub enum Error {
    EDOM = 33,
    ERANGE = 34,
}

pub type Result<T> = std::result::Result<T, Error>;

impl TryFrom<libc::c_int> for Error {
    type Error = libc::c_int;

    fn try_from(value: libc::c_int) -> std::result::Result<Self, Self::Error> {
        match value {
            33 => Ok(Error::EDOM),
            34 => Ok(Error::ERANGE),
            _ => Err(value),
        }
    }
}
