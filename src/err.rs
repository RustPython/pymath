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

/// Set errno to the given value.
#[inline]
pub(crate) fn set_errno(value: i32) {
    unsafe {
        #[cfg(target_os = "linux")]
        {
            *libc::__errno_location() = value;
        }
        #[cfg(target_os = "macos")]
        {
            *libc::__error() = value;
        }
        #[cfg(target_os = "windows")]
        {
            unsafe extern "C" {
                safe fn _errno() -> *mut i32;
            }
            *_errno() = value;
        }
        #[cfg(all(unix, not(any(target_os = "linux", target_os = "macos"))))]
        {
            // FreeBSD, NetBSD, OpenBSD, etc. use __error()
            *libc::__error() = value;
        }
    }
}

/// Get the current errno value.
#[inline]
pub(crate) fn get_errno() -> i32 {
    unsafe {
        #[cfg(target_os = "linux")]
        {
            *libc::__errno_location()
        }
        #[cfg(target_os = "macos")]
        {
            *libc::__error()
        }
        #[cfg(target_os = "windows")]
        {
            unsafe extern "C" {
                safe fn _errno() -> *mut i32;
            }
            *_errno()
        }
        #[cfg(all(unix, not(any(target_os = "linux", target_os = "macos"))))]
        {
            // FreeBSD, NetBSD, OpenBSD, etc. use __error()
            *libc::__error()
        }
    }
}

/// Check errno after libm call and convert to Result.
#[inline]
pub(crate) fn is_error(x: f64) -> Result<f64> {
    match get_errno() {
        0 => Ok(x),
        libc::EDOM => Err(Error::EDOM),
        libc::ERANGE => {
            // Underflow to zero is not an error.
            // Use 1.5 threshold to handle subnormal results that don't underflow to zero
            // (e.g., on Ubuntu/ia64) and to correctly detect underflows in expm1()
            // which may underflow toward -1.0 rather than 0.0. (bpo-46018)
            if x.abs() < 1.5 {
                Ok(x)
            } else {
                Err(Error::ERANGE)
            }
        }
        _ => Ok(x), // Unknown errno, just return the value
    }
}
