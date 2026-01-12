//! Safe wrappers for libm functions.
//! Uses native libm via FFI on unix/windows, falls back to libm crate on other targets (WASM).

#[cfg(any(unix, windows))]
use crate::m_sys;

// Trigonometric functions

#[inline(always)]
pub fn acos(n: f64) -> f64 {
    #[cfg(any(unix, windows))]
    unsafe {
        m_sys::acos(n)
    }
    #[cfg(not(any(unix, windows)))]
    {
        libm::acos(n)
    }
}

#[inline(always)]
pub fn asin(n: f64) -> f64 {
    #[cfg(any(unix, windows))]
    unsafe {
        m_sys::asin(n)
    }
    #[cfg(not(any(unix, windows)))]
    {
        libm::asin(n)
    }
}

#[inline(always)]
pub fn atan(n: f64) -> f64 {
    #[cfg(any(unix, windows))]
    unsafe {
        m_sys::atan(n)
    }
    #[cfg(not(any(unix, windows)))]
    {
        libm::atan(n)
    }
}

#[inline(always)]
pub fn atan2(y: f64, x: f64) -> f64 {
    #[cfg(any(unix, windows))]
    unsafe {
        m_sys::atan2(y, x)
    }
    #[cfg(not(any(unix, windows)))]
    {
        libm::atan2(y, x)
    }
}

#[inline(always)]
pub fn cos(n: f64) -> f64 {
    #[cfg(any(unix, windows))]
    unsafe {
        m_sys::cos(n)
    }
    #[cfg(not(any(unix, windows)))]
    {
        libm::cos(n)
    }
}

#[inline(always)]
pub fn sin(n: f64) -> f64 {
    #[cfg(any(unix, windows))]
    unsafe {
        m_sys::sin(n)
    }
    #[cfg(not(any(unix, windows)))]
    {
        libm::sin(n)
    }
}

#[inline(always)]
pub fn tan(n: f64) -> f64 {
    #[cfg(any(unix, windows))]
    unsafe {
        m_sys::tan(n)
    }
    #[cfg(not(any(unix, windows)))]
    {
        libm::tan(n)
    }
}

// Hyperbolic functions

#[inline(always)]
pub fn acosh(n: f64) -> f64 {
    #[cfg(any(unix, windows))]
    unsafe {
        m_sys::acosh(n)
    }
    #[cfg(not(any(unix, windows)))]
    {
        libm::acosh(n)
    }
}

#[inline(always)]
pub fn asinh(n: f64) -> f64 {
    #[cfg(any(unix, windows))]
    unsafe {
        m_sys::asinh(n)
    }
    #[cfg(not(any(unix, windows)))]
    {
        libm::asinh(n)
    }
}

#[inline(always)]
pub fn atanh(n: f64) -> f64 {
    #[cfg(any(unix, windows))]
    unsafe {
        m_sys::atanh(n)
    }
    #[cfg(not(any(unix, windows)))]
    {
        libm::atanh(n)
    }
}

#[inline(always)]
pub fn cosh(n: f64) -> f64 {
    #[cfg(any(unix, windows))]
    unsafe {
        m_sys::cosh(n)
    }
    #[cfg(not(any(unix, windows)))]
    {
        libm::cosh(n)
    }
}

#[inline(always)]
pub fn sinh(n: f64) -> f64 {
    #[cfg(any(unix, windows))]
    unsafe {
        m_sys::sinh(n)
    }
    #[cfg(not(any(unix, windows)))]
    {
        libm::sinh(n)
    }
}

#[inline(always)]
pub fn tanh(n: f64) -> f64 {
    #[cfg(any(unix, windows))]
    unsafe {
        m_sys::tanh(n)
    }
    #[cfg(not(any(unix, windows)))]
    {
        libm::tanh(n)
    }
}

// Exponential and logarithmic functions

#[inline(always)]
pub fn exp(n: f64) -> f64 {
    #[cfg(any(unix, windows))]
    unsafe {
        m_sys::exp(n)
    }
    #[cfg(not(any(unix, windows)))]
    {
        libm::exp(n)
    }
}

#[inline(always)]
pub fn exp2(n: f64) -> f64 {
    #[cfg(any(unix, windows))]
    unsafe {
        m_sys::exp2(n)
    }
    #[cfg(not(any(unix, windows)))]
    {
        libm::exp2(n)
    }
}

#[inline(always)]
pub fn expm1(n: f64) -> f64 {
    #[cfg(any(unix, windows))]
    unsafe {
        m_sys::expm1(n)
    }
    #[cfg(not(any(unix, windows)))]
    {
        libm::expm1(n)
    }
}

#[inline(always)]
pub fn log(n: f64) -> f64 {
    #[cfg(any(unix, windows))]
    unsafe {
        m_sys::log(n)
    }
    #[cfg(not(any(unix, windows)))]
    {
        libm::log(n)
    }
}

#[inline(always)]
pub fn log10(n: f64) -> f64 {
    #[cfg(any(unix, windows))]
    unsafe {
        m_sys::log10(n)
    }
    #[cfg(not(any(unix, windows)))]
    {
        libm::log10(n)
    }
}

#[inline(always)]
pub fn log1p(n: f64) -> f64 {
    #[cfg(any(unix, windows))]
    unsafe {
        m_sys::log1p(n)
    }
    #[cfg(not(any(unix, windows)))]
    {
        libm::log1p(n)
    }
}

#[inline(always)]
pub fn log2(n: f64) -> f64 {
    #[cfg(any(unix, windows))]
    unsafe {
        m_sys::log2(n)
    }
    #[cfg(not(any(unix, windows)))]
    {
        libm::log2(n)
    }
}

// Power functions

#[inline(always)]
pub fn cbrt(n: f64) -> f64 {
    #[cfg(any(unix, windows))]
    unsafe {
        m_sys::cbrt(n)
    }
    #[cfg(not(any(unix, windows)))]
    {
        libm::cbrt(n)
    }
}

#[inline(always)]
pub fn hypot(x: f64, y: f64) -> f64 {
    #[cfg(any(unix, windows))]
    unsafe {
        m_sys::hypot(x, y)
    }
    #[cfg(not(any(unix, windows)))]
    {
        libm::hypot(x, y)
    }
}

#[inline(always)]
pub fn pow(x: f64, y: f64) -> f64 {
    #[cfg(any(unix, windows))]
    unsafe {
        m_sys::pow(x, y)
    }
    #[cfg(not(any(unix, windows)))]
    {
        libm::pow(x, y)
    }
}

#[inline(always)]
pub fn sqrt(n: f64) -> f64 {
    #[cfg(any(unix, windows))]
    unsafe {
        m_sys::sqrt(n)
    }
    #[cfg(not(any(unix, windows)))]
    {
        libm::sqrt(n)
    }
}

// Floating-point manipulation functions

#[inline(always)]
pub fn ceil(n: f64) -> f64 {
    #[cfg(any(unix, windows))]
    unsafe {
        m_sys::ceil(n)
    }
    #[cfg(not(any(unix, windows)))]
    {
        libm::ceil(n)
    }
}

#[inline(always)]
pub fn copysign(x: f64, y: f64) -> f64 {
    #[cfg(any(unix, windows))]
    unsafe {
        m_sys::copysign(x, y)
    }
    #[cfg(not(any(unix, windows)))]
    {
        libm::copysign(x, y)
    }
}

#[inline(always)]
pub fn fabs(n: f64) -> f64 {
    #[cfg(any(unix, windows))]
    unsafe {
        m_sys::fabs(n)
    }
    #[cfg(not(any(unix, windows)))]
    {
        libm::fabs(n)
    }
}

#[inline(always)]
pub fn floor(n: f64) -> f64 {
    #[cfg(any(unix, windows))]
    unsafe {
        m_sys::floor(n)
    }
    #[cfg(not(any(unix, windows)))]
    {
        libm::floor(n)
    }
}

#[inline(always)]
pub fn fmod(x: f64, y: f64) -> f64 {
    #[cfg(any(unix, windows))]
    unsafe {
        m_sys::fmod(x, y)
    }
    #[cfg(not(any(unix, windows)))]
    {
        libm::fmod(x, y)
    }
}

#[inline(always)]
pub fn frexp(n: f64, exp: &mut i32) -> f64 {
    #[cfg(any(unix, windows))]
    unsafe {
        m_sys::frexp(n, exp)
    }
    #[cfg(not(any(unix, windows)))]
    {
        let (mantissa, exponent) = libm::frexp(n);
        *exp = exponent;
        mantissa
    }
}

#[inline(always)]
pub fn ldexp(x: f64, n: i32) -> f64 {
    #[cfg(any(unix, windows))]
    unsafe {
        m_sys::ldexp(x, n)
    }
    #[cfg(not(any(unix, windows)))]
    {
        libm::ldexp(x, n)
    }
}

#[inline(always)]
pub fn modf(n: f64, iptr: &mut f64) -> f64 {
    #[cfg(any(unix, windows))]
    unsafe {
        m_sys::modf(n, iptr)
    }
    #[cfg(not(any(unix, windows)))]
    {
        let (frac, int) = libm::modf(n);
        *iptr = int;
        frac
    }
}

#[inline(always)]
pub fn nextafter(x: f64, y: f64) -> f64 {
    #[cfg(any(unix, windows))]
    unsafe {
        m_sys::nextafter(x, y)
    }
    #[cfg(not(any(unix, windows)))]
    {
        libm::nextafter(x, y)
    }
}

#[inline(always)]
pub fn remainder(x: f64, y: f64) -> f64 {
    #[cfg(any(unix, windows))]
    unsafe {
        m_sys::remainder(x, y)
    }
    #[cfg(not(any(unix, windows)))]
    {
        libm::remainder(x, y)
    }
}

#[inline(always)]
pub fn trunc(n: f64) -> f64 {
    #[cfg(any(unix, windows))]
    unsafe {
        m_sys::trunc(n)
    }
    #[cfg(not(any(unix, windows)))]
    {
        libm::trunc(n)
    }
}

// Special functions

#[inline(always)]
pub fn erf(n: f64) -> f64 {
    #[cfg(any(unix, windows))]
    unsafe {
        m_sys::erf(n)
    }
    #[cfg(not(any(unix, windows)))]
    {
        libm::erf(n)
    }
}

#[inline(always)]
pub fn erfc(n: f64) -> f64 {
    #[cfg(any(unix, windows))]
    unsafe {
        m_sys::erfc(n)
    }
    #[cfg(not(any(unix, windows)))]
    {
        libm::erfc(n)
    }
}

// Platform-specific sincos

/// Result type for sincos function (matches Apple's __double2)
#[cfg(all(feature = "complex", target_os = "macos"))]
#[repr(C)]
struct SinCosResult {
    sin: f64,
    cos: f64,
}

#[cfg(all(feature = "complex", target_os = "macos"))]
unsafe extern "C" {
    #[link_name = "__sincos_stret"]
    fn sincos_stret(x: f64) -> SinCosResult;
}

/// Compute sin and cos together using Apple's optimized sincos.
/// This matches Python's cmath behavior on macOS.
#[cfg(all(feature = "complex", target_os = "macos"))]
#[inline(always)]
pub fn sincos(x: f64) -> (f64, f64) {
    let sc = unsafe { sincos_stret(x) };
    (sc.sin, sc.cos)
}

/// Fallback for non-macOS: call sin and cos separately
#[cfg(all(feature = "complex", not(target_os = "macos")))]
#[inline(always)]
pub fn sincos(x: f64) -> (f64, f64) {
    (sin(x), cos(x))
}
