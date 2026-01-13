// Raw FFI declarations for system libm.
// Only available on unix/windows where native libm exists.

#[cfg(any(unix, windows))]
#[cfg_attr(unix, link(name = "m"))]
#[allow(dead_code)]
unsafe extern "C" {
    // Trigonometric functions
    pub fn acos(n: f64) -> f64;
    pub fn asin(n: f64) -> f64;
    pub fn atan(n: f64) -> f64;
    pub fn atan2(a: f64, b: f64) -> f64;
    pub fn cos(n: f64) -> f64;
    pub fn sin(n: f64) -> f64;
    pub fn tan(n: f64) -> f64;

    // Hyperbolic functions
    pub fn acosh(n: f64) -> f64;
    pub fn asinh(n: f64) -> f64;
    pub fn atanh(n: f64) -> f64;
    pub fn cosh(n: f64) -> f64;
    pub fn sinh(n: f64) -> f64;
    pub fn tanh(n: f64) -> f64;

    // Exponential and logarithmic functions
    pub fn exp(n: f64) -> f64;
    pub fn exp2(n: f64) -> f64;
    pub fn expm1(n: f64) -> f64;
    pub fn expm1f(n: f32) -> f32;
    pub fn log(n: f64) -> f64;
    pub fn log10(n: f64) -> f64;
    pub fn log1p(n: f64) -> f64;
    pub fn log1pf(n: f32) -> f32;
    pub fn log2(n: f64) -> f64;

    // Power functions
    pub fn cbrt(n: f64) -> f64;
    pub fn cbrtf(n: f32) -> f32;
    #[cfg_attr(target_env = "msvc", link_name = "_hypot")]
    pub fn hypot(x: f64, y: f64) -> f64;
    #[cfg_attr(target_env = "msvc", link_name = "_hypotf")]
    pub fn hypotf(x: f32, y: f32) -> f32;
    pub fn pow(x: f64, y: f64) -> f64;
    pub fn sqrt(n: f64) -> f64;

    // Floating-point manipulation functions
    pub fn ceil(n: f64) -> f64;
    pub fn copysign(x: f64, y: f64) -> f64;
    pub fn fabs(n: f64) -> f64;
    pub fn fdimf(a: f32, b: f32) -> f32;
    pub fn floor(n: f64) -> f64;
    pub fn fmod(x: f64, y: f64) -> f64;
    pub fn frexp(n: f64, exp: *mut i32) -> f64;
    pub fn ldexp(x: f64, n: i32) -> f64;
    pub fn modf(n: f64, iptr: *mut f64) -> f64;
    pub fn nextafter(x: f64, y: f64) -> f64;
    pub fn remainder(x: f64, y: f64) -> f64;
    pub fn trunc(n: f64) -> f64;

    // Special functions
    pub fn erf(n: f64) -> f64;
    pub fn erfc(n: f64) -> f64;
    pub fn erff(n: f32) -> f32;
    pub fn erfcf(n: f32) -> f32;
    #[cfg(not(target_os = "aix"))]
    pub fn lgammaf_r(n: f32, s: &mut i32) -> f32;
    pub fn tgammaf(n: f32) -> f32;
}
