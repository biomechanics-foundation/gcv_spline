# gcv_spline - Generalized Cross-Validated Splines for Interpolation and Derivation in Pure Rust

![Crates.io](https://img.shields.io/crates/v/gcv_spline.svg)

This Rust crate implements the GCV spline, a versatile and easy-to-use spline structure for interpolating data at unknown points and taking accurate derivatives of smooth data.

GCV splines were first developed by Herman J. Woltring. The modules inside the private woltring module are based on his [FORTRAN package](https://isbweb.org/software/sigproc/gcvspl/gcvspl.f) and a [C translation](https://isbweb.org/software/sigproc/gcvspl/Twisk/gcvspl.c) by D. Twisk. Comments from these versions are included in this implementation.
