#![feature(generic_const_exprs)]
#![feature(associated_const_equality)]
#![feature(associated_type_bounds)]

pub mod dft;
pub mod fft;
pub mod fp17;
pub mod coeffs;
pub mod evals;

use std::usize;
use core::ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Sub, SubAssign};
use core::fmt::Debug;

use p3_field::{TwoAdicField,AbstractField};


pub trait AbstractPolynomial<F: AbstractField>:
    Default
    + Clone
    + Add<Output = Self>
    + AddAssign
    + Sub<Output = Self>
    + SubAssign
    + Neg<Output = Self>
    + Mul<Output = Self>
    + MulAssign
    + Debug
{
}

pub trait AbstractCyclicPolynomial<F: TwoAdicField,const N:usize>:
    AbstractPolynomial<F> + Sized
{}

pub trait AbstractPolynomialCoefficients<F: AbstractField,APE:AbstractPolynomialEvaluations<F,Self>>:
    AbstractPolynomial<F>
{
    fn eval(&self, x: F) -> F;
    fn fft(&self)->APE;
}

pub trait AbstractPolynomialEvaluations<F:AbstractField,APC:AbstractPolynomialCoefficients<F,Self>>:
    AbstractPolynomial<F>
{
    fn ifft(&self)->APC;
}
