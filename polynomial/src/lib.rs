#![feature(generic_const_exprs)]
#![feature(associated_const_equality)]

mod dft;
mod fft;

use std::usize;

use p3_field::{AbstractField};

pub trait AbstractPolynomial<BaseField: AbstractField>:
{}

pub struct CoefficientBasis<F: AbstractField,const N:usize> {
    coeffs: [F;N],
}
pub struct EvaluationBasis<F: AbstractField,const N:usize> {
    evals: [F;N],
}
