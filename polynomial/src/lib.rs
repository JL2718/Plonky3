#![feature(generic_const_exprs)]
#![feature(associated_const_equality)]

pub mod dft;
pub mod fft;
pub mod fp17;

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
{   
    fn new(vals: [F;N]) -> Self;
    fn default() -> Self;
    fn coefficients(self) -> CyclicPolynomialCoefficients<F,N>;
    fn evaluations(self) -> CyclicPolynomialEvaluations<F,N>;
    fn add(self, rhs: Self) -> Self {
        let mut result = self.clone();
        result += rhs;
        result
    }
    fn sub(self, rhs: Self) -> Self {
        let mut result = self.clone();
        result -= rhs;
        result
    }
}


//#[derive(PartialEq, Eq, Hash, Copy, Clone, Debug, Default)]
pub struct CyclicPolynomialCoefficients<F: TwoAdicField,const N:usize>
{
    vals: [F;N]
}

impl <F: TwoAdicField,const N:usize> CyclicPolynomialCoefficients<F,N> {
    fn new(vals: [F;N]) -> Self {
        Self { vals }
    }
    fn default() -> Self {
        Self::new([F::ZERO;N])
    }
    fn clone(&self) -> Self {
        Self::new(self.vals.clone())
    }
    fn from_vec(vals:Vec<F>) -> Self {
        let mut result = Self::default();
        for (i,v) in vals.iter().enumerate() {
            result.vals[i%N] += *v;
        }
        return result;
    }
    fn add_assign(&mut self, rhs: Self) {
        for i in 0..N {
            self.vals[i] += rhs.vals[i];
        }
    }
    fn sub_assign(&mut self, rhs: Self) {
        for i in 0..N {
            self.vals[i] -= rhs.vals[i];
        }
    }
    fn neg(self) -> Self {
        let mut result = self.clone();
        for i in 0..N {
            result.vals[i] = -result.vals[i];
        }
        result
    }

    fn eval(&self, x: F) -> F {
        let mut result = F::ZERO;
        let mut x_pow = F::ONE;
        for i in 0..N {
            result += self.vals[i] * x_pow;
            x_pow *= x;
        }
        return result;
    }

    fn coefficients(&self) -> &Self {
        return self;
    }

    fn evaluations(&self) -> CyclicPolynomialEvaluations<F,N> {
        let mut vals = self.vals.clone();
        fft::fft::<F,N,true>(&mut vals);
        let result = CyclicPolynomialEvaluations::new(vals);
        return result;
    }

}

impl<F: TwoAdicField,const N:usize> Mul for CyclicPolynomialCoefficients<F,N> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        let mut result = Self::default();
        for i in 0..N {
            for j in 0..N {
                result.vals[(i + j) % N] += self.vals[i]*rhs.vals[j];
            }
        }
        result
    }
}

impl<F: TwoAdicField,const N:usize> MulAssign for CyclicPolynomialCoefficients<F,N> {
    fn mul_assign(&mut self, rhs: Self) {
        *self = self.clone()*rhs;
    }
}

//#[derive(PartialEq, Eq, Hash, Copy, Clone, Debug, Default)]
pub struct CyclicPolynomialEvaluations<F: TwoAdicField,const N:usize>
{
    vals: [F;N],
}
impl <F: TwoAdicField,const N:usize> CyclicPolynomialEvaluations<F,N>
{
    fn vals(&self) -> &[F;N] {
        return &self.vals;
    }
    fn new(vals: [F;N]) -> Self {
        Self { vals }
    }
    fn default() -> Self {
        Self::new([F::ZERO;N])
    }
    fn clone(&self) -> Self {
        Self::new(self.vals.clone())
    }
    fn evaluations(&self) -> &Self {
        return self;
    }
    fn coefficients(&self) -> CyclicPolynomialCoefficients<F,N> {
        let mut vals = self.vals.clone();
        fft::fft::<F,N,false>(&mut vals);
        let result = CyclicPolynomialCoefficients::new(vals);
        return result;
    }
    fn add_assign(&mut self, rhs: Self) {
        for i in 0..N {
            self.vals[i] += rhs.vals[i];
        }
    }
    fn sub_assign(&mut self, rhs: Self) {
        for i in 0..N {
            self.vals[i] -= rhs.vals[i];
        }
    }
    fn neg(self) -> Self {
        let mut result = self.clone();
        for i in 0..N {
            result.vals[i] = -result.vals[i];
        }
        result
    }
}

impl<F: TwoAdicField,const N:usize> Mul for CyclicPolynomialEvaluations<F,N> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        let mut result = Self::default();
        for i in 0..N {
            result.vals[i] = self.vals[i] * rhs.vals[i];
        }
        result
    }
}

impl<F: TwoAdicField,const N:usize> MulAssign for CyclicPolynomialEvaluations<F,N> {
    fn mul_assign(&mut self, rhs: Self) {
        for i in 0..N {
            self.vals[i] *= rhs.vals[i];
        }
    }
}
