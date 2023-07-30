#![feature(generic_const_exprs)]
#![feature(associated_const_equality)]

mod dft;
mod fft;

use std::usize;
use core::ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Sub, SubAssign};
use core::fmt::Debug;

use p3_field::{AbstractField};


pub trait AbstractPolynomial<F: AbstractField,N:usize>:
    + Default
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

pub trait AbstractCyclicPolynomial<F: AbstractField,N:usize>:
    AbstractPolynomial<F,N> + Sized
{
    vals: [F;N]
    fn new(vals: [F;N]) -> Self {
        Self { vals }
    }
    fn coefficients(self) -> CyclicPolynomialCoefficients<F,N>;
    fn evaluations(self) -> CyclicPolynomialEvaluations<F,N>;
}

impl <F: AbstractField,N:usize> Default for AbstractCyclicPolynomial<F,N> {
    fn default() -> Self {
        Self::new([F::ZERO;N])
    }
}

impl<F: AbstractField,N:usize> Add for AbstractCyclicPolynomial<F,N> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        let mut result = self.clone();
        result += rhs;
        result
    }
}

impl<F: AbstractField,N:usize> AddAssign for AbstractCyclicPolynomial<F,N> {
    fn add_assign(&mut self, rhs: Self) {
        for i in 0..N {
            self.vals[i] += rhs.vals[i];
        }
    }
}

impl<F: AbstractField,N:usize> Sub for AbstractCyclicPolynomial<F,N> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        let mut result = self.clone();
        result -= rhs;
        result
    }
}

impl<F: AbstractField,N:usize> SubAssign for AbstractCyclicPolynomial<F,N> {
    fn sub_assign(&mut self, rhs: Self) {
        for i in 0..N {
            self.vals[i] -= rhs.vals[i];
        }
    }
}

impl<F: AbstractField,N:usize> Neg for AbstractCyclicPolynomial<F,N> {
    type Output = Self;

    fn neg(self) -> Self {
        let mut result = self.clone();
        for i in 0..N {
            result.vals[i] = -result.vals[i];
        }
        result
    }
}


#[derive(PartialEq, Eq, Hash, Copy, Clone, Debug, Default)]
pub struct CyclicPolynomialCoefficients<F: AbstractField,N:usize>:
    AbstractCyclicPolynomial<F,N>
{
    fn new(vals: !vec<F>) -> Self {
        let mut result = Self::default();
        for i,v in vals.iter().enumerate() {
            result.vals[i%N] += v;
        }
        return result;
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
        fft.fft(&mut vals);
        let result = CyclicPolynomialEvaluations::new(vals);
        return result;
    }

}

impl<F: AbstractField,N:usize> Mul for CyclicPolynomialCoefficients<F,N> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        let mut result = Self::zero();
        for i in 0..N {
            for j in 0..N {
                result.vals[(i + j) % N] += self.vals[i]*rhs.vals[j];
            }
        }
        result
    }
}

impl<F: AbstractField,N:usize> MulAssign for CyclicPolynomialCoefficients<F,N> {
    fn mul_assign(&mut self, rhs: Self) {
        let mut result = self.clone();
        *self = result*rhs;
    }
}

#[derive(PartialEq, Eq, Hash, Copy, Clone, Debug, Default)]
pub struct CyclicPolynomialEvaluations<F: AbstractField,N:usize>:
    AbstractCyclicPolynomial<F,N>
{
    fn evaluations(&self) -> &Self {
        return self;
    }

    fn coefficients(&self) -> CyclicPolynomialCoefficients<F,N> {
        let mut vals = self.vals.clone();
        fft.fft(&mut vals,inverse=true);
        let result = CyclicPolynomialCoefficients::new(vals);
        return result;
    }
}

impl<F: AbstractField,N:usize> Mul for CyclicPolynomialEvaluations<F,N> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        let mut result = Self::default();
        for i in 0..N {
            result.vals[i] = self.vals[i] * rhs.vals[i];
        }
        result
    }
}

impl<F: AbstractField,N:usize> MulAssign for CyclicPolynomialEvaluations<F,N> {
    fn mul_assign(&mut self, rhs: Self) {
        for i in 0..N {
            self.vals[i] *= rhs.vals[i];
        }
    }
}
