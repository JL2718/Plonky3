use std::usize;
use core::ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Sub, SubAssign};
use core::fmt::Debug;

use p3_field::{TwoAdicField,AbstractField};

use crate::evals::CyclicPolynomialEvaluations;
use crate::fft::fft;
use crate::{AbstractPolynomial,AbstractCyclicPolynomial,AbstractPolynomialCoefficients,AbstractPolynomialEvaluations};

#[derive(PartialEq, Eq, Hash, Copy, Clone, Debug)]
pub struct CyclicPolynomialCoefficients<F: TwoAdicField,const N:usize>
{
    vals: [F;N]
}

impl <F: TwoAdicField,const N:usize> CyclicPolynomialCoefficients<F,N> {
    pub fn new(vals: [F;N]) -> Self {
        Self { vals }
    }
    pub fn from_vec(vals:Vec<F>) -> Self {
        let mut result = Self::default();
        for (i,v) in vals.iter().enumerate() {
            result.vals[i%N] += *v;
        }
        return result;
    }
}
impl <F: TwoAdicField,const N:usize> Default for CyclicPolynomialCoefficients<F,N> {
    fn default() -> Self {
        Self::new([F::ZERO;N])
    }
}
/*
impl <F: TwoAdicField,const N:usize> Clone for CyclicPolynomialCoefficients<F,N> {
    fn clone(&self) -> Self {
        Self::new(self.vals.clone())
    }
}
impl <F: TwoAdicField,const N:usize> Debug for CyclicPolynomialCoefficients<F,N> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut first = true;
        for i in 0..N {
            if !first {
                write!(f, " + ")?;
            }
            first = false;
            write!(f, "{}", self.vals[i])?;
            if i > 0 {
                write!(f, "x^{}", i)?;
            }
        }
        Ok(())
    }
} 
*/
impl <F: TwoAdicField,const N:usize> AddAssign for CyclicPolynomialCoefficients<F,N> {
    fn add_assign(&mut self, rhs: Self) {
        for i in 0..N {
            self.vals[i] += rhs.vals[i];
        }
    }
}
impl <F: TwoAdicField,const N:usize> Add for CyclicPolynomialCoefficients<F,N> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        let mut result = self.clone();
        result += rhs;
        result
    }
}
impl <F: TwoAdicField,const N:usize> SubAssign for CyclicPolynomialCoefficients<F,N> {
    fn sub_assign(&mut self, rhs: Self) {
        for i in 0..N {
            self.vals[i] -= rhs.vals[i];
        }
    }
}
impl <F: TwoAdicField,const N:usize> Sub for CyclicPolynomialCoefficients<F,N> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        let mut result = self.clone();
        result -= rhs;
        result
    }
}
impl <F: TwoAdicField,const N:usize> Neg for CyclicPolynomialCoefficients<F,N> {
    type Output = Self;
    fn neg(self) -> Self {
        let mut result = self.clone();
        for i in 0..N {
            result.vals[i] = -result.vals[i];
        }
        result
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
impl <F: TwoAdicField,const N:usize> AbstractPolynomial<F> for CyclicPolynomialCoefficients<F,N> {
}

impl <F: TwoAdicField,const N:usize> AbstractPolynomialCoefficients<F,CyclicPolynomialEvaluations<F,N>> for CyclicPolynomialCoefficients<F,N> {
    fn eval(&self, x: F) -> F {
        let mut result = F::ZERO;
        let mut x_pow = F::ONE;
        for i in 0..N {
            result += self.vals[i] * x_pow;
            x_pow *= x;
        }
        return result;
    }
    fn fft(&self) -> CyclicPolynomialEvaluations<F,N> {
        let mut vals = self.vals.clone();
        fft::<F,N,true>(&mut vals);
        let result = CyclicPolynomialEvaluations::new(vals);
        return result;
    }

}
