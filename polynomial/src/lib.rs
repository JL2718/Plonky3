#![feature(generic_const_exprs)]

use std::usize;
use std::iter::zip;

use p3_field::{AbstractExtensionField, AbstractField, AbstractionOf, Field, TwoAdicField};


pub trait FFTField: AbstractField + TwoAdicField
{
}

pub struct CoefficientBasis<F: FFTField,const N:usize> {
    coeffs: [F;N],
}
impl <F: FFTField,const N:usize> CoefficientBasis<F,N> {
    fn fft(&self,root:F) -> EvaluationBasis<F,N> {
        let mut evals = [F::ZERO;N];
        let root:F = F::power_of_two_generator().exp_power_of_2(F::TWO_ADICITY - N);
        EvaluationBasis { evals }
    }
}
 pub struct EvaluationBasis<F: FFTField,const N:usize> {
    evals: [F;N],
}
impl <F: FFTField,const N:usize> EvaluationBasis<F,N> {
    fn ifft(&self) -> CoefficientBasis<F,N> {
        let mut coeffs = [F::ZERO;N];
        let root:F = F::power_of_two_generator().exp_power_of_2(F::TWO_ADICITY - N);
        CoefficientBasis{ coeffs }
    }
}

pub fn fft<F:FFTField,const P2:usize>(input:[F;1<<P2]) -> [F;1<<P2] where [F]: Sized {
    let root:F = F::power_of_two_generator().exp_power_of_2(F::TWO_ADICITY - P2);
    let mut output:[F;1<<P2] = [F::ZERO;1<<P2];
    if P2==0 {
        output[0] = input[0];
        return output;
    } else {
        let even_in:[F;1<<(P2-1)] = input.iter().step_by(2).copied().collect::<[F;1<<(P2-1)]>();
        let even:[F;1<<(P2-1)] = fft::<F,{P2-1}>(
            input.iter().step_by(2)
            .collect::<[F;1<<(P2-1)]>()
        ); 

        let odd:[F;1<<(P2-1)] = fft::<F,{P2-1}>(
            input.iter().skip(1).step_by(2)
        ).iter().scan(F::ONE,|w:&F,o:F| -> F {
            let v = *w*o;
            *w *= root;
            v
        }).collect::<[F;1<<(P2-1)]>();
        output[..(1<<(P2-1))] = zip(even.iter(),odd.iter()).map(|(e,o)|e+o).collect::<[F;1<<(P2-1)]>();
        output[(1<<(P2-1))..] = zip(even.iter(),odd.iter()).map(|(e,o)|e-o).collect::<[F;1<<(P2-1)]>();
        return output;
    } 
}
