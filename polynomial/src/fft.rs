
#![feature(generic_const_exprs)]
#![feature(associated_const_equality)]

use p3_field::{TwoAdicField};

pub fn fft
<F:TwoAdicField,const N:usize,const INV:bool>
(vals:&mut [F;N])
{
    /*
    In-place FFT
    Cooley-Tukey FFT Algorithm on the input 'vals' using the root of unity 'root'
    no bit-reveral is performed
     */
    debug_assert!(N.is_power_of_two());
    // get 2^Nth root of unity or inverse root if inverse
    let root:F = (|r:F|if INV {r.inverse()} else {r})(
        F::power_of_two_generator()
        .exp_power_of_2((F::TWO_ADICITY as isize  - N.trailing_zeros() as isize) as usize)
    );
    debug_assert!(root.exp_power_of_2(N.trailing_zeros().try_into().unwrap()) == F::ONE);
    // Cooley-Tukey FFT Algorithm 
    let mut w = root;
    for i in 0..N.trailing_zeros(){
        for j in (0..N).step_by(1<<(i+1)) {
            for k in j..j + (1<<i) {
                let u = vals[k];
                let v = vals[k + (1<<i)] * w;
                vals[k] = u + v;
                vals[k + (1<<i)] = u - v;
            }
        }
        w *= w;        
    }
    // divide by N if inverse
    if INV {
        let inv = F::TWO.exp_u64(N.trailing_zeros().try_into().unwrap()).inverse();
        for i in 0..N {
            vals[i] *= inv;
        }
    }
}

pub fn bit_reverse
<F:TwoAdicField,const N:usize>
(vals:& mut [F;N]){
    /*
        Bit-reversal permutation
        Rearrange the input 'vals' in bit-reversal order
     */
    debug_assert!(N.is_power_of_two());
    for i in 0..N {
        let j = i.reverse_bits() >> (usize::BITS - N.trailing_zeros());
        if i < j {
            vals.swap(i,j);
        }
    }
}

#[cfg(test)]
mod tests_mersenne {
    use p3_field::AbstractField;
    use p3_mersenne_31::{Mersenne31,Mersenne31Complex};
    use rand::Rng;
    type B = Mersenne31;
    type F = Mersenne31Complex<Mersenne31>;
    use crate::fft::*;
    use crate::dft::dft;
    #[test]
    fn test_fft_simple(){
        const N:usize = 8;
        let mut aa = [F::ZERO;N];
        aa[0] = F::ONE;
        let aa_0 = aa.clone();
        fft::<F,8,false>(& mut aa);
        fft::<F,8,true>(& mut aa);
        assert_eq!(aa,aa_0);
    }
    #[test]
    fn test_fft_random(){
        const N:usize = 8;
        let mut rng = rand::thread_rng();
        let mut aa = [F::default();N];
        aa.iter_mut().for_each(|a| *a=F::new_real(rng.gen::<B>())); 
        let aa_0 = aa.clone();
        fft::<F,8,false>(& mut aa);
        fft::<F,8,true>(& mut aa);
        assert_eq!(aa,aa_0);
    }
    #[test]
    fn test_fft_dft_simple(){
        const N:usize = 8;
        let mut aa = [F::ZERO;N];
        aa[0] = F::ONE;
        let aa_0 = aa.clone();
        bit_reverse(& mut aa);
        fft::<F,8,false>(& mut aa);
        let aa_1 = aa.clone();
        let aa_2 = dft::<F,8,false>(aa_0);
        assert_eq!(aa_1,aa_2);
    }
    #[test]
    fn test_fft_dft_random(){
        const N:usize = 8;
        let mut rng = rand::thread_rng();
        let mut aa = [F::default();N];
        aa.iter_mut().for_each(|a| *a=F::new_real(rng.gen::<B>())); 
        let aa_0 = aa.clone();
        bit_reverse(& mut aa);
        fft::<F,8,false>(& mut aa);
        let aa_1 = aa.clone();
        let aa_2 = dft::<F,8,false>(aa_0);
        assert_eq!(aa_1,aa_2);
    }
}