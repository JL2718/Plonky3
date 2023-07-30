use p3_field::TwoAdicField;

pub fn fft
<F:TwoAdicField,const N:usize,const INV:bool>
(vals:&mut [F;N])
{
    /*
    In-place FFT
    Cooley-Tukey FFT Algorithm on the input 'vals' using the root of unity 'root'
    TODO: look into eliminating the bit-reversal requirement
    */
    static_assert_power_of_two::<N>();
    // root: 2^Nth root of unity or inverse root if inverse
    let root:F = (|r:F|if INV {r.inverse()} else {r})(
        F::power_of_two_generator().exp_power_of_2((F::TWO_ADICITY as isize  - N.trailing_zeros() as isize) as usize)
    );
    debug_assert!(root.exp_power_of_2(N.trailing_zeros() as usize) == F::ONE);
    // rr: sequence of root squares from {2^-1, 2^-2, 2^-4, ..., root=2^-N}
    let rr:Vec<F>=(0..N.trailing_zeros()).scan(root,|ri,_|{*ri*=*ri;Some(*ri)}).collect();
    // bit-reversal permutation
    bit_reverse_permutation(vals);
    // Cooley-Tukey FFT Algorithm (in-place)
    //for i in 0..(N.trailing_zeros() as usize) {
        //let r = root.exp_power_of_2(i);
    for (i,r) in rr.iter().rev().enumerate(){
        for j in (0..N).step_by(1<<(i+1)) {
            let mut s = F::ONE;
            for k in j..j + (1<<i) {
                let u = vals[k];
                let v = vals[k + (1<<i)] * s;
                vals[k] = u + v;
                vals[k + (1<<i)] = u - v;
                s *= *r;
            }
        }
    }
    // divide by N if inverse
    if INV {
        let inv = F::TWO.exp_u64(N.trailing_zeros() as u64).inverse();
        for i in 0..N {
            vals[i] *= inv;
        }
    }
}

pub fn root_of_unity<F:TwoAdicField,const N:usize>()->F{
    static_assert_power_of_two::<N>();
    F::power_of_two_generator().exp_power_of_2((F::TWO_ADICITY as isize  - N.trailing_zeros() as isize) as usize)
}

fn reverse_bits<const N:usize>(v:usize)->usize{
    return v.reverse_bits() >> (usize::BITS - N.trailing_zeros());
}

const fn static_assert_power_of_two<const N:usize>(){
    assert!(N.is_power_of_two(), "array size must be a power of two");
}

pub fn bit_reverse_permutation
<F:TwoAdicField,const N:usize>
(vals:& mut [F;N]){
    /*
        Bit-reversal permutation
        Rearrange the input 'vals' in bit-reversal order
     */
    static_assert_power_of_two::<N>();
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
        fft::<F,8,false>(& mut aa);
        let aa_1 = aa.clone();
        let aa_2 = dft::<F,8,false>(aa_0);
        assert_eq!(aa_1,aa_2);
    }
    #[test]
    fn test_ifft_idft_random(){
        const N:usize = 8;
        let mut rng = rand::thread_rng();
        let mut aa = [F::default();N];
        aa.iter_mut().for_each(|a| *a=F::new_real(rng.gen::<B>())); 
        let aa_0 = aa.clone();
        fft::<F,8,true>(& mut aa);
        let aa_1 = aa.clone();
        let aa_2 = dft::<F,8,true>(aa_0);
        assert_eq!(aa_1,aa_2);
    }
}