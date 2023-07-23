#![feature(generic_const_exprs)]
#![feature(associated_const_equality)]

use p3_field::{TwoAdicField,AbstractField};

pub fn add(left: usize, right: usize) -> usize {
    left + right
}

pub fn dft
<F:TwoAdicField,const LN:usize,const INV:bool>
(vals:[F;1<<LN])->[F;1<<LN]
{
    let root:F = (|r:F|if INV {r.inverse()} else {r})(
        F::power_of_two_generator()
        .exp_power_of_2(F::TWO_ADICITY - LN)
    );
    debug_assert!(root.exp_power_of_2(LN) == F::ONE);
    
    //roots = (0..1<<N).scan(F::ONE,|ri,i|{ri*=root;ri}).map(|ri|{(0..1<<N).scan(F::ONE,|rij,j|{rij*=ri;rij})});
    
    let mut ret:[F;1<<LN] = [F::default();1<<LN];
    let mut ri = F::ONE; // root^i
    for i in 0..(1<<LN) {
        let mut rij = F::ONE;// root^(i*j)
        let mut sum = F::ZERO;
        for j in 0..(1<<LN) {
            sum += vals[j] * rij;
            rij *= ri;
        }
        ri *= root;
        debug_assert!(rij==F::ONE);
        ret[i] = sum * (if INV {F::TWO.exp_u64(LN.try_into().unwrap()).inverse()} else {F::ONE});
    }
    debug_assert!(ri==F::ONE);
    ret
}

pub fn fft<F:TwoAdicField,const N:usize>(mut vals:&[F;N])-> &[F;N]{
    vals
}

pub fn ifft<F:TwoAdicField, const N:usize>(mut vals:&[F;N])->&[F;N]{
    return vals;
}

#[cfg(test)]
mod tests_mersenne {
    use p3_mersenne_31::{Mersenne31,Mersenne31Complex};
    use rand::Rng;
    type B = Mersenne31;
    type F = Mersenne31Complex<Mersenne31>;

    use super::*;
    #[test]
    fn it_works() {
        let result = add(2, 2);
        assert_eq!(result, 4);
    }
    #[test]
    fn test_fft_inversion(){
        /*
        Tests whether fft and ifft are inverse transformations
        */
        type B = Mersenne31;
        type F = Mersenne31Complex<B>;
        const N:usize = 8;
        let mut rng = rand::thread_rng();
        let mut aa = [F::default();N];
        aa.iter_mut().for_each(|a| *a=F::new_real(rng.gen::<B>())); 
        let aa_0 = aa.clone();
        fft(&aa);
        ifft(&aa);
        assert_eq!(aa,aa_0);
    }

    #[test]
    fn test_dft_inversion(){
        const N:usize = 8;
        let mut rng = rand::thread_rng();
        let mut aa = [F::default();N];
        aa.iter_mut().for_each(|a| *a=F::new_real(rng.gen::<B>())); 
        //let mut aa = [F::ZERO;N];
        //aa[0] = F::ONE;
        let tt = dft::<F,3,false>(aa);
        let aa_t = dft::<F,3,true>(tt);
        assert_eq!(aa,aa_t);
    }

    #[test]
    fn test_dft_polymul(){
        const N:usize=8;
        let mut rng = rand::thread_rng();
        let mut aa = [F::default();N];
        let mut bb = [F::default();N];
        aa.iter_mut().for_each(|v| *v=F::new_real(rng.gen::<B>())); 
        bb.iter_mut().for_each(|v| *v=F::new_real(rng.gen::<B>())); 
        // TODO: finish this, maybe need a full polynomial implementation
    }

}