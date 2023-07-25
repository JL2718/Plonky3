use p3_field::{TwoAdicField};

pub fn dft
<F:TwoAdicField,const N:usize,const INV:bool>
(vals:[F;N])->[F;N]
{
    let root:F = (|r:F|if INV {r.inverse()} else {r})(
        F::power_of_two_generator()
        .exp_power_of_2((F::TWO_ADICITY as isize  - N.trailing_zeros() as isize) as usize)
        //.exp_power_of_2(F::TWO_ADICITY - N.trailing_zeros())
    );
    debug_assert!(root.exp_power_of_2(N.trailing_zeros().try_into().unwrap()) == F::ONE);
    
    //roots = (0..1<<N).scan(F::ONE,|ri,i|{ri*=root;ri}).map(|ri|{(0..1<<N).scan(F::ONE,|rij,j|{rij*=ri;rij})});
    
    let mut ret:[F;N] = [F::default();N];
    let mut ri = F::ONE; // root^i
    for i in 0..(N) {
        let mut rij = F::ONE;// root^(i*j)
        let mut sum = F::ZERO;
        for j in 0..(N) {
            sum += vals[j] * rij;
            rij *= ri;
        }
        ri *= root;
        debug_assert!(rij==F::ONE);
        ret[i] = sum * (if INV {F::TWO.exp_u64(N.trailing_zeros().try_into().unwrap()).inverse()} else {F::ONE});
    }
    debug_assert!(ri==F::ONE);
    ret
}


#[cfg(test)]
mod tests_mersenne {
    use p3_mersenne_31::{Mersenne31,Mersenne31Complex};
    use rand::Rng;
    use crate::dft::*;

    type B = Mersenne31;
    type F = Mersenne31Complex<Mersenne31>;

    #[test]
    fn test_dft_inversion(){
        const N:usize = 8;
        let mut rng = rand::thread_rng();
        let mut aa = [F::default();N];
        aa.iter_mut().for_each(|a| *a=F::new_real(rng.gen::<B>())); 
        //let mut aa = [F::ZERO;N];
        //aa[0] = F::ONE;
        let tt = dft::<F,8,false>(aa);
        let aa_t = dft::<F,8,true>(tt);
        assert_eq!(aa,aa_t);
    }
}