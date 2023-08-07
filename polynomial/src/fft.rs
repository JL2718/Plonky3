use p3_field::TwoAdicField;

pub fn fft<F: TwoAdicField, const N: usize>(vals: &mut [F; N]) {
    let root = F::primitive_root_of_unity(N.trailing_zeros() as usize);
    rpermute(vals);
    rfft(vals, root);
    //rfft_dit(vals,root);
}
pub fn ifft<F: TwoAdicField, const N: usize>(vals: &mut [F; N]) {
    let root = F::primitive_root_of_unity(N.trailing_zeros() as usize).inverse();
    rpermute(vals);
    rfft(vals, root);
    //rfft_dif(vals,root);
    vals.iter_mut()
        .for_each(|v| *v *= F::from_canonical_usize(N).inverse());
}

pub fn _fft<F: TwoAdicField, const N: usize>(vals: &mut [F; N], root: F) {
    /*
    In-place FFT
    Cooley-Tukey FFT Algorithm on the input 'vals' using the root of unity 'root'
    */
    assert!(N.is_power_of_two());
    // root: 2^Nth root of unity or inverse root if inverse
    debug_assert!(root.exp_power_of_2(N.trailing_zeros() as usize) == F::ONE);
    // rr: sequence of root squares from {2^-1, 2^-2, 2^-4, ..., root=2^-N}
    let rr: Vec<F> = (0..N.trailing_zeros())
        .scan(root, |ri, _| {
            let ret = *ri;
            *ri *= *ri;
            Some(ret)
        })
        .collect();
    // Cooley-Tukey FFT Algorithm (in-place)
    for (i, r) in rr.iter().rev().enumerate() {
        for j in (0..N).step_by(1 << (i + 1)) {
            let mut s = F::ONE;
            for k in j..j + (1 << i) {
                let u = vals[k];
                let v = vals[k + (1 << i)] * s;
                vals[k] = u + v;
                vals[k + (1 << i)] = u - v;
                s *= *r;
            }
        }
    }
}

pub fn rfft<F: TwoAdicField>(vals: &mut [F], root: F) {
    /*
    in-place FFT DIT recursion without re-ordering transformation
    */
    let n = vals.len();
    assert!(n.is_power_of_two());
    // base case
    if n == 1 {
        assert!(root == F::ONE);
        return;
    }
    // split problem into even and odd halves
    let (aa, bb) = vals.split_at_mut(n / 2);
    // recurse on halves
    let rootsq = root * root;
    rfft(aa, rootsq);
    rfft(bb, rootsq);
    // compute phases
    let phases = (0..n / 2).scan(F::ONE, |ri, _| {
        let ret = *ri;
        *ri *= root;
        Some(ret)
    });
    // combine halves with phase transformation
    aa.iter_mut()
        .zip(bb.iter_mut())
        .zip(phases)
        .for_each(|((a, b), phase)| {
            let t = *b * phase;
            (*a, *b) = (*a + t, *a - t);
        });
}

pub fn rfft_dit<F: TwoAdicField>(vals: &mut [F], root: F) {
    /*
    in-place FFT DIT recursion
    */
    let n = vals.len();
    assert!(n.is_power_of_two());
    // base case
    if n == 1 {
        assert!(root == F::ONE);
        return;
    }
    // split problem into even and odd halves
    let (even, odd) = vals.split_at_mut(n / 2);
    // de-interleave evens and odds
    skipswap(even, odd);
    // recurse on halves
    let rootsq = root * root;
    rfft_dit(even, rootsq);
    rfft_dit(odd, rootsq);
    // compute phases
    let phases = (0..n / 2).scan(F::ONE, |ri, _| {
        let ret = *ri;
        *ri *= root;
        Some(ret)
    });
    // compute phase transformation on halves
    even.iter_mut()
        .zip(odd.iter_mut())
        .zip(phases)
        .for_each(|((e, o), phase)| {
            let t = *o * phase;
            (*e, *o) = (*e + t, *e - t);
        });
}

pub fn rfft_dif<F: TwoAdicField>(vals: &mut [F], root: F) {
    /*
    in-place FFT DIF recursion
    */
    let n = vals.len();
    assert!(n.is_power_of_two());
    assert!(root.exp_power_of_2(n.trailing_zeros() as usize) == F::ONE);
    // base case
    if n == 1 {
        assert!(root == F::ONE);
        return;
    }
    // split problem into even and odd halves
    let (lo, hi) = vals.split_at_mut(n / 2);
    // compute phases
    let phases = (0..n / 2).scan(F::ONE, |ri, _| {
        let ret = *ri;
        *ri *= root;
        Some(ret)
    });
    // apply phase transformation on halves
    lo.iter_mut()
        .zip(hi.iter_mut())
        .zip(phases)
        .for_each(|((l, h), phase)| {
            (*l, *h) = (*l + *h, (*l - *h) * phase);
        });
    // recurse on halves
    let rootsq = root * root;
    rfft_dif(lo, rootsq); // root not modified
    rfft_dif(hi, rootsq); // root modified
                          // interleave frequencies
    skipswap(lo, hi);
}

fn skipswap<F: TwoAdicField>(even: &mut [F], odd: &mut [F]) {
    even.iter_mut()
        .skip(1)
        .step_by(2)
        .zip(odd.iter_mut().step_by(2))
        .for_each(|(e, o)| std::mem::swap(e, o));
}

pub fn permute<F: TwoAdicField, const N: usize>(vals: &mut [F; N]) {
    /*
       Bit-reversal permutation
       Rearrange the input 'vals' in bit-reversal order
       This implementation is operation-efficient, in-place, but not cache efficient
    */
    assert!(N.is_power_of_two());
    for i in 0..N {
        let j = i.reverse_bits() >> (usize::BITS - N.trailing_zeros());
        if i < j {
            vals.swap(i, j);
        }
    }
}

pub fn rpermute<F: TwoAdicField>(vals: &mut [F]) {
    /*
    Recursive bit-reversal permutation
    This implementation is in-place, cache efficient, but not operation-efficient
    */

    let n = vals.len();
    assert!(n.is_power_of_two());
    // base case
    if n == 2 {
        return;
    }
    // split problem into even and odd halves
    let (even, odd) = vals.split_at_mut(n / 2);
    // recurse on halves
    rpermute(even);
    rpermute(odd);
    // exchange middle quarters
    even.iter_mut()
        .skip(n / 4)
        .zip(odd.iter_mut())
        .for_each(|(e, o)| {
            std::mem::swap(e, o);
        });
    interleave(even);
    interleave(odd);
}

pub fn interleave<F: TwoAdicField>(val: &mut [F]) {
    // recursive interleave of top and bottom half
    let n = val.len();
    if n == 2 {
        return;
    }
    let (top, bot) = val.split_at_mut(n / 2);
    interleave(top);
    interleave(bot);
    top.iter_mut()
        .skip(1)
        .zip(bot.iter_mut())
        .step_by(2)
        .for_each(|(t, b)| {
            std::mem::swap(t, b);
        });
}

#[cfg(test)]
mod test {
    use p3_field::AbstractField;

    use super::*;
    use crate::dft::{dft, idft};
    type F = crate::fp17::Fp17;
    use rand::Rng;
    use unroll::unroll_for_loops;
    #[test]
    #[unroll_for_loops]
    fn test_fft_all_simple() {
        for p2 in 1..4 {
            let mut vals = [F::ZERO; 1 << p2];
            vals[0] = F::ONE;
            let vals_0 = vals.clone();
            let vals_dft = dft(vals);
            fft(&mut vals);
            assert_eq!(vals, [F::ONE; 1 << p2]);
            assert_eq!(vals, vals_dft);
            let vals_idft = idft(vals);
            ifft(&mut vals);
            assert_eq!(vals, vals_0);
            assert_eq!(vals, vals_idft);
        }
    }
    #[test]
    #[unroll_for_loops]
    fn test_fft_all_random() {
        for p2 in 1..4 {
            let mut rng = rand::thread_rng();
            let mut vals: [F; 1 << p2] = rng.gen();
            let vals_0 = vals.clone();
            let vals_dft = dft(vals);
            fft(&mut vals);
            assert_eq!(vals, vals_dft);
            let vals_idft = idft(vals);
            ifft(&mut vals);
            assert_eq!(vals, vals_0);
            assert_eq!(vals, vals_idft);
        }
    }
    #[test]
    #[unroll_for_loops]
    fn test_dft_rfft() {
        for p2 in 1..4 {
            let mut vals = [F::ZERO; 1 << p2];
            vals[0] = F::ONE;
            let vals2 = dft(vals);
            fft(&mut vals);
            assert_eq!(vals, vals2);
        }
    }
    #[test]
    #[unroll_for_loops]
    fn test_permute() {
        for p2 in 1..8 {
            let mut vals: [F; 1 << p2] = core::array::from_fn(|i| F::from_canonical_usize(i));
            let mut vals2 = vals.clone();
            permute(&mut vals);
            rpermute(&mut vals2);
            assert_eq!(vals, vals2);
        }
    }
    #[test]
    fn test_permute_n() {
        let mut vals: [F; 8] = core::array::from_fn(|i| F::from_canonical_usize(i));
        let mut vals2 = vals.clone();
        permute(&mut vals);
        rpermute(&mut vals2);
        assert_eq!(vals, vals2);
    }
    #[test]
    #[unroll_for_loops]
    fn test_fft_ifft_random() {
        for p2 in 1..4 {
            for _ in 0..100 {
                let mut rng = rand::thread_rng();
                let mut aa = [F::default(); 1 << p2];
                aa.iter_mut().for_each(|a| *a = F::new(rng.gen()));
                let aa_0 = aa;
                fft(&mut aa);
                ifft(&mut aa);
                assert_eq!(aa, aa_0);
            }
        }
    }
}
