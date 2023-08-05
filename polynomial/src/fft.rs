use p3_field::TwoAdicField;

pub fn _fft<F: TwoAdicField, const N: usize, const INV: bool>(vals: &mut [F; N]) {
    /*
    In-place FFT
    Cooley-Tukey FFT Algorithm on the input 'vals' using the root of unity 'root'
    TODO: look into eliminating the bit-reversal requirement
    */
    assert!(N.is_power_of_two());
    // root: 2^Nth root of unity or inverse root if inverse
    let root: F = if INV {
        F::primitive_root_of_unity(N.trailing_zeros() as usize).inverse()
    } else {
        F::primitive_root_of_unity(N.trailing_zeros() as usize)
    };
    debug_assert!(root.exp_power_of_2(N.trailing_zeros() as usize) == F::ONE);
    // rr: sequence of root squares from {2^-1, 2^-2, 2^-4, ..., root=2^-N}
    let rr: Vec<F> = (0..N.trailing_zeros())
        .scan(root, |ri, _| {
            let ret = *ri;
            *ri *= *ri;
            Some(ret)
        })
        .collect();
    // bit-reversal permutation
    permute(vals);
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
    // divide by N if inverse
    if INV {
        let inv = F::TWO
            .exp_u64(N.trailing_zeros().try_into().unwrap())
            .inverse();
        for v in vals {
            *v *= inv;
        }
    }
}

pub fn rfft<F: TwoAdicField>(vals: &mut [F], root: F) {
    /*
    in-place FFT recursion, includes bit-reversal permutation
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
    // exchange values for bit revesal
    even.iter_mut()
        .skip(1)
        .zip(odd.iter_mut())
        .step_by(2)
        .for_each(|(e, o)| {
            std::mem::swap(e, o);
        });
    // recurse on halves
    let rootsq = root * root;
    rfft(even, rootsq);
    rfft(odd, rootsq);
    // compute phases
    let phases = (0..n / 2).scan(F::ONE, |ri, _| {
        let ret = *ri;
        *ri *= root;
        Some(ret)
    });
    // combine halves with phase transformation
    even.iter_mut()
        .zip(odd)
        .zip(phases)
        .for_each(|((e, o), phase)| {
            let t = *o * phase;
            *o = *e - t;
            *e += t;
        });
}

pub fn fft<F: TwoAdicField, const N: usize>(vals: &mut [F; N]) {
    let root = F::primitive_root_of_unity(N.trailing_zeros() as usize);
    rfft(vals, root);
}
pub fn ifft<F: TwoAdicField, const N: usize>(vals: &mut [F; N]) {
    let root = F::primitive_root_of_unity(N.trailing_zeros() as usize).inverse();
    rfft(vals, root);
    vals.iter_mut()
        .for_each(|v| *v *= F::from_canonical_usize(N).inverse());
}

pub fn permute<F: TwoAdicField, const N: usize>(vals: &mut [F; N]) {
    /*
       Bit-reversal permutation
       Rearrange the input 'vals' in bit-reversal order
    */
    assert!(N.is_power_of_two());
    for i in 0..N {
        let j = i.reverse_bits() >> (usize::BITS - N.trailing_zeros());
        if i < j {
            vals.swap(i, j);
        }
    }
}
