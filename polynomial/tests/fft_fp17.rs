#[cfg(test)]
mod tests_fft_fp17 {
    use p3_field::AbstractField;
    use rand::Rng;
    type F = polynomial::fp17::Fp17;
    use polynomial::dft::{dft, idft};
    use polynomial::fft::{fft, ifft};
    use unroll::unroll_for_loops;
    #[test]
    fn test_fft_ifft_simple() {

        const N: usize = 8;
        let mut aa = [F::ZERO; N];
        aa[0] = F::ONE;
        let aa_0 = aa;
        fft(&mut aa);
        ifft(&mut aa);
        assert_eq!(aa, aa_0);
    }
    #[test]
    #[unroll_for_loops]
    fn test_fft_ifft_random() {
        for p2 in 1..4 {
            let mut rng = rand::thread_rng();
            let mut aa = [F::default(); 1<<p2];
            aa.iter_mut().for_each(|a| *a = F::new(rng.gen()));
            let aa_0 = aa;
            fft(&mut aa);
            ifft(&mut aa);
            assert_eq!(aa, aa_0);
        }
    }
    #[test]
    fn test_fft_dft_simple_0() {
        const N: usize = 8;
        let mut aa = [F::ZERO; N];
        let aa_0 = dft(aa);
        fft(&mut aa);
        assert_eq!(aa, aa_0);
    }
    #[test]
    fn test_fft_dft_simple_1() {
        const N: usize = 8;
        let mut aa = [F::ZERO; N];
        aa[0] = F::ONE;
        let aa_0 = dft(aa);
        fft(&mut aa);
        assert_eq!(aa, aa_0);
    }
    #[test]
    fn test_fft_dft_simple_2() {
        const N: usize = 8;
        let mut aa = [F::ONE; N];
        let aa_0 = dft(aa);
        fft(&mut aa);
        assert_eq!(aa, aa_0);
    }
    #[test]
    fn test_fft_dft_random() {
        const N: usize = 8;
        let mut rng = rand::thread_rng();
        let mut aa = [F::default(); N];
        aa.iter_mut().for_each(|a| *a = F::new(rng.gen()));
        let aa_0 = dft(aa);
        fft(&mut aa);
        assert_eq!(aa, aa_0);
    }
    #[test]
    fn test_ifft_idft_random() {
        const N: usize = 8;
        let mut rng = rand::thread_rng();
        let mut aa = [F::default(); N];
        aa.iter_mut().for_each(|a| *a = F::new(rng.gen()));
        let aa_0 = idft(aa);
        ifft(&mut aa);
        assert_eq!(aa, aa_0);
    }
}
