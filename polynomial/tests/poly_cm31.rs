#[cfg(test)]
mod tests_poly_cm31 {
    use p3_mersenne_31::{Mersenne31Complex, Mersenne31};
    use rand::Rng;
    use p3_field::AbstractField;
    use polynomial::{AbstractPolynomialCoefficients,AbstractPolynomialEvaluations,AbstractCyclicPolynomial,AbstractPolynomial};
    type F = Mersenne31Complex<Mersenne31>;
    type PC = polynomial::coeffs::CyclicPolynomialCoefficients<F,8>;
    type PE = polynomial::evals::CyclicPolynomialEvaluations<F,8>;
    #[test]
    fn test_coeff_add_sub(){
        let mut rng = rand::thread_rng();
        let aa:PC = rng.gen();
        let bb:PC = rng.gen();
        let cc = aa+bb;
        let dd = cc-bb;
        assert_eq!(aa,dd);
    }
    #[test]
    fn test_fft_mul_random_half(){
        let mut rng = rand::thread_rng();
        let mut aa:PC = rng.gen();
        let mut bb:PC = rng.gen();
        for i in 4..8 {
            aa.vals[i] = F::ZERO;
            bb.vals[i] = F::ZERO;
        }
        let cc = aa*bb;
        let dd = (aa.fft()*bb.fft()).ifft();
        assert_eq!(cc,dd);
    }
    #[test]
    fn test_fft_mul_random_full(){
        let mut rng = rand::thread_rng();
        let aa:PC = rng.gen();
        let bb:PC = rng.gen();
        let cc = aa*bb;
        let dd = (aa.fft()*bb.fft()).ifft();
        assert_eq!(cc,dd);
    }
    #[test]
    fn test_dft_mul_random_half(){
        let mut rng = rand::thread_rng();
        let mut aa:PC = rng.gen();
        let mut bb:PC = rng.gen();
        for i in 4..8 {
            aa.vals[i] = F::ZERO;
            bb.vals[i] = F::ZERO;
        }
        let cc = aa*bb;
        let dd = (aa.dft()*bb.dft()).idft();
        assert_eq!(cc,dd);
    }
}