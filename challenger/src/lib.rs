//! Utilities for generating Fiat-Shamir challenges based on an IOP's transcript.

#![no_std]

extern crate alloc;

mod duplex_challenger;
mod hash_challenger;

pub use duplex_challenger::*;
pub use hash_challenger::*;

use alloc::vec::Vec;
use p3_field::{Field, FieldExtension};

/// Observes prover messages during an IOP, and generates Fiat-Shamir challenges in response.
pub trait Challenger<F: Field> {
    fn observe_element(&mut self, element: F);

    fn observe_elements(&mut self, elements: &[F]) {
        for &elt in elements {
            self.observe_element(elt);
        }
    }

    fn observe_ext_element<FE: FieldExtension<Base = F>>(&mut self, ext: FE::Extension) {
        self.observe_elements(FE::as_base_slice(&ext));
    }

    fn random_element(&mut self) -> F;

    fn random_ext_element<FE: FieldExtension<Base = F>>(&mut self) {
        let vec = self.random_vec(FE::DEGREE);
        FE::from_base_slice(&vec);
    }

    fn random_array<const N: usize>(&mut self) -> [F; N] {
        core::array::from_fn(|_| self.random_element())
    }

    fn random_vec(&mut self, n: usize) -> Vec<F> {
        (0..n).map(|_| self.random_element()).collect()
    }
}
