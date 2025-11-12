//! Data structures that represent detected resonance systems.

use crate::core::atom::AtomId;
use crate::core::bond::BondId;

/// Connected conjugated component identified by the resonance detector.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ResonanceSystem {
    /// Stable atom identifiers that form the resonance system.
    pub atoms: Vec<AtomId>,
    /// Stable bond identifiers that form the resonance system.
    pub bonds: Vec<BondId>,
}

impl ResonanceSystem {
    /// Creates a new resonance system while de-duplicating inputs.
    pub fn new(mut atoms: Vec<AtomId>, mut bonds: Vec<BondId>) -> Self {
        atoms.sort_unstable();
        atoms.dedup();
        bonds.sort_unstable();
        bonds.dedup();
        Self { atoms, bonds }
    }
}
