use crate::core::atom::AtomId;
use crate::core::bond::BondId;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ResonanceSystem {
    pub atoms: Vec<AtomId>,
    pub bonds: Vec<BondId>,
}

impl ResonanceSystem {
    pub fn new(mut atoms: Vec<AtomId>, mut bonds: Vec<BondId>) -> Self {
        atoms.sort_unstable();
        atoms.dedup();
        bonds.sort_unstable();
        bonds.dedup();
        Self { atoms, bonds }
    }
}
