pub type BondId = usize;

#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub enum BondOrder {
    Single,
    Double,
    Triple,
    Aromatic,
}

impl BondOrder {
    pub fn multiplicity(self) -> u8 {
        match self {
            BondOrder::Single => 1,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            BondOrder::Aromatic => 1, // Placeholder value; effective multiplicity is handled after Kekulization.
        }
    }
}
