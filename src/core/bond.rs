//! Core bond primitives shared across the perception pipeline.
//!
//! The module defines `BondId` and `BondOrder`, which are reused by the graph
//! abstraction layer and every perception stage.

/// Unique identifier for a bond inside a molecular graph.
pub type BondId = usize;

/// Classification of a bond's order as used by the perception pipeline.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub enum BondOrder {
    /// Single bond containing one shared electron pair.
    Single,
    /// Double bond containing two shared electron pairs.
    Double,
    /// Triple bond containing three shared electron pairs.
    Triple,
    /// Aromatic bond flagged by the input or detected by perception stages.
    Aromatic,
}

impl BondOrder {
    /// Returns the valence contribution represented by this bond order.
    ///
    /// Aromatic bonds yield a multiplicity of 1 because the explicit electron
    /// counting is deferred to the Kekulé resonance model.
    pub fn multiplicity(self) -> u8 {
        match self {
            BondOrder::Single => 1,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            BondOrder::Aromatic => 1, // Placeholder value; effective multiplicity is handled after kekulization.
        }
    }
}
