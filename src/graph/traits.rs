//! Traits that describe the minimum molecular graph interface required by Pauling.
//!
//! User-defined graph types can participate in the perception pipeline by
//! implementing these views. The traits deliberately avoid mutability so that
//! perception can operate on borrowed data without copying.

use crate::core::atom::{AtomId, Element};
use crate::core::bond::{BondId, BondOrder};

/// Read-only view over an atom supplied by a user-defined molecular graph.
pub trait AtomView {
    /// Returns the stable identifier of this atom.
    fn id(&self) -> AtomId;

    /// Reports the chemical element associated with the atom.
    fn element(&self) -> Element;

    /// Returns the formal charge stored on the atom.
    fn formal_charge(&self) -> i8;
}

/// Read-only view over a bond supplied by a user-defined molecular graph.
pub trait BondView {
    /// Returns the stable identifier of this bond.
    fn id(&self) -> BondId;

    /// Reports the bond order recorded in the source graph.
    fn order(&self) -> BondOrder;

    /// Returns the identifier of the atom at the start of the bond.
    fn start_atom_id(&self) -> AtomId;

    /// Returns the identifier of the atom at the end of the bond.
    fn end_atom_id(&self) -> AtomId;
}

/// Lightweight abstraction that allows perception to operate on any graph implementation.
pub trait MoleculeGraph {
    /// Atom view type provided by the graph.
    type Atom: AtomView;
    /// Bond view type provided by the graph.
    type Bond: BondView;

    /// Returns an iterator over all atoms in the graph.
    fn atoms(&self) -> impl Iterator<Item = &Self::Atom>;

    /// Returns an iterator over all bonds in the graph.
    fn bonds(&self) -> impl Iterator<Item = &Self::Bond>;
}
