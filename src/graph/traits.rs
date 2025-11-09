use crate::core::atom::{AtomId, Element};
use crate::core::bond::{BondId, BondOrder};

pub trait AtomView {
    fn id(&self) -> AtomId;

    fn element(&self) -> Element;

    fn formal_charge(&self) -> i8;
}

pub trait BondView {
    fn id(&self) -> BondId;

    fn order(&self) -> BondOrder;

    fn start_atom_id(&self) -> AtomId;

    fn end_atom_id(&self) -> AtomId;
}

pub trait MoleculeGraph {
    type Atom: AtomView;
    type Bond: BondView;

    fn atoms(&self) -> impl Iterator<Item = &Self::Atom>;

    fn bonds(&self) -> impl Iterator<Item = &Self::Bond>;
}
