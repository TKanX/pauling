//! Simple in-memory molecule for quick experiments and doctests.
//!
//! The `Molecule` type implements [`MoleculeGraph`](crate::traits::MoleculeGraph)
//! so it can be fed directly into the perception pipeline. It is deliberately
//! minimal and performs only basic validation on insertions.

use crate::core::atom::{AtomId, Element};
use crate::core::bond::{BondId, BondOrder};
use crate::graph::traits::{AtomView, BondView, MoleculeGraph};
use thiserror::Error;

/// Error emitted when an invalid atom or bond is added to a [`Molecule`].
#[derive(Debug, Error, PartialEq, Eq)]
pub enum MoleculeBuildError {
    /// Referenced atom identifier is absent from the molecule.
    #[error("atom ID {0} is out of bounds (highest ID is {1})")]
    AtomNotFound(AtomId, AtomId),

    /// Attempted to create a second bond between the same atom pair.
    #[error("duplicate bond: a bond already exists between atoms {0} and {1}")]
    DuplicateBond(AtomId, AtomId),

    /// Attempted to connect an atom to itself.
    #[error("self-loop bond is not allowed on atom {0}")]
    SelfLoopBond(AtomId),
}

/// Concrete atom data stored inside [`Molecule`].
#[derive(Clone, Debug)]
pub struct Atom {
    id: AtomId,
    element: Element,
    formal_charge: i8,
}

impl AtomView for Atom {
    fn id(&self) -> AtomId {
        self.id
    }
    fn element(&self) -> Element {
        self.element
    }
    fn formal_charge(&self) -> i8 {
        self.formal_charge
    }
}

/// Concrete bond data stored inside [`Molecule`].
#[derive(Clone, Debug)]
pub struct Bond {
    id: BondId,
    order: BondOrder,
    start: AtomId,
    end: AtomId,
}

impl BondView for Bond {
    fn id(&self) -> BondId {
        self.id
    }
    fn order(&self) -> BondOrder {
        self.order
    }
    fn start_atom_id(&self) -> AtomId {
        self.start
    }
    fn end_atom_id(&self) -> AtomId {
        self.end
    }
}

/// Lightweight adjacency-based molecule that implements [`MoleculeGraph`].
#[derive(Clone, Debug, Default)]
pub struct Molecule {
    atoms: Vec<Atom>,
    bonds: Vec<Bond>,
    adjacency: Vec<Vec<BondId>>,
}

impl Molecule {
    /// Creates an empty molecule with no atoms or bonds.
    pub fn new() -> Self {
        Self::default()
    }

    /// Inserts a new atom and returns its [`AtomId`].
    ///
    /// # Arguments
    ///
    /// * `element` - Chemical element of the atom.
    /// * `formal_charge` - Formal charge in units of `e`.
    ///
    /// # Returns
    ///
    /// The newly assigned [`AtomId`], which equals the insertion index.
    pub fn add_atom(&mut self, element: Element, formal_charge: i8) -> AtomId {
        let id = self.atoms.len();
        self.atoms.push(Atom {
            id,
            element,
            formal_charge,
        });
        self.adjacency.push(Vec::new());
        id
    }

    /// Connects two atoms with a bond of the given order.
    ///
    /// # Arguments
    ///
    /// * `start_id` - Identifier of the first atom.
    /// * `end_id` - Identifier of the second atom.
    /// * `order` - Desired [`BondOrder`].
    ///
    /// # Returns
    ///
    /// The [`BondId`] assigned to the new bond on success.
    ///
    /// # Errors
    ///
    /// * [`MoleculeBuildError::SelfLoopBond`] if `start_id == end_id`.
    /// * [`MoleculeBuildError::AtomNotFound`] if either atom is missing.
    /// * [`MoleculeBuildError::DuplicateBond`] if a bond already connects the atoms.
    pub fn add_bond(
        &mut self,
        start_id: AtomId,
        end_id: AtomId,
        order: BondOrder,
    ) -> Result<BondId, MoleculeBuildError> {
        if start_id == end_id {
            return Err(MoleculeBuildError::SelfLoopBond(start_id));
        }

        let max_id = self.atoms.len().saturating_sub(1);
        if start_id >= self.atoms.len() {
            return Err(MoleculeBuildError::AtomNotFound(start_id, max_id));
        }
        if end_id >= self.atoms.len() {
            return Err(MoleculeBuildError::AtomNotFound(end_id, max_id));
        }

        let check_atom = if self.adjacency[start_id].len() < self.adjacency[end_id].len() {
            start_id
        } else {
            end_id
        };

        for bond_id in &self.adjacency[check_atom] {
            let bond = &self.bonds[*bond_id];
            if (bond.start == start_id && bond.end == end_id)
                || (bond.start == end_id && bond.end == start_id)
            {
                return Err(MoleculeBuildError::DuplicateBond(start_id, end_id));
            }
        }

        let id = self.bonds.len();
        self.bonds.push(Bond {
            id,
            order,
            start: start_id,
            end: end_id,
        });

        self.adjacency[start_id].push(id);
        self.adjacency[end_id].push(id);

        Ok(id)
    }

    /// Returns an immutable view of an atom when the identifier is valid.
    ///
    /// # Arguments
    ///
    /// * `id` - Atom identifier obtained from this molecule.
    ///
    /// # Returns
    ///
    /// `Some(&Atom)` when the identifier exists, otherwise `None`.
    pub fn atom(&self, id: AtomId) -> Option<&Atom> {
        self.atoms.get(id)
    }

    /// Returns an immutable view of a bond when the identifier is valid.
    ///
    /// # Arguments
    ///
    /// * `id` - Bond identifier obtained from this molecule.
    ///
    /// # Returns
    ///
    /// `Some(&Bond)` when the identifier exists, otherwise `None`.
    pub fn bond(&self, id: BondId) -> Option<&Bond> {
        self.bonds.get(id)
    }

    /// Iterates over all bonds incident to the specified atom.
    ///
    /// # Arguments
    ///
    /// * `id` - Atom identifier obtained from this molecule.
    ///
    /// # Returns
    ///
    /// An iterator producing the [`BondId`] values of all adjacent bonds.
    pub fn bonds_of_atom(&self, id: AtomId) -> impl Iterator<Item = BondId> + '_ {
        self.adjacency.get(id).into_iter().flatten().copied()
    }
}

impl MoleculeGraph for Molecule {
    type Atom = Atom;
    type Bond = Bond;

    fn atoms(&self) -> impl Iterator<Item = &Self::Atom> {
        self.atoms.iter()
    }

    fn bonds(&self) -> impl Iterator<Item = &Self::Bond> {
        self.bonds.iter()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::atom::Element;
    use crate::core::bond::BondOrder;

    #[test]
    fn add_atom_assigns_incrementing_ids_and_stores_properties() {
        let mut molecule = Molecule::new();

        let carbon_id = molecule.add_atom(Element::C, 0);
        let oxygen_id = molecule.add_atom(Element::O, -1);

        assert_eq!(carbon_id, 0);
        assert_eq!(oxygen_id, 1);

        let carbon = molecule.atom(carbon_id).expect("carbon atom missing");
        assert_eq!(carbon.id(), carbon_id);
        assert_eq!(carbon.element(), Element::C);
        assert_eq!(carbon.formal_charge(), 0);

        assert!(molecule.atom(2).is_none());
    }

    #[test]
    fn add_bond_connects_atoms_and_updates_adjacency() {
        let mut molecule = Molecule::new();
        let carbon_id = molecule.add_atom(Element::C, 0);
        let oxygen_id = molecule.add_atom(Element::O, 0);

        let bond_id = molecule
            .add_bond(carbon_id, oxygen_id, BondOrder::Double)
            .expect("bond creation failed");

        assert_eq!(bond_id, 0);

        let bond = molecule.bond(bond_id).expect("bond missing");
        assert_eq!(bond.id(), bond_id);
        assert_eq!(bond.start_atom_id(), carbon_id);
        assert_eq!(bond.end_atom_id(), oxygen_id);
        assert_eq!(bond.order(), BondOrder::Double);

        let bonds_of_carbon: Vec<_> = molecule.bonds_of_atom(carbon_id).collect();
        assert_eq!(bonds_of_carbon, vec![bond_id]);

        let bonds_of_oxygen: Vec<_> = molecule.bonds_of_atom(oxygen_id).collect();
        assert_eq!(bonds_of_oxygen, vec![bond_id]);
    }

    #[test]
    fn add_bond_returns_error_for_missing_atoms() {
        let mut molecule = Molecule::new();
        let carbon_id = molecule.add_atom(Element::C, 0);

        let err = molecule
            .add_bond(carbon_id, carbon_id + 1, BondOrder::Single)
            .expect_err("expected missing atom error");

        assert_eq!(
            err,
            MoleculeBuildError::AtomNotFound(carbon_id + 1, carbon_id)
        );
    }

    #[test]
    fn add_bond_returns_error_for_duplicate_edges() {
        let mut molecule = Molecule::new();
        let carbon_id = molecule.add_atom(Element::C, 0);
        let oxygen_id = molecule.add_atom(Element::O, 0);

        molecule
            .add_bond(carbon_id, oxygen_id, BondOrder::Single)
            .expect("first bond creation failed");

        let err = molecule
            .add_bond(carbon_id, oxygen_id, BondOrder::Double)
            .expect_err("expected duplicate bond error");

        assert_eq!(err, MoleculeBuildError::DuplicateBond(carbon_id, oxygen_id));

        let err_reversed = molecule
            .add_bond(oxygen_id, carbon_id, BondOrder::Single)
            .expect_err("expected duplicate bond error");

        assert_eq!(
            err_reversed,
            MoleculeBuildError::DuplicateBond(oxygen_id, carbon_id)
        );
    }

    #[test]
    fn bonds_of_atom_collects_all_incident_bonds() {
        let mut molecule = Molecule::new();
        let carbon_id = molecule.add_atom(Element::C, 0);
        let oxygen_id = molecule.add_atom(Element::O, 0);
        let hydrogen_id = molecule.add_atom(Element::H, 0);

        let first_bond = molecule
            .add_bond(carbon_id, oxygen_id, BondOrder::Double)
            .expect("first bond creation failed");
        let second_bond = molecule
            .add_bond(carbon_id, hydrogen_id, BondOrder::Single)
            .expect("second bond creation failed");

        let mut bonds_of_carbon: Vec<_> = molecule.bonds_of_atom(carbon_id).collect();
        bonds_of_carbon.sort_unstable();
        assert_eq!(bonds_of_carbon, vec![first_bond, second_bond]);

        let bonds_of_oxygen: Vec<_> = molecule.bonds_of_atom(oxygen_id).collect();
        assert_eq!(bonds_of_oxygen, vec![first_bond]);

        let bonds_of_hydrogen: Vec<_> = molecule.bonds_of_atom(hydrogen_id).collect();
        assert_eq!(bonds_of_hydrogen, vec![second_bond]);
    }
}
