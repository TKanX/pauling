//! A cheminformatics library in pure Rust for the perception of chemical resonance
//! from molecular graphs.
//!
//! The theory of chemical resonance was first developed by Linus **Pauling** in the 1930s.
//!
//! # Overview
//!
//! `pauling` perceives resonance by executing a robust, multi-stage pipeline that
//! analyzes a molecular graph's topology. The primary entry point is the
//! [`find_resonance_systems`] function, which performs the full workflow:
//!
//! 1.  **Ring Perception:** Identifies the Smallest Set of Smallest Rings (SSSR).
//! 2.  **Aromaticity Perception:** Determines aromatic systems using Hückel's rule.
//! 3.  **Kekulization:** Assigns a valid Kekulé structure to aromatic rings.
//! 4.  **Atom State Perception:** Calculates valence, lone pairs, and hybridization.
//! 5.  **Resonance Identification:** Groups all conjugated atoms and bonds into
//!     distinct resonance systems.
//!
//! The library is designed to be flexible. It operates on any data structure that
//! implements the [`traits::MoleculeGraph`] trait, allowing seamless integration with
//! existing molecular modeling projects. For convenience, a simple [`Molecule`]
//! implementation is provided for testing and examples.
//!
//! # Quick Start
//!
//! ```
//! use pauling::{find_resonance_systems, BondOrder, Element, Molecule, PerceptionError};
//!
//! // 1. Build a benzene molecule (C6H6) with all atoms explicitly defined.
//! let mut benzene = Molecule::new();
//!
//! // First, create the 6 carbon atoms for the ring.
//! let carbons: Vec<_> = (0..6).map(|_| benzene.add_atom(Element::C, 0)).collect();
//!
//! // Then, create the alternating bonds of the carbon backbone.
//! let mut ring_bonds = Vec::new();
//! ring_bonds.push(benzene.add_bond(carbons[0], carbons[1], BondOrder::Double).unwrap());
//! ring_bonds.push(benzene.add_bond(carbons[1], carbons[2], BondOrder::Single).unwrap());
//! ring_bonds.push(benzene.add_bond(carbons[2], carbons[3], BondOrder::Double).unwrap());
//! ring_bonds.push(benzene.add_bond(carbons[3], carbons[4], BondOrder::Single).unwrap());
//! ring_bonds.push(benzene.add_bond(carbons[4], carbons[5], BondOrder::Double).unwrap());
//! ring_bonds.push(benzene.add_bond(carbons[5], carbons[0], BondOrder::Single).unwrap());
//!
//! // Finally, create and attach the 6 hydrogen atoms.
//! for &carbon_id in &carbons {
//!     let hydrogen_id = benzene.add_atom(Element::H, 0);
//!     benzene.add_bond(carbon_id, hydrogen_id, BondOrder::Single).unwrap();
//! }
//!
//! // 2. Run the perception pipeline.
//! let systems = find_resonance_systems(&benzene)?;
//!
//! // 3. Analyze the results.
//! // Benzene has one resonance system, which consists of only the carbon ring.
//! assert_eq!(systems.len(), 1);
//! let system = &systems[0];
//!
//! // Verify that the system contains exactly the 6 carbon atoms.
//! let mut system_atoms = system.atoms.clone();
//! system_atoms.sort();
//! assert_eq!(system_atoms, carbons);
//!
//! // Verify that the system contains exactly the 6 carbon-carbon bonds.
//! let mut system_bonds = system.bonds.clone();
//! system_bonds.sort();
//! ring_bonds.sort();
//! assert_eq!(system_bonds, ring_bonds);
//! # Ok::<(), PerceptionError>(())
//! ```

mod core;
mod errors;
mod graph;
mod molecule;
mod perception;
mod resonance;

/// The primary entry point to the `pauling` perception pipeline.
pub use crate::find_resonance_systems_impl::find_resonance_systems;

/// A stable, user-facing identifier for an atom.
pub use core::atom::AtomId;
/// An enumeration of chemical elements.
pub use core::atom::Element;
/// A stable, user-facing identifier for a bond.
pub use core::bond::BondId;
/// An enumeration of bond orders (Single, Double, etc.).
pub use core::bond::BondOrder;

/// The error type for all fallible perception operations.
pub use errors::PerceptionError;
/// Represents a single, connected network of conjugated atoms and bonds.
pub use resonance::ResonanceSystem;

/// A simple, in-memory molecular graph implementation for examples and testing.
pub use molecule::Molecule;
/// Errors that can occur during the construction of a [`Molecule`].
pub use molecule::MoleculeBuildError;

/// The core traits (`MoleculeGraph`, `AtomView`, `BondView`) for graph abstraction.
pub use crate::graph::traits;

mod find_resonance_systems_impl {
    use super::*;
    use crate::graph::traits::MoleculeGraph;
    use crate::perception::ChemicalPerception;

    /// Finds all distinct resonance systems present in a molecular graph.
    ///
    /// This function serves as the main entry point to the `pauling` library.
    /// It builds a `ChemicalPerception` from the provided graph and executes
    /// the full perception pipeline, including ring detection, aromaticity
    /// perception, Kekulé assignment, atomic state perception, and resonance
    /// grouping.
    ///
    /// # Arguments
    ///
    /// * `graph` - A reference to any type that implements the [`MoleculeGraph`]
    ///   trait. The graph is treated as read-only.
    ///
    /// # Returns
    ///
    /// On success, returns a `Vec<ResonanceSystem>`. Each system contains the
    /// atom and bond identifiers (as provided by the input graph) that belong
    /// to a single, continuous conjugated network.
    ///
    /// # Errors
    ///
    /// Returns a [`PerceptionError`] if the input graph is structurally
    /// inconsistent (e.g., contains dangling bonds or duplicate bonds), or if a
    /// perception stage fails (e.g., Kekulization does not converge for an
    /// aromatic system).
    ///
    /// # Examples
    ///
    /// Build the formate anion (`HCOO⁻`) with all atoms explicitly defined.
    ///
    /// ```
    /// use pauling::{find_resonance_systems, BondOrder, Element, Molecule, PerceptionError};
    ///
    /// let mut molecule = Molecule::new();
    ///
    /// // Define all atoms for one of its resonance structures: H-C(=O)O⁻
    /// let c = molecule.add_atom(Element::C, 0);
    /// let h = molecule.add_atom(Element::H, 0);
    /// let o_double = molecule.add_atom(Element::O, 0);    // The double-bonded oxygen
    /// let o_single = molecule.add_atom(Element::O, -1);   // The single-bonded, negatively charged oxygen
    ///
    /// // Define all bonds
    /// molecule.add_bond(c, h, BondOrder::Single).unwrap();
    /// let bond_co_double = molecule.add_bond(c, o_double, BondOrder::Double).unwrap();
    /// let bond_co_single = molecule.add_bond(c, o_single, BondOrder::Single).unwrap();
    ///
    /// // Run perception
    /// let systems = find_resonance_systems(&molecule)?;
    ///
    /// // There is one resonance system corresponding to the O-C-O backbone.
    /// assert_eq!(systems.len(), 1);
    /// let system = &systems[0];
    ///
    /// // The hydrogen atom is NOT part of the resonance system.
    /// let mut found_atoms = system.atoms.clone();
    /// found_atoms.sort();
    /// assert_eq!(found_atoms, vec![c, o_double, o_single]);
    ///
    /// // The C-H bond is NOT part of the resonance system.
    /// let mut found_bonds = system.bonds.clone();
    /// found_bonds.sort();
    /// assert_eq!(found_bonds, vec![bond_co_double, bond_co_single]);
    ///
    /// # Ok::<(), PerceptionError>(())
    /// ```
    pub fn find_resonance_systems<G: MoleculeGraph>(
        graph: &G,
    ) -> Result<Vec<ResonanceSystem>, PerceptionError> {
        let perception = ChemicalPerception::from_graph(graph)?;

        let systems = resonance::find_systems(&perception);

        Ok(systems)
    }
}
