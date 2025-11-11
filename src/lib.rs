mod core;
mod errors;
mod graph;
mod molecule;
mod perception;
mod resonance;

pub use crate::find_resonance_systems_impl::find_resonance_systems;

pub use core::atom::{AtomId, Element};
pub use core::bond::{BondId, BondOrder};

pub use errors::PerceptionError;
pub use resonance::ResonanceSystem;

pub use molecule::{Molecule, MoleculeBuildError};

pub use crate::graph::traits;

mod find_resonance_systems_impl {
    use super::*;
    use crate::graph::traits::MoleculeGraph;
    use crate::perception::ChemicalPerception;

    pub fn find_resonance_systems<G: MoleculeGraph>(
        graph: &G,
    ) -> Result<Vec<ResonanceSystem>, PerceptionError> {
        let perception = ChemicalPerception::from_graph(graph)?;

        let systems = resonance::find_systems(&perception);

        Ok(systems)
    }
}
