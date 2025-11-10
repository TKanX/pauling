use crate::core::atom::{AtomId, Element};
use crate::core::bond::{BondId, BondOrder};
use crate::errors::PerceptionError;
use crate::graph::traits::{AtomView, BondView, MoleculeGraph};
use crate::perception::ring::RingInfo;
use std::collections::{HashMap, HashSet};

mod ring;

// Placeholder for future perception modules.
// mod aromaticity;
// mod kekulize;
// mod state;

#[derive(Clone, Debug)]
pub struct PerceivedAtom {
    pub id: AtomId,
    pub element: Element,
    pub formal_charge: i8,
    pub total_degree: u8,
    pub total_valence: u8,
    pub is_in_ring: bool,
    pub is_aromatic: bool,
    // --- Fields for Resonance Analysis ---
    // pub hybridization: Hybridization,
    // pub is_conjugation_candidate: bool,
    // pub lone_pairs: u8,
}

impl PerceivedAtom {
    fn new(id: AtomId, element: Element, formal_charge: i8, total_degree: u8) -> Self {
        Self {
            id,
            element,
            formal_charge,
            total_degree,
            total_valence: 0,
            is_in_ring: false,
            is_aromatic: false,
            // hybridization: Hybridization::Unknown,
            // is_conjugation_candidate: false,
            // lone_pairs: 0,
        }
    }
}

#[derive(Clone, Debug)]
pub struct PerceivedBond {
    pub id: BondId,
    pub order: BondOrder,
    pub start_atom_id: AtomId,
    pub end_atom_id: AtomId,
    pub is_in_ring: bool,
    pub is_aromatic: bool,
    pub kekule_order: Option<BondOrder>,
}

impl PerceivedBond {
    fn new(id: BondId, order: BondOrder, start_atom_id: AtomId, end_atom_id: AtomId) -> Self {
        Self {
            id,
            order,
            start_atom_id,
            end_atom_id,
            is_in_ring: false,
            is_aromatic: false,
            kekule_order: None,
        }
    }

    pub fn other_end(&self, atom_id: AtomId) -> AtomId {
        if self.start_atom_id == atom_id {
            self.end_atom_id
        } else {
            self.start_atom_id
        }
    }
}

pub struct ChemicalPerception {
    pub atoms: Vec<PerceivedAtom>,
    pub bonds: Vec<PerceivedBond>,
    pub adjacency: Vec<Vec<(usize, BondId)>>,

    pub atom_id_to_index: HashMap<AtomId, usize>,
    pub bond_id_to_index: HashMap<BondId, usize>,

    pub ring_info: RingInfo,
}

impl ChemicalPerception {
    pub fn from_graph<G>(graph: &G) -> Result<Self, PerceptionError>
    where
        G: MoleculeGraph,
    {
        let mut atom_id_to_index = HashMap::new();
        let source_atoms: Vec<_> = graph.atoms().collect();
        for (idx, atom_view) in source_atoms.iter().enumerate() {
            atom_id_to_index.insert(atom_view.id(), idx);
        }

        let num_atoms = source_atoms.len();
        let mut adjacency: Vec<Vec<(usize, BondId)>> = vec![Vec::new(); num_atoms];
        let mut bonds = Vec::new();
        let mut bond_id_to_index = HashMap::new();
        let mut seen_bonds = HashSet::new();

        for bond_view in graph.bonds() {
            let start_id = bond_view.start_atom_id();
            let end_id = bond_view.end_atom_id();

            let canonical_start = start_id.min(end_id);
            let canonical_end = start_id.max(end_id);

            if !seen_bonds.insert((canonical_start, canonical_end)) {
                return Err(PerceptionError::DuplicateBond {
                    start: start_id,
                    end: end_id,
                });
            }

            let start_idx = *atom_id_to_index
                .get(&start_id)
                .ok_or(PerceptionError::InconsistentGraph(start_id))?;
            let end_idx = *atom_id_to_index
                .get(&end_id)
                .ok_or(PerceptionError::InconsistentGraph(end_id))?;

            adjacency[start_idx].push((end_idx, bond_view.id()));
            adjacency[end_idx].push((start_idx, bond_view.id()));

            bond_id_to_index.insert(bond_view.id(), bonds.len());
            bonds.push(PerceivedBond::new(
                bond_view.id(),
                bond_view.order(),
                start_id,
                end_id,
            ));
        }

        let mut perceived_atoms = Vec::with_capacity(num_atoms);
        for atom_view in source_atoms {
            let idx = atom_id_to_index[&atom_view.id()];
            perceived_atoms.push(PerceivedAtom::new(
                atom_view.id(),
                atom_view.element(),
                atom_view.formal_charge(),
                adjacency[idx].len() as u8,
            ));
        }

        let mut perception = Self {
            atoms: perceived_atoms,
            bonds,
            adjacency,
            atom_id_to_index,
            bond_id_to_index,
            ring_info: RingInfo::default(),
        };

        let ring_info = ring::find_sssr(&perception);

        for ring in &ring_info.rings {
            for &atom_id in &ring.atom_ids {
                if let Some(&idx) = perception.atom_id_to_index.get(&atom_id) {
                    perception.atoms[idx].is_in_ring = true;
                }
            }
            for &bond_id in &ring.bond_ids {
                if let Some(&idx) = perception.bond_id_to_index.get(&bond_id) {
                    perception.bonds[idx].is_in_ring = true;
                }
            }
        }
        perception.ring_info = ring_info;

        // 2. Perceive aromaticity (placeholder).
        // aromaticity::perceive(&mut perception);

        // 3. Kekulize aromatic systems (placeholder).
        // kekulize::kekulize(&mut perception)?;

        // 4. Perceive atomic states (valence, hybridization) (placeholder).
        // state::perceive(&mut perception);

        Ok(perception)
    }
}
