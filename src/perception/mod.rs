//! Internal perception pipeline that enriches molecular graphs with chemical metadata.
//!
//! The module orchestrates ring perception, aromaticity detection, Kekulé
//! assignments, atomic state inference, and resonance candidate identification.

use crate::core::atom::{AtomId, Element};
use crate::core::bond::{BondId, BondOrder};
use crate::errors::PerceptionError;
use crate::graph::traits::{AtomView, BondView, MoleculeGraph};
use crate::perception::ring::RingInfo;
use crate::resonance;
use std::collections::{HashMap, HashSet};

mod aromaticity;
mod kekulize;
mod ring;
mod state;

/// Hybridization states assigned to perceived atoms.
pub use state::Hybridization;

/// Atom annotated with metadata derived from the perception pipeline.
#[derive(Clone, Debug)]
pub struct PerceivedAtom {
    /// Stable identifier matching the input graph.
    pub id: AtomId,
    /// Chemical element provided by the input graph.
    pub element: Element,
    /// Formal charge supplied by the source molecule.
    pub formal_charge: i8,
    /// Number of adjacent bonds in the original graph.
    pub total_degree: u8,
    /// Sum of bond multiplicities including Kekulé adjustments.
    pub total_valence: u8,
    /// Indicates whether the atom belongs to at least one ring in the SSSR set.
    pub is_in_ring: bool,
    /// Indicates whether the atom participates in an aromatic system.
    pub is_aromatic: bool,
    /// Hybridization state inferred during the perception pipeline.
    pub hybridization: Hybridization,
    /// Flag denoting participation eligibility in conjugation/resonance searches.
    pub is_conjugation_candidate: bool,
    /// Estimated number of lone pairs according to valence heuristics.
    pub lone_pairs: u8,
}

impl PerceivedAtom {
    /// Creates a perceived atom with default perception metadata.
    fn new(id: AtomId, element: Element, formal_charge: i8, total_degree: u8) -> Self {
        Self {
            id,
            element,
            formal_charge,
            total_degree,
            total_valence: 0,
            is_in_ring: false,
            is_aromatic: false,
            hybridization: Hybridization::Unknown,
            is_conjugation_candidate: false,
            lone_pairs: 0,
        }
    }
}

/// Bond annotated with metadata derived from the perception pipeline.
#[derive(Clone, Debug)]
pub struct PerceivedBond {
    /// Stable identifier matching the input graph.
    pub id: BondId,
    /// Bond order supplied by the input graph.
    pub order: BondOrder,
    /// Identifier of one endpoint in the input graph.
    pub start_atom_id: AtomId,
    /// Identifier of the other endpoint in the input graph.
    pub end_atom_id: AtomId,
    /// Indicates whether the bond belongs to at least one ring in the SSSR set.
    pub is_in_ring: bool,
    /// Indicates whether the bond participates in an aromatic system.
    pub is_aromatic: bool,
    /// Kekulé order assigned during Kekulization, when applicable.
    pub kekule_order: Option<BondOrder>,
}

impl PerceivedBond {
    /// Creates a perceived bond with default perception metadata.
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

    /// Returns the opposite atom identifier given one endpoint.
    ///
    /// # Arguments
    ///
    /// * `atom_id` - The known endpoint of the bond.
    ///
    /// # Returns
    ///
    /// The identifier of the other atom connected by the bond. The function
    /// assumes `atom_id` matches either `start_atom_id` or `end_atom_id`.
    pub fn other_end(&self, atom_id: AtomId) -> AtomId {
        if self.start_atom_id == atom_id {
            self.end_atom_id
        } else {
            self.start_atom_id
        }
    }
}

/// Comprehensive snapshot of all perception-derived metadata.
pub struct ChemicalPerception {
    /// Atom-centric perception data.
    pub atoms: Vec<PerceivedAtom>,
    /// Bond-centric perception data.
    pub bonds: Vec<PerceivedBond>,
    /// Adjacency list referencing indices and bond identifiers.
    pub adjacency: Vec<Vec<(usize, BondId)>>,

    /// Maps external `AtomId` values to internal vector indices.
    pub atom_id_to_index: HashMap<AtomId, usize>,
    /// Maps external `BondId` values to internal vector indices.
    pub bond_id_to_index: HashMap<BondId, usize>,

    /// Ring data detected during the perception pipeline.
    pub ring_info: RingInfo,
}

impl ChemicalPerception {
    /// Builds a `ChemicalPerception` from any [`MoleculeGraph`].
    ///
    /// The function copies core topology data, enriches it with ring and
    /// aromaticity annotations, assigns Kekulé resonance orders, infers atomic
    /// states, and finally flags resonance candidates.
    ///
    /// # Arguments
    ///
    /// * `graph` - An implementation of [`MoleculeGraph`].
    ///
    /// # Returns
    ///
    /// A fully populated `ChemicalPerception` ready for downstream resonance
    /// identification.
    ///
    /// # Errors
    ///
    /// Propagates [`PerceptionError`] variants when the input graph contains
    /// structural inconsistencies or when intermediate perception stages fail.
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

        aromaticity::perceive(&mut perception);

        kekulize::kekulize(&mut perception)?;

        state::perceive(&mut perception);

        resonance::candidate::determine(&mut perception);

        Ok(perception)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::atom::Element;
    use crate::core::bond::BondId;
    use crate::core::bond::BondOrder;
    use crate::molecule::Molecule;

    fn attach_hydrogen(molecule: &mut Molecule, atom: AtomId) {
        let h = molecule.add_atom(Element::H, 0);
        molecule
            .add_bond(atom, h, BondOrder::Single)
            .expect("attach hydrogen");
    }

    fn build_benzene() -> (Molecule, Vec<BondId>, Vec<AtomId>) {
        let mut molecule = Molecule::new();
        let atoms: Vec<AtomId> = (0..6).map(|_| molecule.add_atom(Element::C, 0)).collect();
        let mut ring_bonds = Vec::new();

        let orders = [
            BondOrder::Double,
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
        ];

        for i in 0..6 {
            let next = (i + 1) % 6;
            let bond_id = molecule
                .add_bond(atoms[i], atoms[next], orders[i])
                .expect("add aromatic bond");
            ring_bonds.push(bond_id);
        }

        for &carbon in &atoms {
            attach_hydrogen(&mut molecule, carbon);
        }

        (molecule, ring_bonds, atoms)
    }

    fn build_biphenyl() -> Molecule {
        let mut molecule = Molecule::new();
        let left: Vec<AtomId> = (0..6).map(|_| molecule.add_atom(Element::C, 0)).collect();
        let right: Vec<AtomId> = (0..6).map(|_| molecule.add_atom(Element::C, 0)).collect();

        let orders = [
            BondOrder::Double,
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
        ];

        for i in 0..6 {
            let next = (i + 1) % 6;
            molecule
                .add_bond(left[i], left[next], orders[i])
                .expect("left ring bond");
            molecule
                .add_bond(right[i], right[next], orders[i])
                .expect("right ring bond");
        }

        molecule
            .add_bond(left[1], right[4], BondOrder::Single)
            .expect("biaryl link");

        for (i, &carbon) in left.iter().enumerate() {
            if i == 1 {
                continue;
            }
            attach_hydrogen(&mut molecule, carbon);
        }

        for (i, &carbon) in right.iter().enumerate() {
            if i == 4 {
                continue;
            }
            attach_hydrogen(&mut molecule, carbon);
        }

        molecule
    }

    fn build_acetamide() -> Molecule {
        let mut molecule = Molecule::new();
        let carbonyl_c = molecule.add_atom(Element::C, 0);
        let oxygen = molecule.add_atom(Element::O, 0);
        let nitrogen = molecule.add_atom(Element::N, 0);
        let methyl_carbon = molecule.add_atom(Element::C, 0);

        molecule
            .add_bond(carbonyl_c, oxygen, BondOrder::Double)
            .expect("C=O");
        molecule
            .add_bond(carbonyl_c, nitrogen, BondOrder::Single)
            .expect("C-N");
        molecule
            .add_bond(carbonyl_c, methyl_carbon, BondOrder::Single)
            .expect("C-C");

        for _ in 0..3 {
            attach_hydrogen(&mut molecule, methyl_carbon);
        }

        for _ in 0..2 {
            attach_hydrogen(&mut molecule, nitrogen);
        }

        molecule
    }

    fn build_linear_butane() -> Molecule {
        let mut molecule = Molecule::new();
        let carbons: Vec<AtomId> = (0..4).map(|_| molecule.add_atom(Element::C, 0)).collect();

        for i in 0..3 {
            molecule
                .add_bond(carbons[i], carbons[i + 1], BondOrder::Single)
                .expect("C-C");
        }

        for (i, &carbon) in carbons.iter().enumerate() {
            let hydrogens = if i == 0 || i == 3 { 3 } else { 2 };
            for _ in 0..hydrogens {
                attach_hydrogen(&mut molecule, carbon);
            }
        }

        molecule
    }

    #[test]
    fn benzene_pipeline_sets_aromatic_and_kekule_metadata_consistently() {
        let (molecule, ring_bonds, atoms) = build_benzene();
        let perception = ChemicalPerception::from_graph(&molecule).expect("perception");

        assert_eq!(perception.ring_info.rings.len(), 1, "expected single ring");
        let ring = &perception.ring_info.rings[0];
        assert_eq!(ring.atom_ids.len(), 6);
        assert_eq!(ring.bond_ids.len(), 6);

        for &carbon in &atoms {
            let idx = perception.atom_id_to_index[&carbon];
            let atom = &perception.atoms[idx];
            assert!(atom.is_aromatic, "carbon should be aromatic");
            assert_eq!(atom.hybridization, Hybridization::SP2);
            assert!(atom.is_conjugation_candidate, "carbon should conjugate");
            assert_eq!(atom.total_valence, 4, "carbon valence");
        }

        let mut double_count = 0;
        for &bond_id in &ring_bonds {
            let idx = perception.bond_id_to_index[&bond_id];
            let bond = &perception.bonds[idx];
            assert!(bond.is_aromatic, "bond should be aromatic");
            let kekule = bond
                .kekule_order
                .expect("aromatic bond must have kekule order");
            if kekule == BondOrder::Double {
                double_count += 1;
            }
        }
        assert_eq!(
            double_count, 3,
            "expected three double bonds in kekule pattern"
        );
    }

    #[test]
    fn acetamide_pipeline_combines_state_and_resonance_inference() {
        let molecule = build_acetamide();
        let perception = ChemicalPerception::from_graph(&molecule).expect("perception");

        assert!(perception.ring_info.rings.is_empty(), "no rings expected");

        let mut carbonyl_c = None;
        let mut oxygen = None;
        let mut nitrogen = None;
        for atom in &perception.atoms {
            match atom.element {
                Element::C if atom.total_degree == 3 => carbonyl_c = Some(atom.clone()),
                Element::O => oxygen = Some(atom.clone()),
                Element::N => nitrogen = Some(atom.clone()),
                _ => {}
            }
        }

        let carbonyl_c = carbonyl_c.expect("missing carbonyl carbon");
        assert_eq!(carbonyl_c.hybridization, Hybridization::SP2);
        assert_eq!(carbonyl_c.total_valence, 4);
        assert!(carbonyl_c.is_conjugation_candidate);

        let oxygen = oxygen.expect("missing oxygen");
        assert_eq!(oxygen.lone_pairs, 2);
        assert_eq!(oxygen.total_valence, 2);
        assert_eq!(oxygen.hybridization, Hybridization::SP2);

        let nitrogen = nitrogen.expect("missing nitrogen");
        assert_eq!(nitrogen.hybridization, Hybridization::SP2);
        assert_eq!(nitrogen.lone_pairs, 1);
        assert!(
            nitrogen.is_conjugation_candidate,
            "amide nitrogen should conjugate"
        );

        for bond in &perception.bonds {
            if bond.order == BondOrder::Double {
                assert!(
                    bond.kekule_order.is_none(),
                    "non-aromatic carbonyl bond should not have kekule order"
                );
            }
        }
    }

    #[test]
    fn biphenyl_pipeline_identifies_two_independent_ring_components() {
        let molecule = build_biphenyl();
        let perception = ChemicalPerception::from_graph(&molecule).expect("perception");

        assert_eq!(perception.ring_info.rings.len(), 2, "expected two rings");

        let aromatic_atoms: Vec<_> = perception
            .atoms
            .iter()
            .filter(|atom| atom.is_aromatic)
            .collect();
        assert_eq!(
            aromatic_atoms.len(),
            12,
            "two phenyl rings should be aromatic"
        );

        for atom in aromatic_atoms {
            assert_eq!(atom.hybridization, Hybridization::SP2);
            assert!(atom.is_conjugation_candidate);
        }

        for bond in &perception.bonds {
            if bond.order == BondOrder::Single
                && bond.start_atom_id != bond.end_atom_id
                && perception.atoms[perception.atom_id_to_index[&bond.start_atom_id]].is_aromatic
                && perception.atoms[perception.atom_id_to_index[&bond.end_atom_id]].is_aromatic
            {
                if let Some(kekule) = bond.kekule_order {
                    assert_ne!(kekule, BondOrder::Double);
                }
            }
        }
    }

    #[test]
    fn saturated_chain_pipeline_yields_no_resonance_features() {
        let molecule = build_linear_butane();
        let perception = ChemicalPerception::from_graph(&molecule).expect("perception");

        assert!(perception.ring_info.rings.is_empty());

        for atom in &perception.atoms {
            assert!(!atom.is_aromatic);
            assert!(!atom.is_conjugation_candidate);
            match atom.element {
                Element::C => assert_eq!(atom.hybridization, Hybridization::SP3),
                Element::H => assert_eq!(atom.hybridization, Hybridization::Unknown),
                _ => {}
            }
            assert_eq!(atom.lone_pairs, 0);
        }

        for bond in &perception.bonds {
            assert!(!bond.is_aromatic);
            assert!(bond.kekule_order.is_none());
        }
    }
}
