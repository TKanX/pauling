use super::system::ResonanceSystem;
use crate::core::bond::BondOrder;
use crate::perception::ChemicalPerception;
use std::collections::{HashSet, VecDeque};

pub fn find_systems(perception: &ChemicalPerception) -> Vec<ResonanceSystem> {
    if perception.bonds.is_empty() {
        return Vec::new();
    }

    let conjugated_bond_indices = find_and_expand_conjugated_bonds(perception);
    group_systems(perception, &conjugated_bond_indices)
}

fn find_and_expand_conjugated_bonds(perception: &ChemicalPerception) -> HashSet<usize> {
    let mut conjugated: HashSet<usize> = HashSet::new();
    let mut frontier: VecDeque<usize> = VecDeque::new();

    for (bond_idx, bond) in perception.bonds.iter().enumerate() {
        let effective_order = bond.kekule_order.unwrap_or(bond.order);
        if (bond.is_aromatic || matches!(effective_order, BondOrder::Double | BondOrder::Triple))
            && conjugated.insert(bond_idx)
        {
            frontier.push_back(bond_idx);
        }
    }

    while let Some(bond_idx) = frontier.pop_front() {
        let bond = &perception.bonds[bond_idx];

        for atom_id in [bond.start_atom_id, bond.end_atom_id] {
            let atom_idx = perception.atom_id_to_index[&atom_id];
            let atom = &perception.atoms[atom_idx];

            if !atom.is_conjugation_candidate {
                continue;
            }

            for &(_, neighbor_bond_id) in &perception.adjacency[atom_idx] {
                let neighbor_bond_idx = perception.bond_id_to_index[&neighbor_bond_id];

                if !conjugated.contains(&neighbor_bond_idx) {
                    let neighbor_bond = &perception.bonds[neighbor_bond_idx];
                    let other_end_id = neighbor_bond.other_end(atom.id);
                    let other_end_idx = perception.atom_id_to_index[&other_end_id];

                    if perception.atoms[other_end_idx].is_conjugation_candidate
                        && conjugated.insert(neighbor_bond_idx)
                    {
                        frontier.push_back(neighbor_bond_idx);
                    }
                }
            }
        }
    }

    conjugated
}

fn group_systems(
    perception: &ChemicalPerception,
    conjugated_bond_indices: &HashSet<usize>,
) -> Vec<ResonanceSystem> {
    let mut systems = Vec::new();
    let mut visited_bonds = vec![false; perception.bonds.len()];

    for &start_bond_idx in conjugated_bond_indices {
        if visited_bonds[start_bond_idx] {
            continue;
        }

        let mut queue = VecDeque::new();
        let mut system_bond_ids = Vec::new();
        let mut system_atom_ids = HashSet::new();

        queue.push_back(start_bond_idx);
        visited_bonds[start_bond_idx] = true;

        while let Some(bond_idx) = queue.pop_front() {
            let bond = &perception.bonds[bond_idx];
            system_bond_ids.push(bond.id);
            system_atom_ids.insert(bond.start_atom_id);
            system_atom_ids.insert(bond.end_atom_id);

            for atom_id in [bond.start_atom_id, bond.end_atom_id] {
                let atom_idx = perception.atom_id_to_index[&atom_id];
                for &(_, neighbor_bond_id) in &perception.adjacency[atom_idx] {
                    let neighbor_bond_idx = perception.bond_id_to_index[&neighbor_bond_id];
                    if conjugated_bond_indices.contains(&neighbor_bond_idx)
                        && !visited_bonds[neighbor_bond_idx]
                    {
                        visited_bonds[neighbor_bond_idx] = true;
                        queue.push_back(neighbor_bond_idx);
                    }
                }
            }
        }

        if !system_bond_ids.is_empty() {
            let system =
                ResonanceSystem::new(system_atom_ids.into_iter().collect(), system_bond_ids);
            systems.push(system);
        }
    }

    systems.sort_by(|a, b| a.bonds.cmp(&b.bonds));
    systems
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::atom::{AtomId, Element};
    use crate::core::bond::{BondId, BondOrder};
    use crate::perception::{ChemicalPerception, Hybridization, PerceivedAtom, PerceivedBond};
    use std::collections::HashMap;

    #[derive(Clone, Copy)]
    struct AtomSetup {
        element: Element,
        is_candidate: bool,
    }

    #[derive(Clone, Copy)]
    struct BondSetup {
        id: BondId,
        start: AtomId,
        end: AtomId,
        order: BondOrder,
        is_aromatic: bool,
        kekule_order: Option<BondOrder>,
    }

    impl AtomSetup {
        fn candidate(element: Element) -> Self {
            Self {
                element,
                is_candidate: true,
            }
        }

        fn saturated(element: Element) -> Self {
            Self {
                element,
                is_candidate: false,
            }
        }
    }

    impl BondSetup {
        fn new(id: BondId, start: AtomId, end: AtomId, order: BondOrder) -> Self {
            Self {
                id,
                start,
                end,
                order,
                is_aromatic: false,
                kekule_order: None,
            }
        }

        fn aromatic(mut self) -> Self {
            self.is_aromatic = true;
            self
        }

        fn with_kekule(mut self, kekule: BondOrder) -> Self {
            self.kekule_order = Some(kekule);
            self
        }
    }

    fn build_perception(atoms: &[AtomSetup], bonds: &[BondSetup]) -> ChemicalPerception {
        let mut adjacency: Vec<Vec<(usize, BondId)>> = vec![Vec::new(); atoms.len()];
        let mut bond_vec = Vec::with_capacity(bonds.len());
        let mut bond_id_to_index = HashMap::new();

        for (idx, bond) in bonds.iter().enumerate() {
            let start = bond.start as usize;
            let end = bond.end as usize;
            adjacency[start].push((end, bond.id));
            adjacency[end].push((start, bond.id));
            bond_id_to_index.insert(bond.id, idx);
            bond_vec.push(PerceivedBond {
                id: bond.id,
                order: bond.order,
                start_atom_id: bond.start,
                end_atom_id: bond.end,
                is_in_ring: false,
                is_aromatic: bond.is_aromatic,
                kekule_order: bond.kekule_order,
            });
        }

        let mut atom_vec = Vec::with_capacity(atoms.len());
        let mut atom_id_to_index = HashMap::new();
        for (idx, atom) in atoms.iter().enumerate() {
            let hybridization = if atom.is_candidate {
                Hybridization::SP2
            } else {
                Hybridization::SP3
            };

            atom_vec.push(PerceivedAtom {
                id: idx,
                element: atom.element,
                formal_charge: 0,
                total_degree: adjacency[idx].len() as u8,
                total_valence: 0,
                is_in_ring: false,
                is_aromatic: false,
                hybridization,
                is_conjugation_candidate: atom.is_candidate,
                lone_pairs: 0,
            });
            atom_id_to_index.insert(idx, idx);
        }

        ChemicalPerception {
            atoms: atom_vec,
            bonds: bond_vec,
            adjacency,
            atom_id_to_index,
            bond_id_to_index,
            ring_info: Default::default(),
        }
    }

    #[test]
    fn empty_perception_returns_no_systems() {
        let perception = ChemicalPerception {
            atoms: Vec::new(),
            bonds: Vec::new(),
            adjacency: Vec::new(),
            atom_id_to_index: HashMap::new(),
            bond_id_to_index: HashMap::new(),
            ring_info: Default::default(),
        };

        assert!(find_systems(&perception).is_empty());
    }

    #[test]
    fn isolated_double_bond_forms_single_resonance_system() {
        let perception = build_perception(
            &[
                AtomSetup::candidate(Element::C),
                AtomSetup::candidate(Element::C),
            ],
            &[BondSetup::new(0, 0, 1, BondOrder::Double)],
        );

        let systems = find_systems(&perception);
        assert_eq!(systems.len(), 1);
        assert_eq!(systems[0].atoms, vec![0, 1]);
        assert_eq!(systems[0].bonds, vec![0]);
    }

    #[test]
    fn conjugation_expands_through_candidate_single_bonds() {
        let perception = build_perception(
            &[
                AtomSetup::candidate(Element::C),
                AtomSetup::candidate(Element::C),
                AtomSetup::candidate(Element::C),
            ],
            &[
                BondSetup::new(0, 0, 1, BondOrder::Double),
                BondSetup::new(1, 1, 2, BondOrder::Single),
            ],
        );

        let systems = find_systems(&perception);
        assert_eq!(systems.len(), 1);
        assert_eq!(systems[0].atoms, vec![0, 1, 2]);
        assert_eq!(systems[0].bonds, vec![0, 1]);
    }

    #[test]
    fn non_candidate_neighbors_block_conjugation_growth() {
        let perception = build_perception(
            &[
                AtomSetup::candidate(Element::C),
                AtomSetup::candidate(Element::C),
                AtomSetup::saturated(Element::C),
            ],
            &[
                BondSetup::new(0, 0, 1, BondOrder::Double),
                BondSetup::new(1, 1, 2, BondOrder::Single),
            ],
        );

        let systems = find_systems(&perception);
        assert_eq!(systems.len(), 1);
        assert_eq!(systems[0].atoms, vec![0, 1]);
        assert_eq!(systems[0].bonds, vec![0]);
    }

    #[test]
    fn disconnected_sets_remain_distinct_systems() {
        let perception = build_perception(
            &[
                AtomSetup::candidate(Element::C),
                AtomSetup::candidate(Element::C),
                AtomSetup::candidate(Element::C),
                AtomSetup::candidate(Element::C),
            ],
            &[
                BondSetup::new(2, 0, 1, BondOrder::Double),
                BondSetup::new(5, 2, 3, BondOrder::Double),
            ],
        );

        let systems = find_systems(&perception);
        assert_eq!(systems.len(), 2);
        assert_eq!(systems[0].atoms, vec![0, 1]);
        assert_eq!(systems[0].bonds, vec![2]);
        assert_eq!(systems[1].atoms, vec![2, 3]);
        assert_eq!(systems[1].bonds, vec![5]);
    }

    #[test]
    fn aromatic_and_kekule_metadata_seed_conjugation() {
        let perception = build_perception(
            &[
                AtomSetup::candidate(Element::C),
                AtomSetup::candidate(Element::C),
                AtomSetup::candidate(Element::C),
            ],
            &[
                BondSetup::new(0, 0, 1, BondOrder::Single).aromatic(),
                BondSetup::new(3, 1, 2, BondOrder::Single).with_kekule(BondOrder::Double),
            ],
        );

        let systems = find_systems(&perception);
        assert_eq!(systems.len(), 1);
        assert_eq!(systems[0].atoms, vec![0, 1, 2]);
        assert_eq!(systems[0].bonds, vec![0, 3]);
    }
}
