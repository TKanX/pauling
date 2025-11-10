use crate::core::atom::AtomId;
use crate::core::bond::BondOrder;
use crate::errors::PerceptionError;
use crate::perception::ChemicalPerception;
use std::collections::{HashMap, VecDeque};

const KEKULIZATION_ATTEMPT_LIMIT: usize = 1000;

pub fn kekulize(perception: &mut ChemicalPerception) -> Result<(), PerceptionError> {
    let mut visited_bonds = vec![false; perception.bonds.len()];
    let mut total_attempts = 0;

    for bond_idx in 0..perception.bonds.len() {
        if perception.bonds[bond_idx].is_aromatic && !visited_bonds[bond_idx] {
            let component_bond_indices =
                collect_aromatic_component(perception, bond_idx, &mut visited_bonds);

            assign_kekule_orders(perception, &component_bond_indices, &mut total_attempts)?;
        }
    }
    Ok(())
}

fn collect_aromatic_component(
    perception: &ChemicalPerception,
    start_bond_idx: usize,
    visited_bonds: &mut [bool],
) -> Vec<usize> {
    let mut component_indices = Vec::new();
    let mut queue = VecDeque::new();

    visited_bonds[start_bond_idx] = true;
    queue.push_back(start_bond_idx);

    while let Some(bond_idx) = queue.pop_front() {
        component_indices.push(bond_idx);
        let bond = &perception.bonds[bond_idx];

        for atom_id in [bond.start_atom_id, bond.end_atom_id] {
            if let Some(&atom_idx) = perception.atom_id_to_index.get(&atom_id) {
                for &(_, neighbor_bond_id) in &perception.adjacency[atom_idx] {
                    if let Some(&neighbor_bond_idx) =
                        perception.bond_id_to_index.get(&neighbor_bond_id)
                    {
                        if !perception.bonds[neighbor_bond_idx].is_aromatic
                            || visited_bonds[neighbor_bond_idx]
                        {
                            continue;
                        }
                        visited_bonds[neighbor_bond_idx] = true;
                        queue.push_back(neighbor_bond_idx);
                    }
                }
            }
        }
    }
    component_indices
}

fn assign_kekule_orders(
    perception: &mut ChemicalPerception,
    component_bond_indices: &[usize],
    total_attempts: &mut usize,
) -> Result<(), PerceptionError> {
    let unassigned_bond_indices: Vec<usize> = component_bond_indices
        .iter()
        .copied()
        .filter(|&idx| perception.bonds[idx].kekule_order.is_none())
        .collect();

    if unassigned_bond_indices.is_empty() {
        return Ok(());
    }

    let mut atom_double_bond_counts = HashMap::new();

    if !kekule_backtrack(
        perception,
        0,
        &unassigned_bond_indices,
        &mut atom_double_bond_counts,
        total_attempts,
    ) {
        return Err(PerceptionError::KekulizationFailed(*total_attempts));
    }

    for &bond_idx in component_bond_indices {
        if perception.bonds[bond_idx].kekule_order.is_none() {
            perception.bonds[bond_idx].kekule_order = Some(BondOrder::Single);
        }
    }

    Ok(())
}

fn kekule_backtrack(
    perception: &mut ChemicalPerception,
    position: usize,
    unassigned_bonds: &[usize],
    atom_counts: &mut HashMap<AtomId, u8>,
    attempts: &mut usize,
) -> bool {
    if position == unassigned_bonds.len() {
        return true;
    }

    if *attempts >= KEKULIZATION_ATTEMPT_LIMIT {
        return false;
    }
    *attempts += 1;

    let bond_idx = unassigned_bonds[position];
    let (start_id, end_id) = {
        let bond = &perception.bonds[bond_idx];
        (bond.start_atom_id, bond.end_atom_id)
    };

    let can_assign_double = atom_counts.get(&start_id).copied().unwrap_or(0) == 0
        && atom_counts.get(&end_id).copied().unwrap_or(0) == 0;

    if can_assign_double {
        perception.bonds[bond_idx].kekule_order = Some(BondOrder::Double);
        *atom_counts.entry(start_id).or_insert(0) += 1;
        *atom_counts.entry(end_id).or_insert(0) += 1;

        if kekule_backtrack(
            perception,
            position + 1,
            unassigned_bonds,
            atom_counts,
            attempts,
        ) {
            return true;
        }

        *atom_counts.get_mut(&start_id).unwrap() -= 1;
        *atom_counts.get_mut(&end_id).unwrap() -= 1;
        perception.bonds[bond_idx].kekule_order = None;
    }

    perception.bonds[bond_idx].kekule_order = Some(BondOrder::Single);

    if kekule_backtrack(
        perception,
        position + 1,
        unassigned_bonds,
        atom_counts,
        attempts,
    ) {
        return true;
    }

    perception.bonds[bond_idx].kekule_order = None;

    false
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::atom::{AtomId, Element};
    use crate::core::bond::{BondId, BondOrder};
    use crate::molecule::Molecule;
    use std::collections::HashSet;

    fn add_atoms(molecule: &mut Molecule, specs: &[(Element, i8)]) -> Vec<AtomId> {
        specs
            .iter()
            .map(|(element, charge)| molecule.add_atom(*element, *charge))
            .collect()
    }

    fn add_ring_bond(
        molecule: &mut Molecule,
        atoms: &[AtomId],
        start: usize,
        end: usize,
        order: BondOrder,
        ring_bonds: &mut Vec<BondId>,
    ) {
        let bond_id = molecule
            .add_bond(atoms[start], atoms[end], order)
            .expect("failed to add ring bond");
        ring_bonds.push(bond_id);
    }

    fn perceive_and_kekulize(molecule: &Molecule) -> ChemicalPerception {
        let mut perception = ChemicalPerception::from_graph(molecule).expect("perception failed");
        kekulize(&mut perception).expect("kekulization failed");
        perception
    }

    fn verify_kekule_assignments(
        perception: &ChemicalPerception,
        ring_atoms: &[AtomId],
        ring_bonds: &[BondId],
        expected_double_bonds: usize,
    ) {
        let ring_bond_set: HashSet<_> = ring_bonds.iter().copied().collect();

        let double_bond_count = ring_bonds
            .iter()
            .filter(|&&bond_id| {
                let idx = perception.bond_id_to_index[&bond_id];
                perception.bonds[idx].kekule_order == Some(BondOrder::Double)
            })
            .count();
        assert_eq!(double_bond_count, expected_double_bonds);

        for &bond_id in ring_bonds {
            let idx = perception.bond_id_to_index[&bond_id];
            assert!(
                perception.bonds[idx].kekule_order.is_some(),
                "bond {} missing kekule assignment",
                bond_id
            );
        }

        for &atom_id in ring_atoms {
            let idx = perception.atom_id_to_index[&atom_id];
            let double_incident = perception.adjacency[idx]
                .iter()
                .filter(|(_, bond_id)| ring_bond_set.contains(bond_id))
                .filter(|(_, bond_id)| {
                    let bond_idx = perception.bond_id_to_index[bond_id];
                    perception.bonds[bond_idx].kekule_order == Some(BondOrder::Double)
                })
                .count();
            assert_eq!(
                double_incident, 1,
                "atom {} has {} double bonds",
                atom_id, double_incident
            );
        }
    }

    #[test]
    fn benzene_kekulization_assigns_alternating_bonds() {
        let mut molecule = Molecule::new();
        let atoms = add_atoms(&mut molecule, &vec![(Element::C, 0); 6]);
        let mut ring_bonds = Vec::new();

        add_ring_bond(
            &mut molecule,
            &atoms,
            0,
            1,
            BondOrder::Double,
            &mut ring_bonds,
        );
        add_ring_bond(
            &mut molecule,
            &atoms,
            1,
            2,
            BondOrder::Single,
            &mut ring_bonds,
        );
        add_ring_bond(
            &mut molecule,
            &atoms,
            2,
            3,
            BondOrder::Double,
            &mut ring_bonds,
        );
        add_ring_bond(
            &mut molecule,
            &atoms,
            3,
            4,
            BondOrder::Single,
            &mut ring_bonds,
        );
        add_ring_bond(
            &mut molecule,
            &atoms,
            4,
            5,
            BondOrder::Double,
            &mut ring_bonds,
        );
        add_ring_bond(
            &mut molecule,
            &atoms,
            5,
            0,
            BondOrder::Single,
            &mut ring_bonds,
        );

        let perception = perceive_and_kekulize(&molecule);
        verify_kekule_assignments(&perception, &atoms, &ring_bonds, 3);
    }

    #[test]
    fn naphthalene_kekulization_assigns_valid_pattern() {
        let mut molecule = Molecule::new();
        let atoms = add_atoms(&mut molecule, &vec![(Element::C, 0); 10]);
        let mut ring_bonds = Vec::new();

        add_ring_bond(
            &mut molecule,
            &atoms,
            0,
            1,
            BondOrder::Double,
            &mut ring_bonds,
        );
        add_ring_bond(
            &mut molecule,
            &atoms,
            1,
            2,
            BondOrder::Single,
            &mut ring_bonds,
        );
        add_ring_bond(
            &mut molecule,
            &atoms,
            2,
            3,
            BondOrder::Double,
            &mut ring_bonds,
        );
        add_ring_bond(
            &mut molecule,
            &atoms,
            3,
            4,
            BondOrder::Single,
            &mut ring_bonds,
        );
        add_ring_bond(
            &mut molecule,
            &atoms,
            4,
            5,
            BondOrder::Double,
            &mut ring_bonds,
        );
        add_ring_bond(
            &mut molecule,
            &atoms,
            5,
            0,
            BondOrder::Single,
            &mut ring_bonds,
        );

        add_ring_bond(
            &mut molecule,
            &atoms,
            5,
            6,
            BondOrder::Single,
            &mut ring_bonds,
        );
        add_ring_bond(
            &mut molecule,
            &atoms,
            6,
            7,
            BondOrder::Double,
            &mut ring_bonds,
        );
        add_ring_bond(
            &mut molecule,
            &atoms,
            7,
            8,
            BondOrder::Single,
            &mut ring_bonds,
        );
        add_ring_bond(
            &mut molecule,
            &atoms,
            8,
            9,
            BondOrder::Double,
            &mut ring_bonds,
        );
        add_ring_bond(
            &mut molecule,
            &atoms,
            9,
            4,
            BondOrder::Single,
            &mut ring_bonds,
        );

        let perception = perceive_and_kekulize(&molecule);
        verify_kekule_assignments(&perception, &atoms, &ring_bonds, 5);
    }

    #[test]
    fn pyridine_kekulization_assigns_valid_pattern() {
        let mut molecule = Molecule::new();
        let atoms = add_atoms(
            &mut molecule,
            &vec![
                (Element::C, 0),
                (Element::C, 0),
                (Element::C, 0),
                (Element::C, 0),
                (Element::C, 0),
                (Element::N, 0),
            ],
        );
        let mut ring_bonds = Vec::new();

        add_ring_bond(
            &mut molecule,
            &atoms,
            0,
            1,
            BondOrder::Double,
            &mut ring_bonds,
        );
        add_ring_bond(
            &mut molecule,
            &atoms,
            1,
            2,
            BondOrder::Single,
            &mut ring_bonds,
        );
        add_ring_bond(
            &mut molecule,
            &atoms,
            2,
            3,
            BondOrder::Double,
            &mut ring_bonds,
        );
        add_ring_bond(
            &mut molecule,
            &atoms,
            3,
            4,
            BondOrder::Single,
            &mut ring_bonds,
        );
        add_ring_bond(
            &mut molecule,
            &atoms,
            4,
            5,
            BondOrder::Double,
            &mut ring_bonds,
        );
        add_ring_bond(
            &mut molecule,
            &atoms,
            5,
            0,
            BondOrder::Single,
            &mut ring_bonds,
        );

        let perception = perceive_and_kekulize(&molecule);
        verify_kekule_assignments(&perception, &atoms, &ring_bonds, 3);
    }

    #[test]
    fn biphenyl_kekulization_handles_multiple_components() {
        let mut molecule = Molecule::new();
        let atoms = add_atoms(&mut molecule, &vec![(Element::C, 0); 12]);
        let mut ring_bonds = Vec::new();

        add_ring_bond(
            &mut molecule,
            &atoms,
            0,
            1,
            BondOrder::Double,
            &mut ring_bonds,
        );
        add_ring_bond(
            &mut molecule,
            &atoms,
            1,
            2,
            BondOrder::Single,
            &mut ring_bonds,
        );
        add_ring_bond(
            &mut molecule,
            &atoms,
            2,
            3,
            BondOrder::Double,
            &mut ring_bonds,
        );
        add_ring_bond(
            &mut molecule,
            &atoms,
            3,
            4,
            BondOrder::Single,
            &mut ring_bonds,
        );
        add_ring_bond(
            &mut molecule,
            &atoms,
            4,
            5,
            BondOrder::Double,
            &mut ring_bonds,
        );
        add_ring_bond(
            &mut molecule,
            &atoms,
            5,
            0,
            BondOrder::Single,
            &mut ring_bonds,
        );

        add_ring_bond(
            &mut molecule,
            &atoms,
            6,
            7,
            BondOrder::Double,
            &mut ring_bonds,
        );
        add_ring_bond(
            &mut molecule,
            &atoms,
            7,
            8,
            BondOrder::Single,
            &mut ring_bonds,
        );
        add_ring_bond(
            &mut molecule,
            &atoms,
            8,
            9,
            BondOrder::Double,
            &mut ring_bonds,
        );
        add_ring_bond(
            &mut molecule,
            &atoms,
            9,
            10,
            BondOrder::Single,
            &mut ring_bonds,
        );
        add_ring_bond(
            &mut molecule,
            &atoms,
            10,
            11,
            BondOrder::Double,
            &mut ring_bonds,
        );
        add_ring_bond(
            &mut molecule,
            &atoms,
            11,
            6,
            BondOrder::Single,
            &mut ring_bonds,
        );

        molecule
            .add_bond(atoms[5], atoms[6], BondOrder::Single)
            .expect("failed to link phenyl rings");

        let perception = perceive_and_kekulize(&molecule);
        verify_kekule_assignments(&perception, &atoms, &ring_bonds, 6);
    }

    #[test]
    fn purine_kekulization_assigns_valid_pattern() {
        let mut molecule = Molecule::new();

        let atom_specs = vec![
            (Element::N, 0),
            (Element::C, 0),
            (Element::N, 0),
            (Element::C, 0),
            (Element::C, 0),
            (Element::C, 0),
            (Element::N, 0),
            (Element::C, 0),
            (Element::N, 0),
        ];
        let atoms = add_atoms(&mut molecule, &atom_specs);
        let mut ring_bonds = Vec::new();

        add_ring_bond(
            &mut molecule,
            &atoms,
            0,
            1,
            BondOrder::Aromatic,
            &mut ring_bonds,
        );
        add_ring_bond(
            &mut molecule,
            &atoms,
            1,
            2,
            BondOrder::Aromatic,
            &mut ring_bonds,
        );
        add_ring_bond(
            &mut molecule,
            &atoms,
            2,
            3,
            BondOrder::Aromatic,
            &mut ring_bonds,
        );
        add_ring_bond(
            &mut molecule,
            &atoms,
            3,
            4,
            BondOrder::Aromatic,
            &mut ring_bonds,
        );
        add_ring_bond(
            &mut molecule,
            &atoms,
            4,
            5,
            BondOrder::Aromatic,
            &mut ring_bonds,
        );
        add_ring_bond(
            &mut molecule,
            &atoms,
            5,
            0,
            BondOrder::Aromatic,
            &mut ring_bonds,
        );
        add_ring_bond(
            &mut molecule,
            &atoms,
            4,
            6,
            BondOrder::Aromatic,
            &mut ring_bonds,
        );
        add_ring_bond(
            &mut molecule,
            &atoms,
            6,
            7,
            BondOrder::Aromatic,
            &mut ring_bonds,
        );
        add_ring_bond(
            &mut molecule,
            &atoms,
            7,
            8,
            BondOrder::Aromatic,
            &mut ring_bonds,
        );
        add_ring_bond(
            &mut molecule,
            &atoms,
            8,
            3,
            BondOrder::Aromatic,
            &mut ring_bonds,
        );

        let h1 = molecule.add_atom(Element::H, 0);
        molecule
            .add_bond(atoms[6], h1, BondOrder::Single)
            .expect("failed to attach H to N7");
        let h2 = molecule.add_atom(Element::H, 0);
        molecule
            .add_bond(atoms[8], h2, BondOrder::Single)
            .expect("failed to attach H to N9");

        let perception = perceive_and_kekulize(&molecule);

        let double_bond_count = ring_bonds
            .iter()
            .filter(|&&bond_id| {
                let idx = perception.bond_id_to_index[&bond_id];
                perception.bonds[idx].kekule_order == Some(BondOrder::Double)
            })
            .count();
        assert_eq!(
            double_bond_count, 4,
            "Purine should be Kekulized with 4 double bonds"
        );

        for &atom_id in &atoms {
            let idx = perception.atom_id_to_index[&atom_id];
            let double_incident = perception.adjacency[idx]
                .iter()
                .filter(|(_, bond_id)| ring_bonds.contains(bond_id))
                .filter(|(_, bond_id)| {
                    let bond_idx = perception.bond_id_to_index[bond_id];
                    perception.bonds[bond_idx].kekule_order == Some(BondOrder::Double)
                })
                .count();

            assert!(
                double_incident <= 1,
                "atom {} has more than one double bond",
                atom_id
            );
        }
    }
}
