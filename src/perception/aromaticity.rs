use crate::core::atom::Element;
use crate::core::bond::BondOrder;
use crate::perception::ChemicalPerception;
use std::collections::{HashMap, HashSet};

pub fn perceive(perception: &mut ChemicalPerception) {
    apply_explicit_aromaticity(perception);
    apply_topological_aromaticity(perception);
}

/// Phase 1: Handles bonds explicitly marked as `BondOrder::Aromatic`.
/// This step "trusts the input" and seeds the aromaticity perception.
fn apply_explicit_aromaticity(perception: &mut ChemicalPerception) {
    let aromatic_bonds: Vec<usize> = perception
        .bonds
        .iter()
        .enumerate()
        .filter(|(_, bond)| bond.order == BondOrder::Aromatic)
        .map(|(idx, _)| idx)
        .collect();

    for bond_idx in aromatic_bonds {
        let bond = &mut perception.bonds[bond_idx];
        if bond.is_aromatic {
            continue;
        } // Already processed
        bond.is_aromatic = true;

        if let Some(start_idx) = perception.atom_id_to_index.get(&bond.start_atom_id) {
            perception.atoms[*start_idx].is_aromatic = true;
        }
        if let Some(end_idx) = perception.atom_id_to_index.get(&bond.end_atom_id) {
            perception.atoms[*end_idx].is_aromatic = true;
        }
    }
}

/// Phase 2: Detects aromaticity based on topology and Hückel's rule.
fn apply_topological_aromaticity(perception: &mut ChemicalPerception) {
    if perception.ring_info.rings.is_empty() {
        return;
    }

    // Step 2a: Group rings into fused systems.
    let fused_systems = find_fused_ring_systems(perception);

    // Step 2b: Check each fused system for aromaticity.
    for system_indices in fused_systems {
        if is_system_aromatic(perception, &system_indices) {
            // If aromatic, mark all atoms and bonds in the system.
            let mut all_atom_indices = HashSet::new();
            let mut all_bond_indices = HashSet::new();
            for &ring_idx in &system_indices {
                let ring = &perception.ring_info.rings[ring_idx];
                for &atom_id in &ring.atom_ids {
                    if let Some(&idx) = perception.atom_id_to_index.get(&atom_id) {
                        all_atom_indices.insert(idx);
                    }
                }
                for &bond_id in &ring.bond_ids {
                    if let Some(&idx) = perception.bond_id_to_index.get(&bond_id) {
                        all_bond_indices.insert(idx);
                    }
                }
            }

            for idx in all_atom_indices {
                perception.atoms[idx].is_aromatic = true;
            }
            for idx in all_bond_indices {
                perception.bonds[idx].is_aromatic = true;
            }
        }
    }
}

/// Groups rings into connected components based on shared bonds.
fn find_fused_ring_systems(perception: &ChemicalPerception) -> Vec<Vec<usize>> {
    let num_rings = perception.ring_info.rings.len();
    if num_rings == 0 {
        return Vec::new();
    }

    // Create a map from bond ID to the rings it belongs to.
    let mut bond_to_rings = HashMap::new();
    for (ring_idx, ring) in perception.ring_info.rings.iter().enumerate() {
        for &bond_id in &ring.bond_ids {
            bond_to_rings
                .entry(bond_id)
                .or_insert_with(Vec::new)
                .push(ring_idx);
        }
    }

    // Build an adjacency graph of rings.
    let mut ring_adjacency = vec![Vec::new(); num_rings];
    for rings in bond_to_rings.values() {
        if rings.len() > 1 {
            for i in 0..rings.len() {
                for j in (i + 1)..rings.len() {
                    let r1 = rings[i];
                    let r2 = rings[j];
                    ring_adjacency[r1].push(r2);
                    ring_adjacency[r2].push(r1);
                }
            }
        }
    }

    // Find connected components in the ring graph (these are the fused systems).
    let mut visited = vec![false; num_rings];
    let mut components = Vec::new();
    for i in 0..num_rings {
        if !visited[i] {
            let mut component = Vec::new();
            let mut stack = vec![i];
            visited[i] = true;
            while let Some(current) = stack.pop() {
                component.push(current);
                for &neighbor in &ring_adjacency[current] {
                    if !visited[neighbor] {
                        visited[neighbor] = true;
                        stack.push(neighbor);
                    }
                }
            }
            components.push(component);
        }
    }
    components
}

/// Checks if a single fused ring system is aromatic using Hückel's rule.
fn is_system_aromatic(perception: &ChemicalPerception, system_ring_indices: &[usize]) -> bool {
    // Collect all unique atom and bond indices in the system.
    let mut system_atom_indices = HashSet::new();
    let mut system_bond_indices = HashSet::new();

    for &ring_idx in system_ring_indices {
        let ring = &perception.ring_info.rings[ring_idx];
        for &atom_id in &ring.atom_ids {
            system_atom_indices.insert(perception.atom_id_to_index[&atom_id]);
        }
        for &bond_id in &ring.bond_ids {
            system_bond_indices.insert(perception.bond_id_to_index[&bond_id]);
        }
    }

    // An atom must be a potential sp2 hybrid to participate in an aromatic system.
    for &atom_idx in &system_atom_indices {
        if !is_potential_sp2_hybrid(perception, atom_idx) {
            return false;
        }
    }

    // Sum π electrons contributed by each atom in the system.
    let mut pi_electron_count = 0;
    for &atom_idx in &system_atom_indices {
        pi_electron_count += pi_electrons_for_atom(perception, atom_idx, &system_bond_indices);
    }

    // Apply Hückel's rule: 4n + 2 π electrons.
    pi_electron_count > 0 && (pi_electron_count - 2) % 4 == 0
}

/// A heuristic check if an atom can adopt sp2 hybridization for aromaticity.
fn is_potential_sp2_hybrid(perception: &ChemicalPerception, atom_idx: usize) -> bool {
    let atom = &perception.atoms[atom_idx];
    // This rule covers most common cases in organic chemistry.
    // Transition metals and hypervalent atoms are out of scope.
    atom.total_degree <= 3
}

/// Encapsulates the chemical rules for counting π electrons contributed by an atom.
fn pi_electrons_for_atom(
    perception: &ChemicalPerception,
    atom_idx: usize,
    system_bond_indices: &HashSet<usize>,
) -> u32 {
    let atom = &perception.atoms[atom_idx];

    // Case 1: Atom is part of a multiple bond within the ring system.
    // It contributes 1 π electron. This is the most common case (e.g., C in benzene).
    let is_in_multiple_bond_in_system =
        perception.adjacency[atom_idx].iter().any(|&(_, bond_id)| {
            if let Some(&bond_idx) = perception.bond_id_to_index.get(&bond_id) {
                if system_bond_indices.contains(&bond_idx) {
                    let bond = &perception.bonds[bond_idx];
                    matches!(bond.order, BondOrder::Double | BondOrder::Triple)
                } else {
                    false
                }
            } else {
                false
            }
        });

    if is_in_multiple_bond_in_system {
        return 1;
    }

    // Case 2: Atom is NOT part of a multiple bond, contributes via lone pair or empty orbital.
    // These rules are based on common patterns in heterocycles.
    match atom.element {
        // Pyrrole-like Nitrogen
        Element::N if atom.total_degree == 3 => {
            if atom.formal_charge == 1 { 0 } else { 2 } // Positively charged N (e.g., in protonated indole) has no lone pair to donate.
        }
        // Furan-like Oxygen or Thiophene-like Sulfur
        Element::O | Element::S if atom.total_degree == 2 => 2,
        // Carbocation or Carbanion in a ring
        Element::C if atom.total_degree == 3 => match atom.formal_charge {
            -1 => 2, // Carbanion (e.g., cyclopentadienyl anion)
            1 => 0,  // Carbocation (e.g., tropylium cation)
            _ => 0,
        },
        // Boron in a ring
        Element::B if atom.total_degree == 3 => 0, // Has an empty p-orbital, contributes 0 electrons.
        _ => 0,                                    // Default case: atom does not contribute.
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::atom::{AtomId, Element};
    use crate::core::bond::{BondId, BondOrder};
    use crate::molecule::Molecule;

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

    fn assert_aromatic_state(
        perception: &ChemicalPerception,
        atom_ids: &[AtomId],
        bond_ids: &[BondId],
        expected: bool,
    ) {
        for &atom_id in atom_ids {
            let idx = perception.atom_id_to_index[&atom_id];
            assert_eq!(
                perception.atoms[idx].is_aromatic, expected,
                "atom {} aromatic flag mismatch",
                atom_id
            );
        }

        for &bond_id in bond_ids {
            let idx = perception.bond_id_to_index[&bond_id];
            assert_eq!(
                perception.bonds[idx].is_aromatic, expected,
                "bond {} aromatic flag mismatch",
                bond_id
            );
        }
    }

    #[test]
    fn benzene_is_aromatic() {
        let mut molecule = Molecule::new();
        let atom_specs = vec![(Element::C, 0); 6];
        let atoms = add_atoms(&mut molecule, &atom_specs);
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

        let perception = ChemicalPerception::from_graph(&molecule).expect("perception failed");
        assert_aromatic_state(&perception, &atoms, &ring_bonds, true);
    }

    #[test]
    fn cyclohexane_is_non_aromatic() {
        let mut molecule = Molecule::new();
        let atom_specs = vec![(Element::C, 0); 6];
        let atoms = add_atoms(&mut molecule, &atom_specs);
        let mut ring_bonds = Vec::new();

        add_ring_bond(
            &mut molecule,
            &atoms,
            0,
            1,
            BondOrder::Single,
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
            BondOrder::Single,
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
            BondOrder::Single,
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

        let perception = ChemicalPerception::from_graph(&molecule).expect("perception failed");
        assert_aromatic_state(&perception, &atoms, &ring_bonds, false);
    }

    #[test]
    fn cyclobutadiene_is_non_aromatic() {
        let mut molecule = Molecule::new();
        let atom_specs = vec![(Element::C, 0); 4];
        let atoms = add_atoms(&mut molecule, &atom_specs);
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
            0,
            BondOrder::Single,
            &mut ring_bonds,
        );

        let perception = ChemicalPerception::from_graph(&molecule).expect("perception failed");
        assert_aromatic_state(&perception, &atoms, &ring_bonds, false);
    }

    #[test]
    fn pyridine_is_aromatic() {
        let mut molecule = Molecule::new();
        let atom_specs = vec![
            (Element::C, 0),
            (Element::C, 0),
            (Element::C, 0),
            (Element::C, 0),
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

        let perception = ChemicalPerception::from_graph(&molecule).expect("perception failed");
        assert_aromatic_state(&perception, &atoms, &ring_bonds, true);
    }

    #[test]
    fn pyrrole_is_aromatic() {
        let mut molecule = Molecule::new();
        let atom_specs = vec![
            (Element::C, 0),
            (Element::C, 0),
            (Element::C, 0),
            (Element::C, 0),
            (Element::N, 0),
        ];
        let atoms = add_atoms(&mut molecule, &atom_specs);
        let ring_atoms = atoms.clone();
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
            0,
            BondOrder::Single,
            &mut ring_bonds,
        );

        let hydrogen = molecule.add_atom(Element::H, 0);
        molecule
            .add_bond(atoms[4], hydrogen, BondOrder::Single)
            .expect("failed to attach hydrogen");

        let perception = ChemicalPerception::from_graph(&molecule).expect("perception failed");
        assert_aromatic_state(&perception, &ring_atoms, &ring_bonds, true);
    }

    #[test]
    fn furan_is_aromatic() {
        let mut molecule = Molecule::new();
        let atom_specs = vec![
            (Element::C, 0),
            (Element::C, 0),
            (Element::C, 0),
            (Element::C, 0),
            (Element::O, 0),
        ];
        let atoms = add_atoms(&mut molecule, &atom_specs);
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
            0,
            BondOrder::Single,
            &mut ring_bonds,
        );

        let perception = ChemicalPerception::from_graph(&molecule).expect("perception failed");
        assert_aromatic_state(&perception, &atoms, &ring_bonds, true);
    }

    #[test]
    fn imidazole_is_aromatic() {
        let mut molecule = Molecule::new();
        let atom_specs = vec![
            (Element::C, 0),
            (Element::N, 0),
            (Element::C, 0),
            (Element::C, 0),
            (Element::N, 0),
        ];
        let atoms = add_atoms(&mut molecule, &atom_specs);
        let ring_atoms = atoms.clone();
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
            0,
            BondOrder::Single,
            &mut ring_bonds,
        );

        let hydrogen = molecule.add_atom(Element::H, 0);
        molecule
            .add_bond(atoms[4], hydrogen, BondOrder::Single)
            .expect("failed to attach hydrogen");

        let perception = ChemicalPerception::from_graph(&molecule).expect("perception failed");
        assert_aromatic_state(&perception, &ring_atoms, &ring_bonds, true);
    }

    #[test]
    fn naphthalene_is_aromatic() {
        let mut molecule = Molecule::new();
        let atom_specs = vec![(Element::C, 0); 10];
        let atoms = add_atoms(&mut molecule, &atom_specs);
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

        let perception = ChemicalPerception::from_graph(&molecule).expect("perception failed");
        assert_eq!(
            ring_bonds.len(),
            11,
            "expected 11 ring bonds for naphthalene"
        );
        assert_aromatic_state(&perception, &atoms, &ring_bonds, true);
    }

    #[test]
    fn indole_is_aromatic() {
        let mut molecule = Molecule::new();
        let atom_specs = vec![
            (Element::C, 0),
            (Element::C, 0),
            (Element::C, 0),
            (Element::C, 0),
            (Element::C, 0),
            (Element::C, 0),
            (Element::C, 0),
            (Element::N, 0),
            (Element::C, 0),
        ];
        let atoms = add_atoms(&mut molecule, &atom_specs);
        let ring_atoms = atoms.clone();
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
            BondOrder::Single,
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
            BondOrder::Double,
            &mut ring_bonds,
        );
        add_ring_bond(
            &mut molecule,
            &atoms,
            6,
            7,
            BondOrder::Single,
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
            4,
            BondOrder::Double,
            &mut ring_bonds,
        );

        let hydrogen = molecule.add_atom(Element::H, 0);
        molecule
            .add_bond(atoms[7], hydrogen, BondOrder::Single)
            .expect("failed to attach hydrogen");

        let perception = ChemicalPerception::from_graph(&molecule).expect("perception failed");
        assert_eq!(ring_bonds.len(), 10, "expected 10 ring bonds for indole");
        assert_aromatic_state(&perception, &ring_atoms, &ring_bonds, true);
    }

    #[test]
    fn cyclopentadienyl_anion_is_aromatic() {
        let mut molecule = Molecule::new();
        let atom_specs = vec![
            (Element::C, -1),
            (Element::C, 0),
            (Element::C, 0),
            (Element::C, 0),
            (Element::C, 0),
        ];
        let atoms = add_atoms(&mut molecule, &atom_specs);
        let ring_atoms = atoms.clone();
        let mut ring_bonds = Vec::new();

        add_ring_bond(
            &mut molecule,
            &atoms,
            0,
            1,
            BondOrder::Single,
            &mut ring_bonds,
        );
        add_ring_bond(
            &mut molecule,
            &atoms,
            1,
            2,
            BondOrder::Double,
            &mut ring_bonds,
        );
        add_ring_bond(
            &mut molecule,
            &atoms,
            2,
            3,
            BondOrder::Single,
            &mut ring_bonds,
        );
        add_ring_bond(
            &mut molecule,
            &atoms,
            3,
            4,
            BondOrder::Double,
            &mut ring_bonds,
        );
        add_ring_bond(
            &mut molecule,
            &atoms,
            4,
            0,
            BondOrder::Single,
            &mut ring_bonds,
        );

        let hydrogen = molecule.add_atom(Element::H, 0);
        molecule
            .add_bond(atoms[0], hydrogen, BondOrder::Single)
            .expect("failed to attach hydrogen");

        let perception = ChemicalPerception::from_graph(&molecule).expect("perception failed");
        assert_aromatic_state(&perception, &ring_atoms, &ring_bonds, true);
    }

    #[test]
    fn tropylium_cation_is_aromatic() {
        let mut molecule = Molecule::new();
        let atom_specs = vec![
            (Element::C, 1),
            (Element::C, 0),
            (Element::C, 0),
            (Element::C, 0),
            (Element::C, 0),
            (Element::C, 0),
            (Element::C, 0),
        ];
        let atoms = add_atoms(&mut molecule, &atom_specs);
        let ring_atoms = atoms.clone();
        let mut ring_bonds = Vec::new();

        add_ring_bond(
            &mut molecule,
            &atoms,
            0,
            1,
            BondOrder::Single,
            &mut ring_bonds,
        );
        add_ring_bond(
            &mut molecule,
            &atoms,
            1,
            2,
            BondOrder::Double,
            &mut ring_bonds,
        );
        add_ring_bond(
            &mut molecule,
            &atoms,
            2,
            3,
            BondOrder::Single,
            &mut ring_bonds,
        );
        add_ring_bond(
            &mut molecule,
            &atoms,
            3,
            4,
            BondOrder::Double,
            &mut ring_bonds,
        );
        add_ring_bond(
            &mut molecule,
            &atoms,
            4,
            5,
            BondOrder::Single,
            &mut ring_bonds,
        );
        add_ring_bond(
            &mut molecule,
            &atoms,
            5,
            6,
            BondOrder::Double,
            &mut ring_bonds,
        );
        add_ring_bond(
            &mut molecule,
            &atoms,
            6,
            0,
            BondOrder::Single,
            &mut ring_bonds,
        );

        let hydrogen = molecule.add_atom(Element::H, 0);
        molecule
            .add_bond(atoms[0], hydrogen, BondOrder::Single)
            .expect("failed to attach hydrogen");

        let perception = ChemicalPerception::from_graph(&molecule).expect("perception failed");
        assert_aromatic_state(&perception, &ring_atoms, &ring_bonds, true);
    }
}
