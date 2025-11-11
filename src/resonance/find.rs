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

    for (bond_idx, bond) in perception.bonds.iter().enumerate() {
        let effective_order = bond.kekule_order.unwrap_or(bond.order);
        if bond.is_aromatic || matches!(effective_order, BondOrder::Double | BondOrder::Triple) {
            conjugated.insert(bond_idx);
        }
    }

    let mut bonds_to_add_in_pass: Vec<usize> = Vec::new();
    loop {
        for atom_idx in 0..perception.atoms.len() {
            let atom = &perception.atoms[atom_idx];
            if !atom.is_conjugation_candidate {
                continue;
            }

            let is_connected_to_system =
                perception.adjacency[atom_idx].iter().any(|&(_, bond_id)| {
                    let bond_idx = perception.bond_id_to_index[&bond_id];
                    conjugated.contains(&bond_idx)
                });

            if !is_connected_to_system {
                continue;
            }

            for &(_, bond_id) in &perception.adjacency[atom_idx] {
                let bond_idx = perception.bond_id_to_index[&bond_id];
                if conjugated.contains(&bond_idx) {
                    continue;
                }

                let bond = &perception.bonds[bond_idx];
                let other_end_id = bond.other_end(atom.id);
                let other_end_idx = perception.atom_id_to_index[&other_end_id];

                if perception.atoms[other_end_idx].is_conjugation_candidate {
                    bonds_to_add_in_pass.push(bond_idx);
                }
            }
        }

        if bonds_to_add_in_pass.is_empty() {
            break;
        }

        for bond_idx in bonds_to_add_in_pass.drain(..) {
            conjugated.insert(bond_idx);
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
