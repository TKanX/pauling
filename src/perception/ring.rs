use crate::core::atom::AtomId;
use crate::core::bond::BondId;
use crate::perception::ChemicalPerception;
use std::collections::{HashMap, HashSet, VecDeque};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Ring {
    pub atom_ids: Vec<AtomId>,
    pub bond_ids: Vec<BondId>,
}

impl Ring {
    pub fn new(mut atom_ids: Vec<AtomId>, mut bond_ids: Vec<BondId>) -> Self {
        atom_ids.sort_unstable();
        atom_ids.dedup();
        bond_ids.sort_unstable();
        bond_ids.dedup();
        Self { atom_ids, bond_ids }
    }
}

#[derive(Debug, Clone, Default)]
pub struct RingInfo {
    pub rings: Vec<Ring>,
}

pub fn find_sssr(perception: &ChemicalPerception) -> RingInfo {
    let num_components = count_components(perception);
    let cyclomatic_number =
        perception.bonds.len() as isize - perception.atoms.len() as isize + num_components as isize;

    if cyclomatic_number <= 0 {
        return RingInfo::default();
    }

    let candidates = enumerate_cycle_candidates(perception);
    let selected_rings =
        select_minimal_cycle_basis(perception, candidates, cyclomatic_number as usize);

    RingInfo {
        rings: selected_rings,
    }
}

struct PathData {
    atom_ids: Vec<AtomId>,
    bond_ids: Vec<BondId>,
}

fn enumerate_cycle_candidates(perception: &ChemicalPerception) -> Vec<Ring> {
    let mut candidates = Vec::new();
    let mut seen_signatures: HashSet<Vec<BondId>> = HashSet::new();

    for bond in &perception.bonds {
        if let Some(path) =
            shortest_path_excluding_bond(perception, bond.start_atom_id, bond.end_atom_id, bond.id)
        {
            let mut all_bond_ids = path.bond_ids;
            all_bond_ids.push(bond.id);

            all_bond_ids.sort_unstable();
            if seen_signatures.insert(all_bond_ids.clone()) {
                candidates.push(Ring::new(path.atom_ids, all_bond_ids));
            }
        }
    }
    candidates
}

fn shortest_path_excluding_bond(
    perception: &ChemicalPerception,
    start_atom_id: AtomId,
    end_atom_id: AtomId,
    forbidden_bond_id: BondId,
) -> Option<PathData> {
    let start_idx = perception.atom_id_to_index.get(&start_atom_id)?;
    let end_idx = perception.atom_id_to_index.get(&end_atom_id)?;

    let mut queue = VecDeque::new();
    let mut visited = vec![false; perception.atoms.len()];
    let mut parent: Vec<Option<(usize, BondId)>> = vec![None; perception.atoms.len()];

    visited[*start_idx] = true;
    queue.push_back(*start_idx);

    while let Some(current_idx) = queue.pop_front() {
        if current_idx == *end_idx {
            break;
        }
        for &(neighbor_idx, bond_id) in &perception.adjacency[current_idx] {
            if bond_id == forbidden_bond_id {
                continue;
            }
            if !visited[neighbor_idx] {
                visited[neighbor_idx] = true;
                parent[neighbor_idx] = Some((current_idx, bond_id));
                queue.push_back(neighbor_idx);
            }
        }
    }

    if !visited[*end_idx] {
        return None;
    }

    let mut atom_ids = Vec::new();
    let mut bond_ids = Vec::new();
    let mut cursor = *end_idx;

    while let Some((prev_idx, bond_id)) = parent[cursor] {
        atom_ids.push(perception.atoms[cursor].id);
        bond_ids.push(bond_id);
        cursor = prev_idx;
    }
    atom_ids.push(perception.atoms[cursor].id);

    atom_ids.reverse();
    bond_ids.reverse();

    Some(PathData { atom_ids, bond_ids })
}

fn select_minimal_cycle_basis(
    perception: &ChemicalPerception,
    mut candidates: Vec<Ring>,
    cyclomatic_number: usize,
) -> Vec<Ring> {
    candidates.sort_by_key(|r| r.bond_ids.len());

    let mut selected_rings = Vec::new();
    let mut basis: Vec<(BitVec, usize)> = Vec::new();

    for ring in candidates {
        let mut bitvec = BitVec::from_bond_ids(&ring.bond_ids, &perception.bond_id_to_index);

        for (basis_vec, pivot) in &basis {
            if bitvec.test(*pivot) {
                bitvec.xor(basis_vec);
            }
        }

        if !bitvec.is_zero() {
            let pivot = bitvec
                .leading_one()
                .expect("A non-zero vector must have a leading one.");

            basis.push((bitvec, pivot));
            basis.sort_by_key(|&(_, p)| p);

            selected_rings.push(ring);

            if selected_rings.len() == cyclomatic_number {
                break;
            }
        }
    }

    selected_rings
}

fn count_components(perception: &ChemicalPerception) -> usize {
    let mut visited = vec![false; perception.atoms.len()];
    let mut components = 0;

    for i in 0..perception.atoms.len() {
        if !visited[i] {
            components += 1;
            let mut stack = vec![i];
            visited[i] = true;
            while let Some(current) = stack.pop() {
                for &(neighbor_idx, _) in &perception.adjacency[current] {
                    if !visited[neighbor_idx] {
                        visited[neighbor_idx] = true;
                        stack.push(neighbor_idx);
                    }
                }
            }
        }
    }
    components
}

#[derive(Clone)]
struct BitVec {
    data: Vec<u64>,
}

impl BitVec {
    fn new(size: usize) -> Self {
        let words = size.div_ceil(64);
        Self {
            data: vec![0; words],
        }
    }

    fn from_bond_ids(bond_ids: &[BondId], bond_id_to_index: &HashMap<BondId, usize>) -> Self {
        let mut bitvec = Self::new(bond_id_to_index.len());
        for bond_id in bond_ids {
            if let Some(&idx) = bond_id_to_index.get(bond_id) {
                let word = idx / 64;
                let bit = idx % 64;
                if word < bitvec.data.len() {
                    bitvec.data[word] |= 1u64 << bit;
                }
            }
        }
        bitvec
    }

    fn xor(&mut self, other: &Self) {
        debug_assert_eq!(
            self.data.len(),
            other.data.len(),
            "BitVec XOR operation requires vectors of the same size."
        );
        for (a, b) in self.data.iter_mut().zip(&other.data) {
            *a ^= *b;
        }
    }

    fn is_zero(&self) -> bool {
        self.data.iter().all(|&word| word == 0)
    }

    fn test(&self, idx: usize) -> bool {
        let word_idx = idx / 64;
        let bit_idx = idx % 64;
        if word_idx >= self.data.len() {
            return false;
        }
        (self.data[word_idx] & (1u64 << bit_idx)) != 0
    }

    fn leading_one(&self) -> Option<usize> {
        for (word_idx, &word) in self.data.iter().enumerate().rev() {
            if word != 0 {
                let leading_zeros = word.leading_zeros() as usize;
                let bit_idx_in_word = 63 - leading_zeros;
                return Some(word_idx * 64 + bit_idx_in_word);
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::atom::Element;
    use crate::core::bond::BondOrder;
    use crate::perception::{ChemicalPerception, PerceivedAtom, PerceivedBond};
    use std::collections::HashMap;

    fn build_perception(edges: &[(BondId, AtomId, AtomId)]) -> ChemicalPerception {
        let all_atom_ids: HashSet<AtomId> = edges
            .iter()
            .flat_map(|(_, start, end)| [*start, *end])
            .collect();

        let mut sorted_atom_ids: Vec<AtomId> = all_atom_ids.into_iter().collect();
        sorted_atom_ids.sort_unstable();

        let num_atoms = sorted_atom_ids.len();

        let atom_id_to_index: HashMap<AtomId, usize> = sorted_atom_ids
            .iter()
            .enumerate()
            .map(|(index, &id)| (id, index))
            .collect();

        let mut adjacency: Vec<Vec<(usize, BondId)>> = vec![Vec::new(); num_atoms];
        let mut bonds = Vec::new();
        let mut bond_id_to_index = HashMap::new();

        for (bond_index, &(bond_id, start_atom_id, end_atom_id)) in edges.iter().enumerate() {
            let start_idx = atom_id_to_index[&start_atom_id];
            let end_idx = atom_id_to_index[&end_atom_id];

            adjacency[start_idx].push((end_idx, bond_id));
            adjacency[end_idx].push((start_idx, bond_id));

            bond_id_to_index.insert(bond_id, bond_index);
            bonds.push(PerceivedBond::new(
                bond_id,
                BondOrder::Single,
                start_atom_id,
                end_atom_id,
            ));
        }

        let atoms: Vec<PerceivedAtom> = sorted_atom_ids
            .iter()
            .map(|&id| {
                let idx = atom_id_to_index[&id];
                PerceivedAtom::new(id, Element::C, 0, adjacency[idx].len() as u8)
            })
            .collect();

        ChemicalPerception {
            atoms,
            bonds,
            adjacency,
            atom_id_to_index,
            bond_id_to_index,
            ring_info: RingInfo::default(),
        }
    }

    #[test]
    fn find_sssr_returns_no_rings_for_acyclic_graph() {
        let perception = build_perception(&[(0, 0, 1), (1, 1, 2)]);
        let ring_info = find_sssr(&perception);
        assert!(ring_info.rings.is_empty());
    }

    #[test]
    fn find_sssr_detects_single_cycle() {
        let perception = build_perception(&[(0, 0, 1), (1, 1, 2), (2, 2, 3), (3, 3, 0)]);

        let ring_info = find_sssr(&perception);
        assert_eq!(ring_info.rings.len(), 1);

        let ring = &ring_info.rings[0];
        assert_eq!(ring.atom_ids, vec![0, 1, 2, 3]);
        assert_eq!(ring.bond_ids, vec![0, 1, 2, 3]);
    }

    #[test]
    fn find_sssr_identifies_multiple_independent_cycles() {
        let perception = build_perception(&[
            (0, 0, 1),
            (1, 1, 2),
            (2, 2, 3),
            (3, 3, 0),
            (4, 3, 4),
            (5, 4, 5),
            (6, 5, 6),
            (7, 6, 3),
        ]);

        let ring_info = find_sssr(&perception);
        assert_eq!(ring_info.rings.len(), 2);

        let mut atom_sets: Vec<Vec<AtomId>> = ring_info
            .rings
            .iter()
            .map(|ring| ring.atom_ids.clone())
            .collect();
        atom_sets.sort();
        assert_eq!(atom_sets, vec![vec![0, 1, 2, 3], vec![3, 4, 5, 6]]);

        let mut bond_sets: Vec<Vec<BondId>> = ring_info
            .rings
            .iter()
            .map(|ring| ring.bond_ids.clone())
            .collect();
        bond_sets.sort();
        assert_eq!(bond_sets, vec![vec![0, 1, 2, 3], vec![4, 5, 6, 7]]);
    }

    #[test]
    fn find_sssr_handles_disconnected_components_with_cycles() {
        let perception = build_perception(&[
            (0, 0, 1),
            (1, 1, 2),
            (2, 2, 0),
            (3, 3, 4),
            (4, 4, 5),
            (5, 5, 3),
        ]);

        let ring_info = find_sssr(&perception);
        assert_eq!(ring_info.rings.len(), 2);

        let mut atom_sets: Vec<Vec<AtomId>> = ring_info
            .rings
            .iter()
            .map(|ring| ring.atom_ids.clone())
            .collect();
        atom_sets.sort();
        assert_eq!(atom_sets, vec![vec![0, 1, 2], vec![3, 4, 5]]);

        let mut bond_sets: Vec<Vec<BondId>> = ring_info
            .rings
            .iter()
            .map(|ring| ring.bond_ids.clone())
            .collect();
        bond_sets.sort();
        assert_eq!(bond_sets, vec![vec![0, 1, 2], vec![3, 4, 5]]);
    }
}
