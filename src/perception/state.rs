use crate::perception::{ChemicalPerception, PerceivedAtom};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Hybridization {
    SP,
    SP2,
    SP3,
    Unknown,
}

pub fn perceive(perception: &mut ChemicalPerception) {
    compute_valence(perception);
    perceive_hybridization(perception);
}

fn compute_valence(perception: &mut ChemicalPerception) {
    for atom in &mut perception.atoms {
        atom.total_valence = 0;
    }

    for bond in &perception.bonds {
        let effective_order = bond.kekule_order.unwrap_or(bond.order);
        let multiplicity = effective_order.multiplicity();

        if let Some(&start_idx) = perception.atom_id_to_index.get(&bond.start_atom_id) {
            perception.atoms[start_idx].total_valence = perception.atoms[start_idx]
                .total_valence
                .saturating_add(multiplicity);
        }
        if let Some(&end_idx) = perception.atom_id_to_index.get(&bond.end_atom_id) {
            perception.atoms[end_idx].total_valence = perception.atoms[end_idx]
                .total_valence
                .saturating_add(multiplicity);
        }
    }
}

fn perceive_hybridization(perception: &mut ChemicalPerception) {
    let lone_pairs: Vec<u8> = perception.atoms.iter().map(estimate_lone_pairs).collect();

    for (idx, atom) in perception.atoms.iter_mut().enumerate() {
        atom.lone_pairs = lone_pairs[idx];

        if atom.is_aromatic {
            atom.hybridization = Hybridization::SP2;
            continue;
        }

        let steric_number = atom.total_degree.saturating_add(atom.lone_pairs);
        atom.hybridization = match steric_number {
            2 => Hybridization::SP,
            3 => Hybridization::SP2,
            4 => Hybridization::SP3,
            _ => Hybridization::Unknown,
        };
    }

    let snapshot: Vec<Hybridization> = perception
        .atoms
        .iter()
        .map(|atom| atom.hybridization)
        .collect();

    for (idx, atom) in perception.atoms.iter_mut().enumerate() {
        if atom.hybridization == Hybridization::SP3 && atom.lone_pairs > 0 {
            let has_sp2_or_sp_neighbor =
                perception.adjacency[idx].iter().any(|&(neighbor_idx, _)| {
                    matches!(
                        snapshot[neighbor_idx],
                        Hybridization::SP | Hybridization::SP2
                    )
                });

            if has_sp2_or_sp_neighbor {
                atom.hybridization = Hybridization::SP2;
            }
        }
    }
}

fn estimate_lone_pairs(atom: &PerceivedAtom) -> u8 {
    let valence_electrons = match atom.element.valence_electrons() {
        Some(e) => e as i16,
        None => return 0,
    };

    let non_bonding_electrons =
        valence_electrons - (atom.formal_charge as i16) - (atom.total_valence as i16);

    (non_bonding_electrons.max(0) / 2) as u8
}
