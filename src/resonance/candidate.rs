use crate::core::atom::Element;
use crate::perception::ChemicalPerception;
use crate::perception::Hybridization;

pub fn determine(perception: &mut ChemicalPerception) {
    for atom in &mut perception.atoms {
        atom.is_conjugation_candidate =
            matches!(atom.hybridization, Hybridization::SP | Hybridization::SP2);
    }

    for atom_idx in 0..perception.atoms.len() {
        if perception.atoms[atom_idx].is_conjugation_candidate {
            continue;
        }

        if is_lone_pair_candidate(perception, atom_idx)
            || is_charged_carbon_candidate(perception, atom_idx)
        {
            perception.atoms[atom_idx].is_conjugation_candidate = true;
        }
    }
}

fn is_lone_pair_candidate(perception: &ChemicalPerception, atom_idx: usize) -> bool {
    let atom = &perception.atoms[atom_idx];

    if atom.lone_pairs > 0 {
        for &(neighbor_idx, _) in &perception.adjacency[atom_idx] {
            if perception.atoms[neighbor_idx].is_conjugation_candidate {
                return true;
            }
        }
    }

    false
}

fn is_charged_carbon_candidate(perception: &ChemicalPerception, atom_idx: usize) -> bool {
    let atom = &perception.atoms[atom_idx];

    if atom.element != Element::C {
        return false;
    }

    match atom.formal_charge {
        1 if atom.total_degree == 3 => true,
        -1 => true,
        _ => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::atom::{AtomId, Element};
    use crate::core::bond::BondOrder;
    use crate::molecule::Molecule;
    use crate::perception::ChemicalPerception;

    fn index(perception: &ChemicalPerception, atom_id: AtomId) -> usize {
        perception.atom_id_to_index[&atom_id]
    }

    fn build_ethene() -> (ChemicalPerception, Vec<AtomId>) {
        let mut molecule = Molecule::new();
        let c0 = molecule.add_atom(Element::C, 0);
        let c1 = molecule.add_atom(Element::C, 0);
        let hydrogens: Vec<AtomId> = (0..4).map(|_| molecule.add_atom(Element::H, 0)).collect();

        molecule
            .add_bond(c0, c1, BondOrder::Double)
            .expect("add C=C");
        for &h in &hydrogens[..2] {
            molecule
                .add_bond(c0, h, BondOrder::Single)
                .expect("add C-H");
        }
        for &h in &hydrogens[2..] {
            molecule
                .add_bond(c1, h, BondOrder::Single)
                .expect("add C-H");
        }

        let mut perception = ChemicalPerception::from_graph(&molecule).expect("perception");
        determine(&mut perception);
        (perception, vec![c0, c1])
    }

    fn build_ethane() -> (ChemicalPerception, Vec<AtomId>) {
        let mut molecule = Molecule::new();
        let c0 = molecule.add_atom(Element::C, 0);
        let c1 = molecule.add_atom(Element::C, 0);
        let hydrogens: Vec<AtomId> = (0..6).map(|_| molecule.add_atom(Element::H, 0)).collect();

        molecule
            .add_bond(c0, c1, BondOrder::Single)
            .expect("add C-C");
        for &h in &hydrogens[..3] {
            molecule
                .add_bond(c0, h, BondOrder::Single)
                .expect("add C-H");
        }
        for &h in &hydrogens[3..] {
            molecule
                .add_bond(c1, h, BondOrder::Single)
                .expect("add C-H");
        }

        let mut perception = ChemicalPerception::from_graph(&molecule).expect("perception");
        determine(&mut perception);
        (perception, vec![c0, c1])
    }

    fn build_acetamide() -> (ChemicalPerception, AtomId) {
        let mut molecule = Molecule::new();
        let carbonyl_c = molecule.add_atom(Element::C, 0);
        let oxygen = molecule.add_atom(Element::O, 0);
        let nitrogen = molecule.add_atom(Element::N, 0);
        let methyl_carbon = molecule.add_atom(Element::C, 0);

        molecule
            .add_bond(carbonyl_c, oxygen, BondOrder::Double)
            .expect("carbonyl");
        molecule
            .add_bond(carbonyl_c, nitrogen, BondOrder::Single)
            .expect("amide N");
        molecule
            .add_bond(carbonyl_c, methyl_carbon, BondOrder::Single)
            .expect("methyl bond");

        for _ in 0..3 {
            let h = molecule.add_atom(Element::H, 0);
            molecule
                .add_bond(methyl_carbon, h, BondOrder::Single)
                .expect("methyl H");
        }

        for _ in 0..2 {
            let h = molecule.add_atom(Element::H, 0);
            molecule
                .add_bond(nitrogen, h, BondOrder::Single)
                .expect("amide NH");
        }

        let mut perception = ChemicalPerception::from_graph(&molecule).expect("perception");
        determine(&mut perception);
        (perception, nitrogen)
    }

    fn build_aniline() -> (ChemicalPerception, AtomId) {
        let mut molecule = Molecule::new();
        let ring: Vec<AtomId> = (0..6).map(|_| molecule.add_atom(Element::C, 0)).collect();
        let nitrogen = molecule.add_atom(Element::N, 0);

        let bond_pairs = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)];
        for (idx, &(a, b)) in bond_pairs.iter().enumerate() {
            let order = if idx % 2 == 0 {
                BondOrder::Double
            } else {
                BondOrder::Single
            };
            molecule
                .add_bond(ring[a], ring[b], order)
                .expect("aromatic bond");
        }

        for &carbon in &ring[1..] {
            let h = molecule.add_atom(Element::H, 0);
            molecule
                .add_bond(carbon, h, BondOrder::Single)
                .expect("aryl H");
        }

        molecule
            .add_bond(ring[0], nitrogen, BondOrder::Single)
            .expect("aniline N");
        for _ in 0..2 {
            let h = molecule.add_atom(Element::H, 0);
            molecule
                .add_bond(nitrogen, h, BondOrder::Single)
                .expect("amine H");
        }

        let mut perception = ChemicalPerception::from_graph(&molecule).expect("perception");
        determine(&mut perception);
        (perception, nitrogen)
    }

    fn build_cyclopentadienyl_anion() -> (ChemicalPerception, Vec<AtomId>) {
        let mut molecule = Molecule::new();
        let atoms: Vec<AtomId> = (0..5)
            .map(|i| {
                if i == 0 {
                    molecule.add_atom(Element::C, -1)
                } else {
                    molecule.add_atom(Element::C, 0)
                }
            })
            .collect();

        let orders = [
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
        ];
        for i in 0..5 {
            let next = (i + 1) % 5;
            molecule
                .add_bond(atoms[i], atoms[next], orders[i])
                .expect("ring bond");
        }

        for &carbon in &atoms {
            let h = molecule.add_atom(Element::H, 0);
            molecule
                .add_bond(carbon, h, BondOrder::Single)
                .expect("ring hydrogen");
        }

        let mut perception = ChemicalPerception::from_graph(&molecule).expect("perception");
        determine(&mut perception);
        (perception, atoms)
    }

    fn build_allyl_cation() -> (ChemicalPerception, AtomId) {
        let mut molecule = Molecule::new();
        let c_plus = molecule.add_atom(Element::C, 1);
        let c1 = molecule.add_atom(Element::C, 0);
        let c2 = molecule.add_atom(Element::C, 0);

        molecule
            .add_bond(c_plus, c1, BondOrder::Single)
            .expect("C+-C");
        molecule.add_bond(c1, c2, BondOrder::Double).expect("C=C");

        for &carbon in &[c_plus, c1, c2] {
            let hydrogens = if carbon == c1 { 1 } else { 2 };
            for _ in 0..hydrogens {
                let h = molecule.add_atom(Element::H, 0);
                molecule
                    .add_bond(carbon, h, BondOrder::Single)
                    .expect("C-H");
            }
        }

        let mut perception = ChemicalPerception::from_graph(&molecule).expect("perception");
        determine(&mut perception);
        (perception, c_plus)
    }

    fn build_dimethyl_ether() -> (ChemicalPerception, AtomId) {
        let mut molecule = Molecule::new();
        let c0 = molecule.add_atom(Element::C, 0);
        let c1 = molecule.add_atom(Element::C, 0);
        let oxygen = molecule.add_atom(Element::O, 0);

        molecule
            .add_bond(c0, oxygen, BondOrder::Single)
            .expect("C-O");
        molecule
            .add_bond(c1, oxygen, BondOrder::Single)
            .expect("C-O");

        for &carbon in &[c0, c1] {
            for _ in 0..3 {
                let h = molecule.add_atom(Element::H, 0);
                molecule
                    .add_bond(carbon, h, BondOrder::Single)
                    .expect("C-H");
            }
        }

        let mut perception = ChemicalPerception::from_graph(&molecule).expect("perception");
        determine(&mut perception);
        (perception, oxygen)
    }

    #[test]
    fn sp2_carbons_from_multiple_bonds_are_candidates() {
        let (perception, carbons) = build_ethene();
        for carbon in carbons {
            let idx = index(&perception, carbon);
            assert!(
                perception.atoms[idx].is_conjugation_candidate,
                "ethene carbon"
            );
        }
    }

    #[test]
    fn saturated_carbons_are_not_candidates() {
        let (perception, carbons) = build_ethane();
        for carbon in carbons {
            let idx = index(&perception, carbon);
            assert!(
                !perception.atoms[idx].is_conjugation_candidate,
                "ethane carbon"
            );
        }
    }

    #[test]
    fn lone_pair_on_adjacent_atom_promotes_candidate_status() {
        let (perception, nitrogen) = build_acetamide();
        let idx = index(&perception, nitrogen);
        assert!(
            perception.atoms[idx].is_conjugation_candidate,
            "acetamide nitrogen"
        );

        let (perception, nitrogen) = build_aniline();
        let idx = index(&perception, nitrogen);
        assert!(
            perception.atoms[idx].is_conjugation_candidate,
            "aniline nitrogen"
        );
    }

    #[test]
    fn charged_carbons_and_delocalized_anions_become_candidates() {
        let (perception, atoms) = build_cyclopentadienyl_anion();
        for atom in atoms {
            let idx = index(&perception, atom);
            assert!(
                perception.atoms[idx].is_conjugation_candidate,
                "cyclopentadienyl carbon"
            );
        }

        let (perception, c_plus) = build_allyl_cation();
        let idx = index(&perception, c_plus);
        assert!(
            perception.atoms[idx].is_conjugation_candidate,
            "allyl cation"
        );
    }

    #[test]
    fn non_conjugating_oxygen_remains_non_candidate() {
        let (perception, oxygen) = build_dimethyl_ether();
        let idx = index(&perception, oxygen);
        assert!(
            !perception.atoms[idx].is_conjugation_candidate,
            "dimethyl ether oxygen"
        );
    }
}
