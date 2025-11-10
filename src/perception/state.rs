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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::atom::{AtomId, Element};
    use crate::core::bond::{BondId, BondOrder};
    use crate::molecule::Molecule;
    use crate::perception::ChemicalPerception;

    fn atom_index(perception: &ChemicalPerception, atom_id: AtomId) -> usize {
        perception.atom_id_to_index[&atom_id]
    }

    fn add_ring_bond(
        molecule: &mut Molecule,
        atoms: &[AtomId],
        start: usize,
        end: usize,
        order: BondOrder,
        bonds: &mut Vec<BondId>,
    ) {
        let bond_id = molecule
            .add_bond(atoms[start], atoms[end], order)
            .expect("failed to add ring bond");
        bonds.push(bond_id);
    }

    fn build_ethane() -> (ChemicalPerception, [AtomId; 2]) {
        let mut molecule = Molecule::new();
        let c0 = molecule.add_atom(Element::C, 0);
        let c1 = molecule.add_atom(Element::C, 0);
        let hydrogens: Vec<AtomId> = (0..6).map(|_| molecule.add_atom(Element::H, 0)).collect();

        molecule
            .add_bond(c0, c1, BondOrder::Single)
            .expect("failed to add C-C bond");
        for &h in &hydrogens[..3] {
            molecule
                .add_bond(c0, h, BondOrder::Single)
                .expect("failed to add C-H bond");
        }
        for &h in &hydrogens[3..] {
            molecule
                .add_bond(c1, h, BondOrder::Single)
                .expect("failed to add C-H bond");
        }

        let perception = ChemicalPerception::from_graph(&molecule).expect("perception failed");
        (perception, [c0, c1])
    }

    fn build_ethene() -> (ChemicalPerception, [AtomId; 2]) {
        let mut molecule = Molecule::new();
        let c0 = molecule.add_atom(Element::C, 0);
        let c1 = molecule.add_atom(Element::C, 0);
        let hydrogens: Vec<AtomId> = (0..4).map(|_| molecule.add_atom(Element::H, 0)).collect();

        molecule
            .add_bond(c0, c1, BondOrder::Double)
            .expect("failed to add C=C bond");
        for &h in &hydrogens[..2] {
            molecule
                .add_bond(c0, h, BondOrder::Single)
                .expect("failed to add C-H bond");
        }
        for &h in &hydrogens[2..] {
            molecule
                .add_bond(c1, h, BondOrder::Single)
                .expect("failed to add C-H bond");
        }

        let perception = ChemicalPerception::from_graph(&molecule).expect("perception failed");
        (perception, [c0, c1])
    }

    fn build_ethyne() -> (ChemicalPerception, [AtomId; 2]) {
        let mut molecule = Molecule::new();
        let c0 = molecule.add_atom(Element::C, 0);
        let c1 = molecule.add_atom(Element::C, 0);
        let h0 = molecule.add_atom(Element::H, 0);
        let h1 = molecule.add_atom(Element::H, 0);

        molecule
            .add_bond(c0, c1, BondOrder::Triple)
            .expect("failed to add C≡C bond");
        molecule
            .add_bond(c0, h0, BondOrder::Single)
            .expect("failed to add C-H bond");
        molecule
            .add_bond(c1, h1, BondOrder::Single)
            .expect("failed to add C-H bond");

        let perception = ChemicalPerception::from_graph(&molecule).expect("perception failed");
        (perception, [c0, c1])
    }

    fn build_water() -> (ChemicalPerception, AtomId) {
        let mut molecule = Molecule::new();
        let o = molecule.add_atom(Element::O, 0);
        let h0 = molecule.add_atom(Element::H, 0);
        let h1 = molecule.add_atom(Element::H, 0);

        molecule
            .add_bond(o, h0, BondOrder::Single)
            .expect("failed to add O-H bond");
        molecule
            .add_bond(o, h1, BondOrder::Single)
            .expect("failed to add O-H bond");

        let perception = ChemicalPerception::from_graph(&molecule).expect("perception failed");
        (perception, o)
    }

    fn build_ammonia() -> (ChemicalPerception, AtomId) {
        let mut molecule = Molecule::new();
        let n = molecule.add_atom(Element::N, 0);
        for _ in 0..3 {
            let h = molecule.add_atom(Element::H, 0);
            molecule
                .add_bond(n, h, BondOrder::Single)
                .expect("failed to add N-H bond");
        }

        let perception = ChemicalPerception::from_graph(&molecule).expect("perception failed");
        (perception, n)
    }

    fn build_methane() -> (ChemicalPerception, AtomId) {
        let mut molecule = Molecule::new();
        let c = molecule.add_atom(Element::C, 0);
        for _ in 0..4 {
            let h = molecule.add_atom(Element::H, 0);
            molecule
                .add_bond(c, h, BondOrder::Single)
                .expect("failed to add C-H bond");
        }

        let perception = ChemicalPerception::from_graph(&molecule).expect("perception failed");
        (perception, c)
    }

    fn build_hydroxide() -> (ChemicalPerception, AtomId) {
        let mut molecule = Molecule::new();
        let o = molecule.add_atom(Element::O, -1);
        let h = molecule.add_atom(Element::H, 0);

        molecule
            .add_bond(o, h, BondOrder::Single)
            .expect("failed to add O-H bond");

        let perception = ChemicalPerception::from_graph(&molecule).expect("perception failed");
        (perception, o)
    }

    fn build_hydronium() -> (ChemicalPerception, AtomId) {
        let mut molecule = Molecule::new();
        let o = molecule.add_atom(Element::O, 1);
        for _ in 0..3 {
            let h = molecule.add_atom(Element::H, 0);
            molecule
                .add_bond(o, h, BondOrder::Single)
                .expect("failed to add O-H bond");
        }

        let perception = ChemicalPerception::from_graph(&molecule).expect("perception failed");
        (perception, o)
    }

    fn add_benzene_core(molecule: &mut Molecule) -> Vec<AtomId> {
        let carbons: Vec<AtomId> = (0..6).map(|_| molecule.add_atom(Element::C, 0)).collect();
        let mut bonds = Vec::new();
        add_ring_bond(molecule, &carbons, 0, 1, BondOrder::Double, &mut bonds);
        add_ring_bond(molecule, &carbons, 1, 2, BondOrder::Single, &mut bonds);
        add_ring_bond(molecule, &carbons, 2, 3, BondOrder::Double, &mut bonds);
        add_ring_bond(molecule, &carbons, 3, 4, BondOrder::Single, &mut bonds);
        add_ring_bond(molecule, &carbons, 4, 5, BondOrder::Double, &mut bonds);
        add_ring_bond(molecule, &carbons, 5, 0, BondOrder::Single, &mut bonds);
        carbons
    }

    fn build_benzene() -> (ChemicalPerception, Vec<AtomId>) {
        let mut molecule = Molecule::new();
        let carbons = add_benzene_core(&mut molecule);
        for &carbon in &carbons {
            let h = molecule.add_atom(Element::H, 0);
            molecule
                .add_bond(carbon, h, BondOrder::Single)
                .expect("failed to add C-H bond");
        }

        let perception = ChemicalPerception::from_graph(&molecule).expect("perception failed");
        (perception, carbons)
    }

    fn build_aniline() -> (ChemicalPerception, AtomId) {
        let mut molecule = Molecule::new();
        let carbons = add_benzene_core(&mut molecule);
        for &carbon in &carbons[1..] {
            let h = molecule.add_atom(Element::H, 0);
            molecule
                .add_bond(carbon, h, BondOrder::Single)
                .expect("failed to add ring C-H bond");
        }

        let nitrogen = molecule.add_atom(Element::N, 0);
        molecule
            .add_bond(carbons[0], nitrogen, BondOrder::Single)
            .expect("failed to connect aniline nitrogen");
        for _ in 0..2 {
            let h = molecule.add_atom(Element::H, 0);
            molecule
                .add_bond(nitrogen, h, BondOrder::Single)
                .expect("failed to add N-H bond");
        }

        let perception = ChemicalPerception::from_graph(&molecule).expect("perception failed");
        (perception, nitrogen)
    }

    fn build_phenol() -> (ChemicalPerception, AtomId) {
        let mut molecule = Molecule::new();
        let carbons = add_benzene_core(&mut molecule);
        for &carbon in &carbons[1..] {
            let h = molecule.add_atom(Element::H, 0);
            molecule
                .add_bond(carbon, h, BondOrder::Single)
                .expect("failed to add ring C-H bond");
        }

        let oxygen = molecule.add_atom(Element::O, 0);
        molecule
            .add_bond(carbons[0], oxygen, BondOrder::Single)
            .expect("failed to connect phenol oxygen");
        let h = molecule.add_atom(Element::H, 0);
        molecule
            .add_bond(oxygen, h, BondOrder::Single)
            .expect("failed to add O-H bond");

        let perception = ChemicalPerception::from_graph(&molecule).expect("perception failed");
        (perception, oxygen)
    }

    fn build_pyrrole() -> (ChemicalPerception, AtomId) {
        let mut molecule = Molecule::new();
        let c0 = molecule.add_atom(Element::C, 0);
        let c1 = molecule.add_atom(Element::C, 0);
        let c2 = molecule.add_atom(Element::C, 0);
        let c3 = molecule.add_atom(Element::C, 0);
        let n = molecule.add_atom(Element::N, 0);
        let atoms = [c0, c1, c2, c3, n];
        let mut bonds = Vec::new();
        add_ring_bond(&mut molecule, &atoms, 0, 1, BondOrder::Double, &mut bonds);
        add_ring_bond(&mut molecule, &atoms, 1, 2, BondOrder::Single, &mut bonds);
        add_ring_bond(&mut molecule, &atoms, 2, 3, BondOrder::Double, &mut bonds);
        add_ring_bond(&mut molecule, &atoms, 3, 4, BondOrder::Single, &mut bonds);
        add_ring_bond(&mut molecule, &atoms, 4, 0, BondOrder::Single, &mut bonds);

        for &carbon in &atoms[..4] {
            let h = molecule.add_atom(Element::H, 0);
            molecule
                .add_bond(carbon, h, BondOrder::Single)
                .expect("failed to add C-H bond");
        }
        let h = molecule.add_atom(Element::H, 0);
        molecule
            .add_bond(n, h, BondOrder::Single)
            .expect("failed to add N-H bond");

        let perception = ChemicalPerception::from_graph(&molecule).expect("perception failed");
        (perception, n)
    }

    #[test]
    fn total_valence_matches_expected_across_bond_orders_and_charges() {
        let (perception, carbons) = build_ethane();
        for carbon in carbons {
            let idx = atom_index(&perception, carbon);
            assert_eq!(
                perception.atoms[idx].total_valence, 4,
                "ethane carbon valence"
            );
        }

        let (perception, carbons) = build_ethene();
        for carbon in carbons {
            let idx = atom_index(&perception, carbon);
            assert_eq!(
                perception.atoms[idx].total_valence, 4,
                "ethene carbon valence"
            );
        }

        let (perception, carbons) = build_ethyne();
        for carbon in carbons {
            let idx = atom_index(&perception, carbon);
            assert_eq!(
                perception.atoms[idx].total_valence, 4,
                "ethyne carbon valence"
            );
        }

        let (perception, oxygen) = build_hydronium();
        let idx = atom_index(&perception, oxygen);
        assert_eq!(
            perception.atoms[idx].total_valence, 3,
            "hydronium oxygen valence"
        );
    }

    #[test]
    fn lone_pairs_are_estimated_correctly_for_common_species() {
        let (perception, oxygen) = build_water();
        let idx = atom_index(&perception, oxygen);
        assert_eq!(
            perception.atoms[idx].lone_pairs, 2,
            "water oxygen lone pairs"
        );

        let (perception, nitrogen) = build_ammonia();
        let idx = atom_index(&perception, nitrogen);
        assert_eq!(
            perception.atoms[idx].lone_pairs, 1,
            "ammonia nitrogen lone pairs"
        );

        let (perception, carbon) = build_methane();
        let idx = atom_index(&perception, carbon);
        assert_eq!(
            perception.atoms[idx].lone_pairs, 0,
            "methane carbon lone pairs"
        );

        let (perception, oxygen) = build_hydroxide();
        let idx = atom_index(&perception, oxygen);
        assert_eq!(
            perception.atoms[idx].lone_pairs, 3,
            "hydroxide oxygen lone pairs"
        );

        let (perception, oxygen) = build_hydronium();
        let idx = atom_index(&perception, oxygen);
        assert_eq!(
            perception.atoms[idx].lone_pairs, 1,
            "hydronium oxygen lone pairs"
        );
    }

    #[test]
    fn hybridization_matches_expected_for_simple_geometries() {
        let (perception, carbon) = build_methane();
        let idx = atom_index(&perception, carbon);
        assert_eq!(perception.atoms[idx].hybridization, Hybridization::SP3);

        let (perception, carbons) = build_ethene();
        for carbon in carbons {
            let idx = atom_index(&perception, carbon);
            assert_eq!(perception.atoms[idx].hybridization, Hybridization::SP2);
        }

        let (perception, carbons) = build_ethyne();
        for carbon in carbons {
            let idx = atom_index(&perception, carbon);
            assert_eq!(perception.atoms[idx].hybridization, Hybridization::SP);
        }

        let (perception, oxygen) = build_water();
        let idx = atom_index(&perception, oxygen);
        assert_eq!(perception.atoms[idx].hybridization, Hybridization::SP3);
    }

    #[test]
    fn aromatic_atoms_are_treated_as_sp2() {
        let (perception, carbons) = build_benzene();
        let idx = atom_index(&perception, carbons[0]);
        assert!(
            perception.atoms[idx].is_aromatic,
            "expected aromatic carbon"
        );
        assert_eq!(perception.atoms[idx].hybridization, Hybridization::SP2);
    }

    #[test]
    fn conjugation_correction_promotes_heteroatoms_to_sp2() {
        let (perception, nitrogen) = build_aniline();
        let n_idx = atom_index(&perception, nitrogen);
        assert_eq!(perception.atoms[n_idx].hybridization, Hybridization::SP2);
        assert_eq!(
            perception.atoms[n_idx].lone_pairs, 1,
            "aniline nitrogen lone pairs"
        );

        let aromatic_neighbor_sp2 = perception.adjacency[n_idx]
            .iter()
            .any(|&(neighbor_idx, _)| {
                perception.atoms[neighbor_idx].hybridization == Hybridization::SP2
            });
        assert!(
            aromatic_neighbor_sp2,
            "aniline nitrogen should see an sp2 neighbor"
        );

        let (perception, oxygen) = build_phenol();
        let o_idx = atom_index(&perception, oxygen);
        assert_eq!(perception.atoms[o_idx].hybridization, Hybridization::SP2);
        assert_eq!(
            perception.atoms[o_idx].lone_pairs, 2,
            "phenol oxygen lone pairs"
        );

        let (perception, nitrogen) = build_pyrrole();
        let n_idx = atom_index(&perception, nitrogen);
        assert_eq!(perception.atoms[n_idx].hybridization, Hybridization::SP2);
        assert!(
            perception.atoms[n_idx].is_aromatic,
            "pyrrole nitrogen should be aromatic"
        );
    }
}
