//! Determines which atoms can participate in conjugation and resonance.

use crate::core::atom::Element;
use crate::core::bond::BondOrder;
use crate::perception::{ChemicalPerception, ConjugationRole, Hybridization};

/// Marks atoms as conjugation candidates based on hybridization, charge, and
/// hypervalent heuristics.
///
/// # Arguments
///
/// * `perception` - Fully populated perception snapshot that will receive the
///   updated `is_conjugation_candidate` flags.
pub fn determine(perception: &mut ChemicalPerception) {
    reset_conjugation_state(perception);
    mark_hypervalent_bridges(perception);
    mark_intrinsic_pi_carriers(perception);
    promote_lone_pair_donors(perception);
    promote_charged_carbons(perception);
    finalize_candidate_flags(perception);
}

/// Resets all conjugation-related flags on atoms.
fn reset_conjugation_state(perception: &mut ChemicalPerception) {
    for atom in &mut perception.atoms {
        atom.is_conjugation_candidate = false;
        atom.conjugation_roles = ConjugationRole::NONE;
    }
}

/// Marks atoms that are intrinsic pi carriers, such as sp2/sp carbons and
/// aromatic atoms.
fn mark_intrinsic_pi_carriers(perception: &mut ChemicalPerception) {
    for atom_idx in 0..perception.atoms.len() {
        if should_skip_intrinsic_pi(perception, atom_idx) {
            continue;
        }

        if perception.atoms[atom_idx].is_aromatic
            || matches!(
                perception.atoms[atom_idx].hybridization,
                Hybridization::SP | Hybridization::SP2
            )
        {
            perception.atoms[atom_idx]
                .conjugation_roles
                .insert(ConjugationRole::PI_CARRIER);
        }
    }
}

/// Identifies hypervalent bridge atoms and marks them accordingly.
fn mark_hypervalent_bridges(perception: &mut ChemicalPerception) {
    for atom_idx in 0..perception.atoms.len() {
        if is_hypervalent_bridge(perception, atom_idx) {
            perception.atoms[atom_idx]
                .conjugation_roles
                .insert(ConjugationRole::HYPERVALENT_BRIDGE);
        }
    }
}

/// Promotes atoms with lone pairs that are adjacent to conjugation-capable
/// atoms to conjugation candidates.
fn promote_lone_pair_donors(perception: &mut ChemicalPerception) {
    for atom_idx in 0..perception.atoms.len() {
        let atom = &perception.atoms[atom_idx];
        if atom.lone_pairs == 0 {
            continue;
        }

        let mut promote = false;

        for &(neighbor_idx, _) in &perception.adjacency[atom_idx] {
            let neighbour = &perception.atoms[neighbor_idx];
            let neighbour_roles = neighbour.conjugation_roles;

            if neighbour_roles.is_empty() {
                continue;
            }

            if neighbour_roles.contains(ConjugationRole::HYPERVALENT_BRIDGE) {
                if atom.formal_charge < 0 {
                    promote = true;
                    break;
                }
                continue;
            }

            promote = true;
            break;
        }

        if !promote {
            continue;
        }

        if atom.element == Element::O && atom.formal_charge == 0 && atom.total_degree > 1 {
            continue;
        }

        perception.atoms[atom_idx]
            .conjugation_roles
            .insert(ConjugationRole::LONE_PAIR_DONOR);
    }
}

/// Promotes charged carbons to conjugation candidates based on their
/// formal charge and degree.
fn promote_charged_carbons(perception: &mut ChemicalPerception) {
    for atom in &mut perception.atoms {
        if atom.element != Element::C {
            continue;
        }

        let is_delocalised =
            (atom.formal_charge == 1 && atom.total_degree == 3) || atom.formal_charge == -1;

        if is_delocalised {
            atom.conjugation_roles
                .insert(ConjugationRole::CHARGE_MEDIATOR);
        }
    }
}

/// Finalizes the candidate flags based on the assigned conjugation roles.
fn finalize_candidate_flags(perception: &mut ChemicalPerception) {
    for atom in &mut perception.atoms {
        atom.is_conjugation_candidate = !atom.conjugation_roles.is_empty();
    }
}

/// Determines if an atom should be skipped when marking intrinsic pi carriers.
fn should_skip_intrinsic_pi(perception: &ChemicalPerception, atom_idx: usize) -> bool {
    let atom = &perception.atoms[atom_idx];

    if atom.element == Element::O && atom.formal_charge == 0 && atom.total_degree > 1 {
        return perception.adjacency[atom_idx]
            .iter()
            .any(|&(neighbor_idx, _)| {
                perception.atoms[neighbor_idx]
                    .conjugation_roles
                    .contains(ConjugationRole::HYPERVALENT_BRIDGE)
            });
    }

    false
}

/// Determines if an atom is a hypervalent bridge based on its element,
/// valence, and bonding environment.
fn is_hypervalent_bridge(perception: &ChemicalPerception, atom_idx: usize) -> bool {
    use Element::{Br, Cl, I, P, S};

    let atom = &perception.atoms[atom_idx];

    if !matches!(atom.element, P | S | Cl | Br | I) {
        return false;
    }

    if atom.total_valence <= 4 {
        return false;
    }

    let mut has_pi_partner = false;
    let mut has_sigma_partner = false;

    for &(neighbor_idx, bond_id) in &perception.adjacency[atom_idx] {
        let bond_idx = perception.bond_id_to_index[&bond_id];
        let bond = &perception.bonds[bond_idx];
        let effective_order = bond.kekule_order.unwrap_or(bond.order);

        if matches!(effective_order, BondOrder::Double | BondOrder::Triple) {
            if perception.atoms[neighbor_idx]
                .element
                .is_common_conjugation_element()
            {
                has_pi_partner = true;
            }
        } else {
            let neighbour = &perception.atoms[neighbor_idx];
            if neighbour.lone_pairs > 0
                || neighbour.formal_charge < 0
                || neighbour.element.is_common_conjugation_element()
            {
                has_sigma_partner = true;
            }
        }
    }

    has_pi_partner && has_sigma_partner
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::atom::{AtomId, Element};
    use crate::core::bond::BondOrder;
    use crate::molecule::Molecule;
    use crate::perception::{ChemicalPerception, ConjugationRole};

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

    fn build_phosphate_fragment() -> (ChemicalPerception, AtomId, AtomId, AtomId) {
        let mut molecule = Molecule::new();
        let phosphorus = molecule.add_atom(Element::P, 0);
        let _oxo = molecule.add_atom(Element::O, 0);
        molecule
            .add_bond(phosphorus, _oxo, BondOrder::Double)
            .expect("P=O");

        let anionic_oxygen = molecule.add_atom(Element::O, -1);
        molecule
            .add_bond(phosphorus, anionic_oxygen, BondOrder::Single)
            .expect("P-O-");

        let bridging_oxygen = molecule.add_atom(Element::O, 0);
        let carbon_stub = molecule.add_atom(Element::C, 0);
        molecule
            .add_bond(phosphorus, bridging_oxygen, BondOrder::Single)
            .expect("P-O (ester)");
        molecule
            .add_bond(bridging_oxygen, carbon_stub, BondOrder::Single)
            .expect("C-O stub");

        let _solvent_oxygen = molecule.add_atom(Element::O, 0);
        molecule
            .add_bond(phosphorus, _solvent_oxygen, BondOrder::Single)
            .expect("P-O (solvent)");

        let hydrogen = molecule.add_atom(Element::H, 0);
        molecule
            .add_bond(carbon_stub, hydrogen, BondOrder::Single)
            .expect("C-H");

        let mut perception = ChemicalPerception::from_graph(&molecule).expect("perception");
        determine(&mut perception);
        (perception, phosphorus, anionic_oxygen, bridging_oxygen)
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

    #[test]
    fn hypervalent_phosphate_promotes_bridge_and_neighbors() {
        let (perception, phosphorus, anionic_oxygen, bridging_oxygen) = build_phosphate_fragment();
        let p_idx = index(&perception, phosphorus);
        let anionic_idx = index(&perception, anionic_oxygen);
        let bridging_idx = index(&perception, bridging_oxygen);

        assert!(
            perception.atoms[p_idx]
                .conjugation_roles
                .contains(ConjugationRole::HYPERVALENT_BRIDGE),
            "phosphorus should register as a hypervalent bridge"
        );
        assert!(perception.atoms[p_idx].is_conjugation_candidate);
        assert!(perception.atoms[anionic_idx].is_conjugation_candidate);
        let bridging_roles = perception.atoms[bridging_idx].conjugation_roles;
        assert!(
            bridging_roles.is_empty(),
            "bridging oxygen roles: {:?}",
            bridging_roles
        );
        assert!(
            !perception.atoms[bridging_idx].is_conjugation_candidate,
            "bridging oxygen should remain outside the conjugated core"
        );
    }
}
