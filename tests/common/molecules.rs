use super::graph::TestMolecule;
use pauling::{BondOrder, Element};

pub fn build_glycine_zwitterion() -> TestMolecule {
    let mut mol = TestMolecule::new();

    mol.add_atom(0, Element::N, 1);
    mol.add_atom(1, Element::C, 0);
    mol.add_atom(2, Element::C, 0);
    mol.add_atom(3, Element::O, -1);
    mol.add_atom(4, Element::O, 0);
    mol.add_atom(5, Element::H, 0);
    mol.add_atom(6, Element::H, 0);
    mol.add_atom(7, Element::H, 0);
    mol.add_atom(8, Element::H, 0);
    mol.add_atom(9, Element::H, 0);

    mol.add_bond(0, 0, 1, BondOrder::Single);
    mol.add_bond(1, 1, 2, BondOrder::Single);
    mol.add_bond(2, 2, 3, BondOrder::Single);
    mol.add_bond(3, 2, 4, BondOrder::Double);
    mol.add_bond(4, 0, 5, BondOrder::Single);
    mol.add_bond(5, 0, 6, BondOrder::Single);
    mol.add_bond(6, 0, 7, BondOrder::Single);
    mol.add_bond(7, 1, 8, BondOrder::Single);
    mol.add_bond(8, 1, 9, BondOrder::Single);

    mol
}

pub fn build_alanine_zwitterion() -> TestMolecule {
    let mut mol = TestMolecule::new();

    mol.add_atom(0, Element::N, 1);
    mol.add_atom(1, Element::C, 0);
    mol.add_atom(2, Element::C, 0);
    mol.add_atom(3, Element::O, -1);
    mol.add_atom(4, Element::O, 0);
    mol.add_atom(5, Element::C, 0);
    mol.add_atom(6, Element::H, 0);
    mol.add_atom(7, Element::H, 0);
    mol.add_atom(8, Element::H, 0);
    mol.add_atom(9, Element::H, 0);
    mol.add_atom(10, Element::H, 0);
    mol.add_atom(11, Element::H, 0);
    mol.add_atom(12, Element::H, 0);

    mol.add_bond(0, 0, 1, BondOrder::Single);
    mol.add_bond(1, 1, 2, BondOrder::Single);
    mol.add_bond(2, 2, 3, BondOrder::Single);
    mol.add_bond(3, 2, 4, BondOrder::Double);
    mol.add_bond(4, 1, 5, BondOrder::Single);
    mol.add_bond(5, 0, 6, BondOrder::Single);
    mol.add_bond(6, 0, 7, BondOrder::Single);
    mol.add_bond(7, 0, 8, BondOrder::Single);
    mol.add_bond(8, 1, 9, BondOrder::Single);
    mol.add_bond(9, 5, 10, BondOrder::Single);
    mol.add_bond(10, 5, 11, BondOrder::Single);
    mol.add_bond(11, 5, 12, BondOrder::Single);

    mol
}

pub fn build_valine_zwitterion() -> TestMolecule {
    let mut mol = TestMolecule::new();

    mol.add_atom(0, Element::N, 1);
    mol.add_atom(1, Element::C, 0);
    mol.add_atom(2, Element::C, 0);
    mol.add_atom(3, Element::O, -1);
    mol.add_atom(4, Element::O, 0);
    mol.add_atom(5, Element::C, 0);
    mol.add_atom(6, Element::C, 0);
    mol.add_atom(7, Element::C, 0);
    mol.add_atom(8, Element::H, 0);
    mol.add_atom(9, Element::H, 0);
    mol.add_atom(10, Element::H, 0);
    mol.add_atom(11, Element::H, 0);
    mol.add_atom(12, Element::H, 0);
    mol.add_atom(13, Element::H, 0);
    mol.add_atom(14, Element::H, 0);
    mol.add_atom(15, Element::H, 0);
    mol.add_atom(16, Element::H, 0);
    mol.add_atom(17, Element::H, 0);
    mol.add_atom(18, Element::H, 0);

    mol.add_bond(0, 0, 1, BondOrder::Single);
    mol.add_bond(1, 1, 2, BondOrder::Single);
    mol.add_bond(2, 2, 3, BondOrder::Single);
    mol.add_bond(3, 2, 4, BondOrder::Double);
    mol.add_bond(4, 1, 5, BondOrder::Single);
    mol.add_bond(5, 5, 6, BondOrder::Single);
    mol.add_bond(6, 5, 7, BondOrder::Single);
    mol.add_bond(7, 0, 8, BondOrder::Single);
    mol.add_bond(8, 0, 9, BondOrder::Single);
    mol.add_bond(9, 0, 10, BondOrder::Single);
    mol.add_bond(10, 1, 11, BondOrder::Single);
    mol.add_bond(11, 5, 12, BondOrder::Single);
    mol.add_bond(12, 6, 13, BondOrder::Single);
    mol.add_bond(13, 6, 14, BondOrder::Single);
    mol.add_bond(14, 6, 15, BondOrder::Single);
    mol.add_bond(15, 7, 16, BondOrder::Single);
    mol.add_bond(16, 7, 17, BondOrder::Single);
    mol.add_bond(17, 7, 18, BondOrder::Single);

    mol
}

pub fn build_leucine_zwitterion() -> TestMolecule {
    let mut mol = TestMolecule::new();

    mol.add_atom(0, Element::N, 1);
    mol.add_atom(1, Element::C, 0);
    mol.add_atom(2, Element::C, 0);
    mol.add_atom(3, Element::O, -1);
    mol.add_atom(4, Element::O, 0);
    mol.add_atom(5, Element::C, 0);
    mol.add_atom(6, Element::C, 0);
    mol.add_atom(7, Element::C, 0);
    mol.add_atom(8, Element::C, 0);
    mol.add_atom(9, Element::H, 0);
    mol.add_atom(10, Element::H, 0);
    mol.add_atom(11, Element::H, 0);
    mol.add_atom(12, Element::H, 0);
    mol.add_atom(13, Element::H, 0);
    mol.add_atom(14, Element::H, 0);
    mol.add_atom(15, Element::H, 0);
    mol.add_atom(16, Element::H, 0);
    mol.add_atom(17, Element::H, 0);
    mol.add_atom(18, Element::H, 0);
    mol.add_atom(19, Element::H, 0);
    mol.add_atom(20, Element::H, 0);
    mol.add_atom(21, Element::H, 0);

    mol.add_bond(0, 0, 1, BondOrder::Single);
    mol.add_bond(1, 1, 2, BondOrder::Single);
    mol.add_bond(2, 2, 3, BondOrder::Single);
    mol.add_bond(3, 2, 4, BondOrder::Double);
    mol.add_bond(4, 1, 5, BondOrder::Single);
    mol.add_bond(5, 5, 6, BondOrder::Single);
    mol.add_bond(6, 6, 7, BondOrder::Single);
    mol.add_bond(7, 6, 8, BondOrder::Single);
    mol.add_bond(8, 0, 9, BondOrder::Single);
    mol.add_bond(9, 0, 10, BondOrder::Single);
    mol.add_bond(10, 0, 11, BondOrder::Single);
    mol.add_bond(11, 1, 12, BondOrder::Single);
    mol.add_bond(12, 5, 13, BondOrder::Single);
    mol.add_bond(13, 5, 14, BondOrder::Single);
    mol.add_bond(14, 6, 15, BondOrder::Single);
    mol.add_bond(15, 7, 16, BondOrder::Single);
    mol.add_bond(16, 7, 17, BondOrder::Single);
    mol.add_bond(17, 7, 18, BondOrder::Single);
    mol.add_bond(18, 8, 19, BondOrder::Single);
    mol.add_bond(19, 8, 20, BondOrder::Single);
    mol.add_bond(20, 8, 21, BondOrder::Single);

    mol
}

pub fn build_isoleucine_zwitterion() -> TestMolecule {
    let mut mol = TestMolecule::new();

    mol.add_atom(0, Element::N, 1);
    mol.add_atom(1, Element::C, 0);
    mol.add_atom(2, Element::C, 0);
    mol.add_atom(3, Element::O, -1);
    mol.add_atom(4, Element::O, 0);
    mol.add_atom(5, Element::C, 0);
    mol.add_atom(6, Element::C, 0);
    mol.add_atom(7, Element::C, 0);
    mol.add_atom(8, Element::C, 0);
    mol.add_atom(9, Element::H, 0);
    mol.add_atom(10, Element::H, 0);
    mol.add_atom(11, Element::H, 0);
    mol.add_atom(12, Element::H, 0);
    mol.add_atom(13, Element::H, 0);
    mol.add_atom(14, Element::H, 0);
    mol.add_atom(15, Element::H, 0);
    mol.add_atom(16, Element::H, 0);
    mol.add_atom(17, Element::H, 0);
    mol.add_atom(18, Element::H, 0);
    mol.add_atom(19, Element::H, 0);
    mol.add_atom(20, Element::H, 0);
    mol.add_atom(21, Element::H, 0);

    mol.add_bond(0, 0, 1, BondOrder::Single);
    mol.add_bond(1, 1, 2, BondOrder::Single);
    mol.add_bond(2, 2, 3, BondOrder::Single);
    mol.add_bond(3, 2, 4, BondOrder::Double);
    mol.add_bond(4, 1, 5, BondOrder::Single);
    mol.add_bond(5, 5, 6, BondOrder::Single);
    mol.add_bond(6, 5, 7, BondOrder::Single);
    mol.add_bond(7, 6, 8, BondOrder::Single);
    mol.add_bond(8, 0, 9, BondOrder::Single);
    mol.add_bond(9, 0, 10, BondOrder::Single);
    mol.add_bond(10, 0, 11, BondOrder::Single);
    mol.add_bond(11, 1, 12, BondOrder::Single);
    mol.add_bond(12, 5, 13, BondOrder::Single);
    mol.add_bond(13, 6, 14, BondOrder::Single);
    mol.add_bond(14, 6, 15, BondOrder::Single);
    mol.add_bond(15, 7, 16, BondOrder::Single);
    mol.add_bond(16, 7, 17, BondOrder::Single);
    mol.add_bond(17, 7, 18, BondOrder::Single);
    mol.add_bond(18, 8, 19, BondOrder::Single);
    mol.add_bond(19, 8, 20, BondOrder::Single);
    mol.add_bond(20, 8, 21, BondOrder::Single);

    mol
}

pub fn build_proline_zwitterion() -> TestMolecule {
    let mut mol = TestMolecule::new();

    mol.add_atom(0, Element::N, 1);
    mol.add_atom(1, Element::C, 0);
    mol.add_atom(2, Element::C, 0);
    mol.add_atom(3, Element::C, 0);
    mol.add_atom(4, Element::C, 0);
    mol.add_atom(5, Element::C, 0);
    mol.add_atom(6, Element::O, -1);
    mol.add_atom(7, Element::O, 0);
    mol.add_atom(8, Element::H, 0);
    mol.add_atom(9, Element::H, 0);
    mol.add_atom(10, Element::H, 0);
    mol.add_atom(11, Element::H, 0);
    mol.add_atom(12, Element::H, 0);
    mol.add_atom(13, Element::H, 0);
    mol.add_atom(14, Element::H, 0);
    mol.add_atom(15, Element::H, 0);
    mol.add_atom(16, Element::H, 0);

    mol.add_bond(0, 0, 1, BondOrder::Single);
    mol.add_bond(1, 1, 2, BondOrder::Single);
    mol.add_bond(2, 2, 3, BondOrder::Single);
    mol.add_bond(3, 3, 4, BondOrder::Single);
    mol.add_bond(4, 4, 0, BondOrder::Single);
    mol.add_bond(5, 1, 5, BondOrder::Single);
    mol.add_bond(6, 5, 6, BondOrder::Single);
    mol.add_bond(7, 5, 7, BondOrder::Double);
    mol.add_bond(8, 0, 8, BondOrder::Single);
    mol.add_bond(9, 0, 9, BondOrder::Single);
    mol.add_bond(10, 1, 10, BondOrder::Single);
    mol.add_bond(11, 2, 11, BondOrder::Single);
    mol.add_bond(12, 2, 12, BondOrder::Single);
    mol.add_bond(13, 3, 13, BondOrder::Single);
    mol.add_bond(14, 3, 14, BondOrder::Single);
    mol.add_bond(15, 4, 15, BondOrder::Single);
    mol.add_bond(16, 4, 16, BondOrder::Single);

    mol
}

pub fn build_serine_zwitterion() -> TestMolecule {
    let mut mol = TestMolecule::new();

    mol.add_atom(0, Element::N, 1);
    mol.add_atom(1, Element::C, 0);
    mol.add_atom(2, Element::C, 0);
    mol.add_atom(3, Element::O, -1);
    mol.add_atom(4, Element::O, 0);
    mol.add_atom(5, Element::C, 0);
    mol.add_atom(6, Element::O, 0);
    mol.add_atom(7, Element::H, 0);
    mol.add_atom(8, Element::H, 0);
    mol.add_atom(9, Element::H, 0);
    mol.add_atom(10, Element::H, 0);
    mol.add_atom(11, Element::H, 0);
    mol.add_atom(12, Element::H, 0);
    mol.add_atom(13, Element::H, 0);

    mol.add_bond(0, 0, 1, BondOrder::Single);
    mol.add_bond(1, 1, 2, BondOrder::Single);
    mol.add_bond(2, 2, 3, BondOrder::Single);
    mol.add_bond(3, 2, 4, BondOrder::Double);
    mol.add_bond(4, 1, 5, BondOrder::Single);
    mol.add_bond(5, 5, 6, BondOrder::Single);
    mol.add_bond(6, 0, 7, BondOrder::Single);
    mol.add_bond(7, 0, 8, BondOrder::Single);
    mol.add_bond(8, 0, 9, BondOrder::Single);
    mol.add_bond(9, 1, 10, BondOrder::Single);
    mol.add_bond(10, 5, 11, BondOrder::Single);
    mol.add_bond(11, 5, 12, BondOrder::Single);
    mol.add_bond(12, 6, 13, BondOrder::Single);

    mol
}

pub fn build_threonine_zwitterion() -> TestMolecule {
    let mut mol = TestMolecule::new();

    mol.add_atom(0, Element::N, 1);
    mol.add_atom(1, Element::C, 0);
    mol.add_atom(2, Element::C, 0);
    mol.add_atom(3, Element::O, -1);
    mol.add_atom(4, Element::O, 0);
    mol.add_atom(5, Element::C, 0);
    mol.add_atom(6, Element::O, 0);
    mol.add_atom(7, Element::C, 0);
    mol.add_atom(8, Element::H, 0);
    mol.add_atom(9, Element::H, 0);
    mol.add_atom(10, Element::H, 0);
    mol.add_atom(11, Element::H, 0);
    mol.add_atom(12, Element::H, 0);
    mol.add_atom(13, Element::H, 0);
    mol.add_atom(14, Element::H, 0);
    mol.add_atom(15, Element::H, 0);
    mol.add_atom(16, Element::H, 0);

    mol.add_bond(0, 0, 1, BondOrder::Single);
    mol.add_bond(1, 1, 2, BondOrder::Single);
    mol.add_bond(2, 2, 3, BondOrder::Single);
    mol.add_bond(3, 2, 4, BondOrder::Double);
    mol.add_bond(4, 1, 5, BondOrder::Single);
    mol.add_bond(5, 5, 6, BondOrder::Single);
    mol.add_bond(6, 5, 7, BondOrder::Single);
    mol.add_bond(7, 0, 8, BondOrder::Single);
    mol.add_bond(8, 0, 9, BondOrder::Single);
    mol.add_bond(9, 0, 10, BondOrder::Single);
    mol.add_bond(10, 1, 11, BondOrder::Single);
    mol.add_bond(11, 5, 12, BondOrder::Single);
    mol.add_bond(12, 6, 13, BondOrder::Single);
    mol.add_bond(13, 7, 14, BondOrder::Single);
    mol.add_bond(14, 7, 15, BondOrder::Single);
    mol.add_bond(15, 7, 16, BondOrder::Single);

    mol
}

pub fn build_cysteine_zwitterion() -> TestMolecule {
    let mut mol = TestMolecule::new();

    mol.add_atom(0, Element::N, 1);
    mol.add_atom(1, Element::C, 0);
    mol.add_atom(2, Element::C, 0);
    mol.add_atom(3, Element::O, -1);
    mol.add_atom(4, Element::O, 0);
    mol.add_atom(5, Element::C, 0);
    mol.add_atom(6, Element::S, 0);
    mol.add_atom(7, Element::H, 0);
    mol.add_atom(8, Element::H, 0);
    mol.add_atom(9, Element::H, 0);
    mol.add_atom(10, Element::H, 0);
    mol.add_atom(11, Element::H, 0);
    mol.add_atom(12, Element::H, 0);
    mol.add_atom(13, Element::H, 0);

    mol.add_bond(0, 0, 1, BondOrder::Single);
    mol.add_bond(1, 1, 2, BondOrder::Single);
    mol.add_bond(2, 2, 3, BondOrder::Single);
    mol.add_bond(3, 2, 4, BondOrder::Double);
    mol.add_bond(4, 1, 5, BondOrder::Single);
    mol.add_bond(5, 5, 6, BondOrder::Single);
    mol.add_bond(6, 0, 7, BondOrder::Single);
    mol.add_bond(7, 0, 8, BondOrder::Single);
    mol.add_bond(8, 0, 9, BondOrder::Single);
    mol.add_bond(9, 1, 10, BondOrder::Single);
    mol.add_bond(10, 5, 11, BondOrder::Single);
    mol.add_bond(11, 5, 12, BondOrder::Single);
    mol.add_bond(12, 6, 13, BondOrder::Single);

    mol
}

pub fn build_methionine_zwitterion() -> TestMolecule {
    let mut mol = TestMolecule::new();

    mol.add_atom(0, Element::N, 1);
    mol.add_atom(1, Element::C, 0);
    mol.add_atom(2, Element::C, 0);
    mol.add_atom(3, Element::O, -1);
    mol.add_atom(4, Element::O, 0);
    mol.add_atom(5, Element::C, 0);
    mol.add_atom(6, Element::C, 0);
    mol.add_atom(7, Element::S, 0);
    mol.add_atom(8, Element::C, 0);
    mol.add_atom(9, Element::H, 0);
    mol.add_atom(10, Element::H, 0);
    mol.add_atom(11, Element::H, 0);
    mol.add_atom(12, Element::H, 0);
    mol.add_atom(13, Element::H, 0);
    mol.add_atom(14, Element::H, 0);
    mol.add_atom(15, Element::H, 0);
    mol.add_atom(16, Element::H, 0);
    mol.add_atom(17, Element::H, 0);
    mol.add_atom(18, Element::H, 0);
    mol.add_atom(19, Element::H, 0);

    mol.add_bond(0, 0, 1, BondOrder::Single);
    mol.add_bond(1, 1, 2, BondOrder::Single);
    mol.add_bond(2, 2, 3, BondOrder::Single);
    mol.add_bond(3, 2, 4, BondOrder::Double);
    mol.add_bond(4, 1, 5, BondOrder::Single);
    mol.add_bond(5, 5, 6, BondOrder::Single);
    mol.add_bond(6, 6, 7, BondOrder::Single);
    mol.add_bond(7, 7, 8, BondOrder::Single);
    mol.add_bond(8, 0, 9, BondOrder::Single);
    mol.add_bond(9, 0, 10, BondOrder::Single);
    mol.add_bond(10, 0, 11, BondOrder::Single);
    mol.add_bond(11, 1, 12, BondOrder::Single);
    mol.add_bond(12, 5, 13, BondOrder::Single);
    mol.add_bond(13, 5, 14, BondOrder::Single);
    mol.add_bond(14, 6, 15, BondOrder::Single);
    mol.add_bond(15, 6, 16, BondOrder::Single);
    mol.add_bond(16, 8, 17, BondOrder::Single);
    mol.add_bond(17, 8, 18, BondOrder::Single);
    mol.add_bond(18, 8, 19, BondOrder::Single);

    mol
}

pub fn build_aspartate_zwitterion() -> TestMolecule {
    let mut mol = TestMolecule::new();

    mol.add_atom(0, Element::N, 1);
    mol.add_atom(1, Element::C, 0);
    mol.add_atom(2, Element::C, 0);
    mol.add_atom(3, Element::O, -1);
    mol.add_atom(4, Element::O, 0);
    mol.add_atom(5, Element::C, 0);
    mol.add_atom(6, Element::C, 0);
    mol.add_atom(7, Element::O, 0);
    mol.add_atom(8, Element::O, -1);
    mol.add_atom(9, Element::H, 0);
    mol.add_atom(10, Element::H, 0);
    mol.add_atom(11, Element::H, 0);
    mol.add_atom(12, Element::H, 0);
    mol.add_atom(13, Element::H, 0);
    mol.add_atom(14, Element::H, 0);

    mol.add_bond(0, 0, 1, BondOrder::Single);
    mol.add_bond(1, 1, 2, BondOrder::Single);
    mol.add_bond(2, 2, 3, BondOrder::Single);
    mol.add_bond(3, 2, 4, BondOrder::Double);
    mol.add_bond(4, 1, 5, BondOrder::Single);
    mol.add_bond(5, 5, 6, BondOrder::Single);
    mol.add_bond(6, 6, 7, BondOrder::Double);
    mol.add_bond(7, 6, 8, BondOrder::Single);
    mol.add_bond(8, 0, 9, BondOrder::Single);
    mol.add_bond(9, 0, 10, BondOrder::Single);
    mol.add_bond(10, 0, 11, BondOrder::Single);
    mol.add_bond(11, 1, 12, BondOrder::Single);
    mol.add_bond(12, 5, 13, BondOrder::Single);
    mol.add_bond(13, 5, 14, BondOrder::Single);

    mol
}

pub fn build_asparagine_zwitterion() -> TestMolecule {
    let mut mol = TestMolecule::new();

    mol.add_atom(0, Element::N, 1);
    mol.add_atom(1, Element::C, 0);
    mol.add_atom(2, Element::C, 0);
    mol.add_atom(3, Element::O, -1);
    mol.add_atom(4, Element::O, 0);
    mol.add_atom(5, Element::C, 0);
    mol.add_atom(6, Element::C, 0);
    mol.add_atom(7, Element::O, 0);
    mol.add_atom(8, Element::N, 0);
    mol.add_atom(9, Element::H, 0);
    mol.add_atom(10, Element::H, 0);
    mol.add_atom(11, Element::H, 0);
    mol.add_atom(12, Element::H, 0);
    mol.add_atom(13, Element::H, 0);
    mol.add_atom(14, Element::H, 0);
    mol.add_atom(15, Element::H, 0);
    mol.add_atom(16, Element::H, 0);

    mol.add_bond(0, 0, 1, BondOrder::Single);
    mol.add_bond(1, 1, 2, BondOrder::Single);
    mol.add_bond(2, 2, 3, BondOrder::Single);
    mol.add_bond(3, 2, 4, BondOrder::Double);
    mol.add_bond(4, 1, 5, BondOrder::Single);
    mol.add_bond(5, 5, 6, BondOrder::Single);
    mol.add_bond(6, 6, 7, BondOrder::Double);
    mol.add_bond(7, 6, 8, BondOrder::Single);
    mol.add_bond(8, 0, 9, BondOrder::Single);
    mol.add_bond(9, 0, 10, BondOrder::Single);
    mol.add_bond(10, 0, 11, BondOrder::Single);
    mol.add_bond(11, 1, 12, BondOrder::Single);
    mol.add_bond(12, 5, 13, BondOrder::Single);
    mol.add_bond(13, 5, 14, BondOrder::Single);
    mol.add_bond(14, 8, 15, BondOrder::Single);
    mol.add_bond(15, 8, 16, BondOrder::Single);

    mol
}

pub fn build_glutamate_zwitterion() -> TestMolecule {
    let mut mol = TestMolecule::new();

    mol.add_atom(0, Element::N, 1);
    mol.add_atom(1, Element::C, 0);
    mol.add_atom(2, Element::C, 0);
    mol.add_atom(3, Element::O, -1);
    mol.add_atom(4, Element::O, 0);
    mol.add_atom(5, Element::C, 0);
    mol.add_atom(6, Element::C, 0);
    mol.add_atom(7, Element::C, 0);
    mol.add_atom(8, Element::O, 0);
    mol.add_atom(9, Element::O, -1);
    mol.add_atom(10, Element::H, 0);
    mol.add_atom(11, Element::H, 0);
    mol.add_atom(12, Element::H, 0);
    mol.add_atom(13, Element::H, 0);
    mol.add_atom(14, Element::H, 0);
    mol.add_atom(15, Element::H, 0);
    mol.add_atom(16, Element::H, 0);
    mol.add_atom(17, Element::H, 0);

    mol.add_bond(0, 0, 1, BondOrder::Single);
    mol.add_bond(1, 1, 2, BondOrder::Single);
    mol.add_bond(2, 2, 3, BondOrder::Single);
    mol.add_bond(3, 2, 4, BondOrder::Double);
    mol.add_bond(4, 1, 5, BondOrder::Single);
    mol.add_bond(5, 5, 6, BondOrder::Single);
    mol.add_bond(6, 6, 7, BondOrder::Single);
    mol.add_bond(7, 7, 8, BondOrder::Double);
    mol.add_bond(8, 7, 9, BondOrder::Single);
    mol.add_bond(9, 0, 10, BondOrder::Single);
    mol.add_bond(10, 0, 11, BondOrder::Single);
    mol.add_bond(11, 0, 12, BondOrder::Single);
    mol.add_bond(12, 1, 13, BondOrder::Single);
    mol.add_bond(13, 5, 14, BondOrder::Single);
    mol.add_bond(14, 5, 15, BondOrder::Single);
    mol.add_bond(15, 6, 16, BondOrder::Single);
    mol.add_bond(16, 6, 17, BondOrder::Single);

    mol
}

pub fn build_glutamine_zwitterion() -> TestMolecule {
    let mut mol = TestMolecule::new();

    mol.add_atom(0, Element::N, 1);
    mol.add_atom(1, Element::C, 0);
    mol.add_atom(2, Element::C, 0);
    mol.add_atom(3, Element::O, -1);
    mol.add_atom(4, Element::O, 0);
    mol.add_atom(5, Element::C, 0);
    mol.add_atom(6, Element::C, 0);
    mol.add_atom(7, Element::C, 0);
    mol.add_atom(8, Element::O, 0);
    mol.add_atom(9, Element::N, 0);
    mol.add_atom(10, Element::H, 0);
    mol.add_atom(11, Element::H, 0);
    mol.add_atom(12, Element::H, 0);
    mol.add_atom(13, Element::H, 0);
    mol.add_atom(14, Element::H, 0);
    mol.add_atom(15, Element::H, 0);
    mol.add_atom(16, Element::H, 0);
    mol.add_atom(17, Element::H, 0);
    mol.add_atom(18, Element::H, 0);
    mol.add_atom(19, Element::H, 0);

    mol.add_bond(0, 0, 1, BondOrder::Single);
    mol.add_bond(1, 1, 2, BondOrder::Single);
    mol.add_bond(2, 2, 3, BondOrder::Single);
    mol.add_bond(3, 2, 4, BondOrder::Double);
    mol.add_bond(4, 1, 5, BondOrder::Single);
    mol.add_bond(5, 5, 6, BondOrder::Single);
    mol.add_bond(6, 6, 7, BondOrder::Single);
    mol.add_bond(7, 7, 8, BondOrder::Double);
    mol.add_bond(8, 7, 9, BondOrder::Single);
    mol.add_bond(9, 0, 10, BondOrder::Single);
    mol.add_bond(10, 0, 11, BondOrder::Single);
    mol.add_bond(11, 0, 12, BondOrder::Single);
    mol.add_bond(12, 1, 13, BondOrder::Single);
    mol.add_bond(13, 5, 14, BondOrder::Single);
    mol.add_bond(14, 5, 15, BondOrder::Single);
    mol.add_bond(15, 6, 16, BondOrder::Single);
    mol.add_bond(16, 6, 17, BondOrder::Single);
    mol.add_bond(17, 9, 18, BondOrder::Single);
    mol.add_bond(18, 9, 19, BondOrder::Single);

    mol
}

pub fn build_lysine_zwitterion() -> TestMolecule {
    let mut mol = TestMolecule::new();

    mol.add_atom(0, Element::N, 1);
    mol.add_atom(1, Element::C, 0);
    mol.add_atom(2, Element::C, 0);
    mol.add_atom(3, Element::O, -1);
    mol.add_atom(4, Element::O, 0);
    mol.add_atom(5, Element::C, 0);
    mol.add_atom(6, Element::C, 0);
    mol.add_atom(7, Element::C, 0);
    mol.add_atom(8, Element::C, 0);
    mol.add_atom(9, Element::N, 1);
    mol.add_atom(10, Element::H, 0);
    mol.add_atom(11, Element::H, 0);
    mol.add_atom(12, Element::H, 0);
    mol.add_atom(13, Element::H, 0);
    mol.add_atom(14, Element::H, 0);
    mol.add_atom(15, Element::H, 0);
    mol.add_atom(16, Element::H, 0);
    mol.add_atom(17, Element::H, 0);
    mol.add_atom(18, Element::H, 0);
    mol.add_atom(19, Element::H, 0);
    mol.add_atom(20, Element::H, 0);
    mol.add_atom(21, Element::H, 0);
    mol.add_atom(22, Element::H, 0);
    mol.add_atom(23, Element::H, 0);
    mol.add_atom(24, Element::H, 0);

    mol.add_bond(0, 0, 1, BondOrder::Single);
    mol.add_bond(1, 1, 2, BondOrder::Single);
    mol.add_bond(2, 2, 3, BondOrder::Single);
    mol.add_bond(3, 2, 4, BondOrder::Double);
    mol.add_bond(4, 1, 5, BondOrder::Single);
    mol.add_bond(5, 5, 6, BondOrder::Single);
    mol.add_bond(6, 6, 7, BondOrder::Single);
    mol.add_bond(7, 7, 8, BondOrder::Single);
    mol.add_bond(8, 8, 9, BondOrder::Single);
    mol.add_bond(9, 0, 10, BondOrder::Single);
    mol.add_bond(10, 0, 11, BondOrder::Single);
    mol.add_bond(11, 0, 12, BondOrder::Single);
    mol.add_bond(12, 1, 13, BondOrder::Single);
    mol.add_bond(13, 5, 14, BondOrder::Single);
    mol.add_bond(14, 5, 15, BondOrder::Single);
    mol.add_bond(15, 6, 16, BondOrder::Single);
    mol.add_bond(16, 6, 17, BondOrder::Single);
    mol.add_bond(17, 7, 18, BondOrder::Single);
    mol.add_bond(18, 7, 19, BondOrder::Single);
    mol.add_bond(19, 8, 20, BondOrder::Single);
    mol.add_bond(20, 8, 21, BondOrder::Single);
    mol.add_bond(21, 9, 22, BondOrder::Single);
    mol.add_bond(22, 9, 23, BondOrder::Single);
    mol.add_bond(23, 9, 24, BondOrder::Single);

    mol
}

pub fn build_phenylalanine_zwitterion_aromatic() -> TestMolecule {
    let mut mol = TestMolecule::new();

    mol.add_atom(0, Element::N, 1);
    mol.add_atom(1, Element::C, 0);
    mol.add_atom(2, Element::C, 0);
    mol.add_atom(3, Element::O, -1);
    mol.add_atom(4, Element::O, 0);
    mol.add_atom(5, Element::C, 0);
    mol.add_atom(6, Element::C, 0);
    mol.add_atom(7, Element::C, 0);
    mol.add_atom(8, Element::C, 0);
    mol.add_atom(9, Element::C, 0);
    mol.add_atom(10, Element::C, 0);
    mol.add_atom(11, Element::C, 0);
    mol.add_atom(12, Element::H, 0);
    mol.add_atom(13, Element::H, 0);
    mol.add_atom(14, Element::H, 0);
    mol.add_atom(15, Element::H, 0);
    mol.add_atom(16, Element::H, 0);
    mol.add_atom(17, Element::H, 0);
    mol.add_atom(18, Element::H, 0);
    mol.add_atom(19, Element::H, 0);
    mol.add_atom(20, Element::H, 0);
    mol.add_atom(21, Element::H, 0);
    mol.add_atom(22, Element::H, 0);

    mol.add_bond(0, 0, 1, BondOrder::Single);
    mol.add_bond(1, 1, 2, BondOrder::Single);
    mol.add_bond(2, 2, 3, BondOrder::Single);
    mol.add_bond(3, 2, 4, BondOrder::Double);
    mol.add_bond(4, 1, 5, BondOrder::Single);
    mol.add_bond(5, 5, 6, BondOrder::Single);
    mol.add_bond(6, 6, 7, BondOrder::Aromatic);
    mol.add_bond(7, 7, 9, BondOrder::Aromatic);
    mol.add_bond(8, 9, 11, BondOrder::Aromatic);
    mol.add_bond(9, 11, 10, BondOrder::Aromatic);
    mol.add_bond(10, 10, 8, BondOrder::Aromatic);
    mol.add_bond(11, 8, 6, BondOrder::Aromatic);
    mol.add_bond(12, 0, 12, BondOrder::Single);
    mol.add_bond(13, 0, 13, BondOrder::Single);
    mol.add_bond(14, 0, 14, BondOrder::Single);
    mol.add_bond(15, 1, 15, BondOrder::Single);
    mol.add_bond(16, 5, 16, BondOrder::Single);
    mol.add_bond(17, 5, 17, BondOrder::Single);
    mol.add_bond(18, 7, 18, BondOrder::Single);
    mol.add_bond(19, 8, 19, BondOrder::Single);
    mol.add_bond(20, 9, 20, BondOrder::Single);
    mol.add_bond(21, 10, 21, BondOrder::Single);
    mol.add_bond(22, 11, 22, BondOrder::Single);

    mol
}

pub fn build_phenylalanine_zwitterion_kekule() -> TestMolecule {
    let mut mol = TestMolecule::new();

    mol.add_atom(0, Element::N, 1);
    mol.add_atom(1, Element::C, 0);
    mol.add_atom(2, Element::C, 0);
    mol.add_atom(3, Element::O, -1);
    mol.add_atom(4, Element::O, 0);
    mol.add_atom(5, Element::C, 0);
    mol.add_atom(6, Element::C, 0);
    mol.add_atom(7, Element::C, 0);
    mol.add_atom(8, Element::C, 0);
    mol.add_atom(9, Element::C, 0);
    mol.add_atom(10, Element::C, 0);
    mol.add_atom(11, Element::C, 0);
    mol.add_atom(12, Element::H, 0);
    mol.add_atom(13, Element::H, 0);
    mol.add_atom(14, Element::H, 0);
    mol.add_atom(15, Element::H, 0);
    mol.add_atom(16, Element::H, 0);
    mol.add_atom(17, Element::H, 0);
    mol.add_atom(18, Element::H, 0);
    mol.add_atom(19, Element::H, 0);
    mol.add_atom(20, Element::H, 0);
    mol.add_atom(21, Element::H, 0);
    mol.add_atom(22, Element::H, 0);

    mol.add_bond(0, 0, 1, BondOrder::Single);
    mol.add_bond(1, 1, 2, BondOrder::Single);
    mol.add_bond(2, 2, 3, BondOrder::Single);
    mol.add_bond(3, 2, 4, BondOrder::Double);
    mol.add_bond(4, 1, 5, BondOrder::Single);
    mol.add_bond(5, 5, 6, BondOrder::Single);
    mol.add_bond(6, 6, 7, BondOrder::Double);
    mol.add_bond(7, 7, 9, BondOrder::Single);
    mol.add_bond(8, 9, 11, BondOrder::Double);
    mol.add_bond(9, 11, 10, BondOrder::Single);
    mol.add_bond(10, 10, 8, BondOrder::Double);
    mol.add_bond(11, 8, 6, BondOrder::Single);
    mol.add_bond(12, 0, 12, BondOrder::Single);
    mol.add_bond(13, 0, 13, BondOrder::Single);
    mol.add_bond(14, 0, 14, BondOrder::Single);
    mol.add_bond(15, 1, 15, BondOrder::Single);
    mol.add_bond(16, 5, 16, BondOrder::Single);
    mol.add_bond(17, 5, 17, BondOrder::Single);
    mol.add_bond(18, 7, 18, BondOrder::Single);
    mol.add_bond(19, 8, 19, BondOrder::Single);
    mol.add_bond(20, 9, 20, BondOrder::Single);
    mol.add_bond(21, 10, 21, BondOrder::Single);
    mol.add_bond(22, 11, 22, BondOrder::Single);

    mol
}

pub fn build_tyrosine_zwitterion_aromatic() -> TestMolecule {
    let mut mol = TestMolecule::new();

    mol.add_atom(0, Element::N, 1);
    mol.add_atom(1, Element::C, 0);
    mol.add_atom(2, Element::C, 0);
    mol.add_atom(3, Element::O, -1);
    mol.add_atom(4, Element::O, 0);
    mol.add_atom(5, Element::C, 0);
    mol.add_atom(6, Element::C, 0);
    mol.add_atom(7, Element::C, 0);
    mol.add_atom(8, Element::C, 0);
    mol.add_atom(9, Element::C, 0);
    mol.add_atom(10, Element::C, 0);
    mol.add_atom(11, Element::C, 0);
    mol.add_atom(12, Element::O, 0);
    mol.add_atom(13, Element::H, 0);
    mol.add_atom(14, Element::H, 0);
    mol.add_atom(15, Element::H, 0);
    mol.add_atom(16, Element::H, 0);
    mol.add_atom(17, Element::H, 0);
    mol.add_atom(18, Element::H, 0);
    mol.add_atom(19, Element::H, 0);
    mol.add_atom(20, Element::H, 0);
    mol.add_atom(21, Element::H, 0);
    mol.add_atom(22, Element::H, 0);
    mol.add_atom(23, Element::H, 0);

    mol.add_bond(0, 0, 1, BondOrder::Single);
    mol.add_bond(1, 1, 2, BondOrder::Single);
    mol.add_bond(2, 2, 3, BondOrder::Single);
    mol.add_bond(3, 2, 4, BondOrder::Double);
    mol.add_bond(4, 1, 5, BondOrder::Single);
    mol.add_bond(5, 5, 6, BondOrder::Single);
    mol.add_bond(6, 6, 7, BondOrder::Aromatic);
    mol.add_bond(7, 7, 9, BondOrder::Aromatic);
    mol.add_bond(8, 9, 11, BondOrder::Aromatic);
    mol.add_bond(9, 11, 10, BondOrder::Aromatic);
    mol.add_bond(10, 10, 8, BondOrder::Aromatic);
    mol.add_bond(11, 8, 6, BondOrder::Aromatic);
    mol.add_bond(12, 11, 12, BondOrder::Single);
    mol.add_bond(13, 0, 13, BondOrder::Single);
    mol.add_bond(14, 0, 14, BondOrder::Single);
    mol.add_bond(15, 0, 15, BondOrder::Single);
    mol.add_bond(16, 1, 16, BondOrder::Single);
    mol.add_bond(17, 5, 17, BondOrder::Single);
    mol.add_bond(18, 5, 18, BondOrder::Single);
    mol.add_bond(19, 7, 19, BondOrder::Single);
    mol.add_bond(20, 8, 20, BondOrder::Single);
    mol.add_bond(21, 9, 21, BondOrder::Single);
    mol.add_bond(22, 10, 22, BondOrder::Single);
    mol.add_bond(23, 12, 23, BondOrder::Single);

    mol
}

pub fn build_tyrosine_zwitterion_kekule() -> TestMolecule {
    let mut mol = TestMolecule::new();

    mol.add_atom(0, Element::N, 1);
    mol.add_atom(1, Element::C, 0);
    mol.add_atom(2, Element::C, 0);
    mol.add_atom(3, Element::O, -1);
    mol.add_atom(4, Element::O, 0);
    mol.add_atom(5, Element::C, 0);
    mol.add_atom(6, Element::C, 0);
    mol.add_atom(7, Element::C, 0);
    mol.add_atom(8, Element::C, 0);
    mol.add_atom(9, Element::C, 0);
    mol.add_atom(10, Element::C, 0);
    mol.add_atom(11, Element::C, 0);
    mol.add_atom(12, Element::O, 0);
    mol.add_atom(13, Element::H, 0);
    mol.add_atom(14, Element::H, 0);
    mol.add_atom(15, Element::H, 0);
    mol.add_atom(16, Element::H, 0);
    mol.add_atom(17, Element::H, 0);
    mol.add_atom(18, Element::H, 0);
    mol.add_atom(19, Element::H, 0);
    mol.add_atom(20, Element::H, 0);
    mol.add_atom(21, Element::H, 0);
    mol.add_atom(22, Element::H, 0);
    mol.add_atom(23, Element::H, 0);

    mol.add_bond(0, 0, 1, BondOrder::Single);
    mol.add_bond(1, 1, 2, BondOrder::Single);
    mol.add_bond(2, 2, 3, BondOrder::Single);
    mol.add_bond(3, 2, 4, BondOrder::Double);
    mol.add_bond(4, 1, 5, BondOrder::Single);
    mol.add_bond(5, 5, 6, BondOrder::Single);
    mol.add_bond(6, 6, 7, BondOrder::Double);
    mol.add_bond(7, 7, 9, BondOrder::Single);
    mol.add_bond(8, 9, 11, BondOrder::Double);
    mol.add_bond(9, 11, 10, BondOrder::Single);
    mol.add_bond(10, 10, 8, BondOrder::Double);
    mol.add_bond(11, 8, 6, BondOrder::Single);
    mol.add_bond(12, 11, 12, BondOrder::Single);
    mol.add_bond(13, 0, 13, BondOrder::Single);
    mol.add_bond(14, 0, 14, BondOrder::Single);
    mol.add_bond(15, 0, 15, BondOrder::Single);
    mol.add_bond(16, 1, 16, BondOrder::Single);
    mol.add_bond(17, 5, 17, BondOrder::Single);
    mol.add_bond(18, 5, 18, BondOrder::Single);
    mol.add_bond(19, 7, 19, BondOrder::Single);
    mol.add_bond(20, 8, 20, BondOrder::Single);
    mol.add_bond(21, 9, 21, BondOrder::Single);
    mol.add_bond(22, 10, 22, BondOrder::Single);
    mol.add_bond(23, 12, 23, BondOrder::Single);

    mol
}

pub fn build_histidine_zwitterion_aromatic() -> TestMolecule {
    let mut mol = TestMolecule::new();

    mol.add_atom(0, Element::N, 1);
    mol.add_atom(1, Element::C, 0);
    mol.add_atom(2, Element::C, 0);
    mol.add_atom(3, Element::O, -1);
    mol.add_atom(4, Element::O, 0);
    mol.add_atom(5, Element::C, 0);
    mol.add_atom(6, Element::C, 0);
    mol.add_atom(7, Element::N, 0);
    mol.add_atom(8, Element::C, 0);
    mol.add_atom(9, Element::N, 0);
    mol.add_atom(10, Element::C, 0);

    for id in 11..=19 {
        mol.add_atom(id, Element::H, 0);
    }

    mol.add_bond(0, 0, 1, BondOrder::Single);
    mol.add_bond(1, 1, 2, BondOrder::Single);
    mol.add_bond(2, 2, 3, BondOrder::Single);
    mol.add_bond(3, 2, 4, BondOrder::Double);
    mol.add_bond(4, 1, 5, BondOrder::Single);
    mol.add_bond(5, 5, 6, BondOrder::Single);

    mol.add_bond(6, 6, 7, BondOrder::Aromatic);
    mol.add_bond(7, 7, 8, BondOrder::Aromatic);
    mol.add_bond(8, 8, 9, BondOrder::Aromatic);
    mol.add_bond(9, 9, 10, BondOrder::Aromatic);
    mol.add_bond(10, 10, 6, BondOrder::Aromatic);

    mol.add_bond(11, 0, 11, BondOrder::Single);
    mol.add_bond(12, 0, 12, BondOrder::Single);
    mol.add_bond(13, 0, 13, BondOrder::Single);
    mol.add_bond(14, 1, 14, BondOrder::Single);
    mol.add_bond(15, 5, 15, BondOrder::Single);
    mol.add_bond(16, 5, 16, BondOrder::Single);
    mol.add_bond(17, 7, 17, BondOrder::Single);
    mol.add_bond(18, 8, 18, BondOrder::Single);
    mol.add_bond(19, 10, 19, BondOrder::Single);

    mol
}

pub fn build_histidine_zwitterion_kekule() -> TestMolecule {
    let mut mol = TestMolecule::new();

    mol.add_atom(0, Element::N, 1);
    mol.add_atom(1, Element::C, 0);
    mol.add_atom(2, Element::C, 0);
    mol.add_atom(3, Element::O, -1);
    mol.add_atom(4, Element::O, 0);
    mol.add_atom(5, Element::C, 0);
    mol.add_atom(6, Element::C, 0);
    mol.add_atom(7, Element::N, 0);
    mol.add_atom(8, Element::C, 0);
    mol.add_atom(9, Element::N, 0);
    mol.add_atom(10, Element::C, 0);

    for id in 11..=19 {
        mol.add_atom(id, Element::H, 0);
    }

    mol.add_bond(0, 0, 1, BondOrder::Single);
    mol.add_bond(1, 1, 2, BondOrder::Single);
    mol.add_bond(2, 2, 3, BondOrder::Single);
    mol.add_bond(3, 2, 4, BondOrder::Double);
    mol.add_bond(4, 1, 5, BondOrder::Single);
    mol.add_bond(5, 5, 6, BondOrder::Single);

    mol.add_bond(6, 6, 7, BondOrder::Single);
    mol.add_bond(7, 7, 8, BondOrder::Double);
    mol.add_bond(8, 8, 9, BondOrder::Single);
    mol.add_bond(9, 9, 10, BondOrder::Double);
    mol.add_bond(10, 10, 6, BondOrder::Single);

    mol.add_bond(11, 0, 11, BondOrder::Single);
    mol.add_bond(12, 0, 12, BondOrder::Single);
    mol.add_bond(13, 0, 13, BondOrder::Single);
    mol.add_bond(14, 1, 14, BondOrder::Single);
    mol.add_bond(15, 5, 15, BondOrder::Single);
    mol.add_bond(16, 5, 16, BondOrder::Single);
    mol.add_bond(17, 7, 17, BondOrder::Single);
    mol.add_bond(18, 8, 18, BondOrder::Single);
    mol.add_bond(19, 10, 19, BondOrder::Single);

    mol
}

pub fn build_tryptophan_zwitterion_aromatic() -> TestMolecule {
    let mut mol = TestMolecule::new();

    mol.add_atom(0, Element::N, 1);
    mol.add_atom(1, Element::C, 0);
    mol.add_atom(2, Element::C, 0);
    mol.add_atom(3, Element::O, -1);
    mol.add_atom(4, Element::O, 0);
    mol.add_atom(5, Element::C, 0);
    mol.add_atom(6, Element::C, 0);
    mol.add_atom(7, Element::C, 0);
    mol.add_atom(8, Element::N, 0);
    mol.add_atom(9, Element::C, 0);
    mol.add_atom(10, Element::C, 0);
    mol.add_atom(11, Element::C, 0);
    mol.add_atom(12, Element::C, 0);
    mol.add_atom(13, Element::C, 0);
    mol.add_atom(14, Element::C, 0);
    mol.add_atom(15, Element::H, 0);
    mol.add_atom(16, Element::H, 0);
    mol.add_atom(17, Element::H, 0);
    mol.add_atom(18, Element::H, 0);
    mol.add_atom(19, Element::H, 0);
    mol.add_atom(20, Element::H, 0);
    mol.add_atom(21, Element::H, 0);
    mol.add_atom(22, Element::H, 0);
    mol.add_atom(23, Element::H, 0);
    mol.add_atom(24, Element::H, 0);
    mol.add_atom(25, Element::H, 0);
    mol.add_atom(26, Element::H, 0);

    mol.add_bond(0, 0, 1, BondOrder::Single);
    mol.add_bond(1, 1, 2, BondOrder::Single);
    mol.add_bond(2, 2, 3, BondOrder::Single);
    mol.add_bond(3, 2, 4, BondOrder::Double);
    mol.add_bond(4, 1, 5, BondOrder::Single);
    mol.add_bond(5, 5, 6, BondOrder::Single);
    mol.add_bond(6, 6, 7, BondOrder::Aromatic);
    mol.add_bond(7, 7, 8, BondOrder::Aromatic);
    mol.add_bond(8, 8, 9, BondOrder::Aromatic);
    mol.add_bond(9, 9, 10, BondOrder::Aromatic);
    mol.add_bond(10, 10, 6, BondOrder::Aromatic);
    mol.add_bond(11, 10, 11, BondOrder::Aromatic);
    mol.add_bond(12, 11, 12, BondOrder::Aromatic);
    mol.add_bond(13, 12, 13, BondOrder::Aromatic);
    mol.add_bond(14, 13, 14, BondOrder::Aromatic);
    mol.add_bond(15, 14, 9, BondOrder::Aromatic);
    mol.add_bond(16, 0, 15, BondOrder::Single);
    mol.add_bond(17, 0, 16, BondOrder::Single);
    mol.add_bond(18, 0, 17, BondOrder::Single);
    mol.add_bond(19, 1, 18, BondOrder::Single);
    mol.add_bond(20, 5, 19, BondOrder::Single);
    mol.add_bond(21, 5, 20, BondOrder::Single);
    mol.add_bond(22, 7, 21, BondOrder::Single);
    mol.add_bond(23, 8, 22, BondOrder::Single);
    mol.add_bond(24, 11, 23, BondOrder::Single);
    mol.add_bond(25, 12, 24, BondOrder::Single);
    mol.add_bond(26, 13, 25, BondOrder::Single);
    mol.add_bond(27, 14, 26, BondOrder::Single);

    mol
}

pub fn build_tryptophan_zwitterion_kekule() -> TestMolecule {
    let mut mol = TestMolecule::new();

    mol.add_atom(0, Element::N, 1);
    mol.add_atom(1, Element::C, 0);
    mol.add_atom(2, Element::C, 0);
    mol.add_atom(3, Element::O, -1);
    mol.add_atom(4, Element::O, 0);
    mol.add_atom(5, Element::C, 0);
    mol.add_atom(6, Element::C, 0);
    mol.add_atom(7, Element::C, 0);
    mol.add_atom(8, Element::N, 0);
    mol.add_atom(9, Element::C, 0);
    mol.add_atom(10, Element::C, 0);
    mol.add_atom(11, Element::C, 0);
    mol.add_atom(12, Element::C, 0);
    mol.add_atom(13, Element::C, 0);
    mol.add_atom(14, Element::C, 0);
    mol.add_atom(15, Element::H, 0);
    mol.add_atom(16, Element::H, 0);
    mol.add_atom(17, Element::H, 0);
    mol.add_atom(18, Element::H, 0);
    mol.add_atom(19, Element::H, 0);
    mol.add_atom(20, Element::H, 0);
    mol.add_atom(21, Element::H, 0);
    mol.add_atom(22, Element::H, 0);
    mol.add_atom(23, Element::H, 0);
    mol.add_atom(24, Element::H, 0);
    mol.add_atom(25, Element::H, 0);
    mol.add_atom(26, Element::H, 0);

    mol.add_bond(0, 0, 1, BondOrder::Single);
    mol.add_bond(1, 1, 2, BondOrder::Single);
    mol.add_bond(2, 2, 3, BondOrder::Single);
    mol.add_bond(3, 2, 4, BondOrder::Double);
    mol.add_bond(4, 1, 5, BondOrder::Single);
    mol.add_bond(5, 5, 6, BondOrder::Single);
    mol.add_bond(6, 6, 7, BondOrder::Double);
    mol.add_bond(7, 7, 8, BondOrder::Single);
    mol.add_bond(8, 8, 9, BondOrder::Single);
    mol.add_bond(9, 9, 10, BondOrder::Double);
    mol.add_bond(10, 10, 6, BondOrder::Single);
    mol.add_bond(11, 10, 11, BondOrder::Single);
    mol.add_bond(12, 11, 12, BondOrder::Double);
    mol.add_bond(13, 12, 13, BondOrder::Single);
    mol.add_bond(14, 13, 14, BondOrder::Double);
    mol.add_bond(15, 14, 9, BondOrder::Single);
    mol.add_bond(16, 0, 15, BondOrder::Single);
    mol.add_bond(17, 0, 16, BondOrder::Single);
    mol.add_bond(18, 0, 17, BondOrder::Single);
    mol.add_bond(19, 1, 18, BondOrder::Single);
    mol.add_bond(20, 5, 19, BondOrder::Single);
    mol.add_bond(21, 5, 20, BondOrder::Single);
    mol.add_bond(22, 7, 21, BondOrder::Single);
    mol.add_bond(23, 8, 22, BondOrder::Single);
    mol.add_bond(24, 11, 23, BondOrder::Single);
    mol.add_bond(25, 12, 24, BondOrder::Single);
    mol.add_bond(26, 13, 25, BondOrder::Single);
    mol.add_bond(27, 14, 26, BondOrder::Single);

    mol
}

pub fn build_uracil_aromatic() -> TestMolecule {
    let mut mol = TestMolecule::new();

    mol.add_atom(0, Element::N, 0);
    mol.add_atom(1, Element::C, 0);
    mol.add_atom(2, Element::N, 0);
    mol.add_atom(3, Element::C, 0);
    mol.add_atom(4, Element::C, 0);
    mol.add_atom(5, Element::C, 0);
    mol.add_atom(6, Element::O, 0);
    mol.add_atom(7, Element::O, 0);
    mol.add_atom(8, Element::H, 0);
    mol.add_atom(9, Element::H, 0);
    mol.add_atom(10, Element::H, 0);
    mol.add_atom(11, Element::H, 0);

    mol.add_bond(0, 0, 1, BondOrder::Aromatic);
    mol.add_bond(1, 1, 2, BondOrder::Aromatic);
    mol.add_bond(2, 2, 3, BondOrder::Aromatic);
    mol.add_bond(3, 3, 4, BondOrder::Aromatic);
    mol.add_bond(4, 4, 5, BondOrder::Aromatic);
    mol.add_bond(5, 5, 0, BondOrder::Aromatic);
    mol.add_bond(6, 1, 6, BondOrder::Double);
    mol.add_bond(7, 3, 7, BondOrder::Double);
    mol.add_bond(8, 0, 8, BondOrder::Single);
    mol.add_bond(9, 2, 9, BondOrder::Single);
    mol.add_bond(10, 4, 10, BondOrder::Single);
    mol.add_bond(11, 5, 11, BondOrder::Single);

    mol
}

pub fn build_uracil_kekule() -> TestMolecule {
    let mut mol = TestMolecule::new();

    mol.add_atom(0, Element::N, 0);
    mol.add_atom(1, Element::C, 0);
    mol.add_atom(2, Element::N, 0);
    mol.add_atom(3, Element::C, 0);
    mol.add_atom(4, Element::C, 0);
    mol.add_atom(5, Element::C, 0);
    mol.add_atom(6, Element::O, 0);
    mol.add_atom(7, Element::O, 0);
    mol.add_atom(8, Element::H, 0);
    mol.add_atom(9, Element::H, 0);
    mol.add_atom(10, Element::H, 0);
    mol.add_atom(11, Element::H, 0);

    mol.add_bond(0, 0, 1, BondOrder::Single);
    mol.add_bond(1, 1, 2, BondOrder::Double);
    mol.add_bond(2, 2, 3, BondOrder::Single);
    mol.add_bond(3, 3, 4, BondOrder::Double);
    mol.add_bond(4, 4, 5, BondOrder::Single);
    mol.add_bond(5, 5, 0, BondOrder::Double);
    mol.add_bond(6, 1, 6, BondOrder::Double);
    mol.add_bond(7, 3, 7, BondOrder::Double);
    mol.add_bond(8, 0, 8, BondOrder::Single);
    mol.add_bond(9, 2, 9, BondOrder::Single);
    mol.add_bond(10, 4, 10, BondOrder::Single);
    mol.add_bond(11, 5, 11, BondOrder::Single);

    mol
}

pub fn build_thymine_aromatic() -> TestMolecule {
    let mut mol = TestMolecule::new();

    mol.add_atom(0, Element::N, 0);
    mol.add_atom(1, Element::C, 0);
    mol.add_atom(2, Element::N, 0);
    mol.add_atom(3, Element::C, 0);
    mol.add_atom(4, Element::C, 0);
    mol.add_atom(5, Element::C, 0);
    mol.add_atom(6, Element::O, 0);
    mol.add_atom(7, Element::O, 0);
    mol.add_atom(8, Element::C, 0);
    mol.add_atom(9, Element::H, 0);
    mol.add_atom(10, Element::H, 0);
    mol.add_atom(11, Element::H, 0);
    mol.add_atom(12, Element::H, 0);
    mol.add_atom(13, Element::H, 0);
    mol.add_atom(14, Element::H, 0);

    mol.add_bond(0, 0, 1, BondOrder::Aromatic);
    mol.add_bond(1, 1, 2, BondOrder::Aromatic);
    mol.add_bond(2, 2, 3, BondOrder::Aromatic);
    mol.add_bond(3, 3, 4, BondOrder::Aromatic);
    mol.add_bond(4, 4, 5, BondOrder::Aromatic);
    mol.add_bond(5, 5, 0, BondOrder::Aromatic);
    mol.add_bond(6, 1, 6, BondOrder::Double);
    mol.add_bond(7, 3, 7, BondOrder::Double);
    mol.add_bond(8, 4, 8, BondOrder::Single);
    mol.add_bond(9, 0, 9, BondOrder::Single);
    mol.add_bond(10, 2, 10, BondOrder::Single);
    mol.add_bond(11, 5, 11, BondOrder::Single);
    mol.add_bond(12, 8, 12, BondOrder::Single);
    mol.add_bond(13, 8, 13, BondOrder::Single);
    mol.add_bond(14, 8, 14, BondOrder::Single);

    mol
}

pub fn build_thymine_kekule() -> TestMolecule {
    let mut mol = TestMolecule::new();

    mol.add_atom(0, Element::N, 0);
    mol.add_atom(1, Element::C, 0);
    mol.add_atom(2, Element::N, 0);
    mol.add_atom(3, Element::C, 0);
    mol.add_atom(4, Element::C, 0);
    mol.add_atom(5, Element::C, 0);
    mol.add_atom(6, Element::O, 0);
    mol.add_atom(7, Element::O, 0);
    mol.add_atom(8, Element::C, 0);
    mol.add_atom(9, Element::H, 0);
    mol.add_atom(10, Element::H, 0);
    mol.add_atom(11, Element::H, 0);
    mol.add_atom(12, Element::H, 0);
    mol.add_atom(13, Element::H, 0);
    mol.add_atom(14, Element::H, 0);

    mol.add_bond(0, 0, 1, BondOrder::Single);
    mol.add_bond(1, 1, 2, BondOrder::Single);
    mol.add_bond(2, 2, 3, BondOrder::Single);
    mol.add_bond(3, 3, 4, BondOrder::Single);
    mol.add_bond(4, 4, 5, BondOrder::Double);
    mol.add_bond(5, 5, 0, BondOrder::Single);
    mol.add_bond(6, 1, 6, BondOrder::Double);
    mol.add_bond(7, 3, 7, BondOrder::Double);
    mol.add_bond(8, 4, 8, BondOrder::Single);
    mol.add_bond(9, 0, 9, BondOrder::Single);
    mol.add_bond(10, 2, 10, BondOrder::Single);
    mol.add_bond(11, 5, 11, BondOrder::Single);
    mol.add_bond(12, 8, 12, BondOrder::Single);
    mol.add_bond(13, 8, 13, BondOrder::Single);
    mol.add_bond(14, 8, 14, BondOrder::Single);

    mol
}

pub fn build_cytosine_aromatic() -> TestMolecule {
    let mut mol = TestMolecule::new();

    mol.add_atom(0, Element::N, 0);
    mol.add_atom(1, Element::C, 0);
    mol.add_atom(2, Element::N, 0);
    mol.add_atom(3, Element::C, 0);
    mol.add_atom(4, Element::C, 0);
    mol.add_atom(5, Element::C, 0);
    mol.add_atom(6, Element::O, 0);
    mol.add_atom(7, Element::N, 0);
    mol.add_atom(8, Element::H, 0);
    mol.add_atom(9, Element::H, 0);
    mol.add_atom(10, Element::H, 0);
    mol.add_atom(11, Element::H, 0);
    mol.add_atom(12, Element::H, 0);

    mol.add_bond(0, 0, 1, BondOrder::Aromatic);
    mol.add_bond(1, 1, 2, BondOrder::Aromatic);
    mol.add_bond(2, 2, 3, BondOrder::Aromatic);
    mol.add_bond(3, 3, 4, BondOrder::Aromatic);
    mol.add_bond(4, 4, 5, BondOrder::Aromatic);
    mol.add_bond(5, 5, 0, BondOrder::Aromatic);
    mol.add_bond(6, 1, 6, BondOrder::Double);
    mol.add_bond(7, 3, 7, BondOrder::Single);
    mol.add_bond(8, 0, 8, BondOrder::Single);
    mol.add_bond(9, 7, 9, BondOrder::Single);
    mol.add_bond(10, 7, 10, BondOrder::Single);
    mol.add_bond(11, 4, 11, BondOrder::Single);
    mol.add_bond(12, 5, 12, BondOrder::Single);

    mol
}

pub fn build_cytosine_kekule() -> TestMolecule {
    let mut mol = TestMolecule::new();

    mol.add_atom(0, Element::N, 0);
    mol.add_atom(1, Element::C, 0);
    mol.add_atom(2, Element::N, 0);
    mol.add_atom(3, Element::C, 0);
    mol.add_atom(4, Element::C, 0);
    mol.add_atom(5, Element::C, 0);
    mol.add_atom(6, Element::O, 0);
    mol.add_atom(7, Element::N, 0);
    mol.add_atom(8, Element::H, 0);
    mol.add_atom(9, Element::H, 0);
    mol.add_atom(10, Element::H, 0);
    mol.add_atom(11, Element::H, 0);
    mol.add_atom(12, Element::H, 0);

    mol.add_bond(0, 0, 1, BondOrder::Single);
    mol.add_bond(1, 1, 2, BondOrder::Single);
    mol.add_bond(2, 2, 3, BondOrder::Double);
    mol.add_bond(3, 3, 4, BondOrder::Single);
    mol.add_bond(4, 4, 5, BondOrder::Double);
    mol.add_bond(5, 5, 0, BondOrder::Single);
    mol.add_bond(6, 1, 6, BondOrder::Double);
    mol.add_bond(7, 3, 7, BondOrder::Single);
    mol.add_bond(8, 0, 8, BondOrder::Single);
    mol.add_bond(9, 7, 9, BondOrder::Single);
    mol.add_bond(10, 7, 10, BondOrder::Single);
    mol.add_bond(11, 4, 11, BondOrder::Single);
    mol.add_bond(12, 5, 12, BondOrder::Single);

    mol
}

pub fn build_adenine_aromatic() -> TestMolecule {
    let mut mol = TestMolecule::new();

    mol.add_atom(0, Element::N, 0);
    mol.add_atom(1, Element::C, 0);
    mol.add_atom(2, Element::N, 0);
    mol.add_atom(3, Element::C, 0);
    mol.add_atom(4, Element::C, 0);
    mol.add_atom(5, Element::C, 0);
    mol.add_atom(6, Element::N, 0);
    mol.add_atom(7, Element::C, 0);
    mol.add_atom(8, Element::N, 0);
    mol.add_atom(9, Element::N, 0);
    mol.add_atom(10, Element::H, 0);
    mol.add_atom(11, Element::H, 0);
    mol.add_atom(12, Element::H, 0);
    mol.add_atom(13, Element::H, 0);
    mol.add_atom(14, Element::H, 0);

    mol.add_bond(0, 8, 3, BondOrder::Aromatic);
    mol.add_bond(1, 3, 2, BondOrder::Aromatic);
    mol.add_bond(2, 2, 1, BondOrder::Aromatic);
    mol.add_bond(3, 1, 0, BondOrder::Aromatic);
    mol.add_bond(4, 0, 5, BondOrder::Aromatic);
    mol.add_bond(5, 5, 4, BondOrder::Aromatic);
    mol.add_bond(6, 4, 3, BondOrder::Aromatic);
    mol.add_bond(7, 4, 6, BondOrder::Aromatic);
    mol.add_bond(8, 6, 7, BondOrder::Aromatic);
    mol.add_bond(9, 7, 8, BondOrder::Aromatic);
    mol.add_bond(10, 5, 9, BondOrder::Single);
    mol.add_bond(11, 1, 10, BondOrder::Single);
    mol.add_bond(12, 7, 11, BondOrder::Single);
    mol.add_bond(13, 8, 12, BondOrder::Single);
    mol.add_bond(14, 9, 13, BondOrder::Single);
    mol.add_bond(15, 9, 14, BondOrder::Single);

    mol
}

pub fn build_adenine_kekule() -> TestMolecule {
    let mut mol = TestMolecule::new();

    mol.add_atom(0, Element::N, 0);
    mol.add_atom(1, Element::C, 0);
    mol.add_atom(2, Element::N, 0);
    mol.add_atom(3, Element::C, 0);
    mol.add_atom(4, Element::C, 0);
    mol.add_atom(5, Element::C, 0);
    mol.add_atom(6, Element::N, 0);
    mol.add_atom(7, Element::C, 0);
    mol.add_atom(8, Element::N, 0);
    mol.add_atom(9, Element::N, 0);
    mol.add_atom(10, Element::H, 0);
    mol.add_atom(11, Element::H, 0);
    mol.add_atom(12, Element::H, 0);
    mol.add_atom(13, Element::H, 0);
    mol.add_atom(14, Element::H, 0);

    mol.add_bond(0, 8, 3, BondOrder::Single);
    mol.add_bond(1, 3, 2, BondOrder::Double);
    mol.add_bond(2, 2, 1, BondOrder::Single);
    mol.add_bond(3, 1, 0, BondOrder::Double);
    mol.add_bond(4, 0, 5, BondOrder::Single);
    mol.add_bond(5, 5, 4, BondOrder::Double);
    mol.add_bond(6, 4, 3, BondOrder::Single);
    mol.add_bond(7, 4, 6, BondOrder::Single);
    mol.add_bond(8, 6, 7, BondOrder::Double);
    mol.add_bond(9, 7, 8, BondOrder::Single);
    mol.add_bond(10, 5, 9, BondOrder::Single);
    mol.add_bond(11, 1, 10, BondOrder::Single);
    mol.add_bond(12, 7, 11, BondOrder::Single);
    mol.add_bond(13, 8, 12, BondOrder::Single);
    mol.add_bond(14, 9, 13, BondOrder::Single);
    mol.add_bond(15, 9, 14, BondOrder::Single);

    mol
}

pub fn build_guanine_aromatic() -> TestMolecule {
    let mut mol = TestMolecule::new();

    mol.add_atom(0, Element::N, 0);
    mol.add_atom(1, Element::C, 0);
    mol.add_atom(2, Element::N, 0);
    mol.add_atom(3, Element::C, 0);
    mol.add_atom(4, Element::C, 0);
    mol.add_atom(5, Element::C, 0);
    mol.add_atom(6, Element::N, 0);
    mol.add_atom(7, Element::C, 0);
    mol.add_atom(8, Element::N, 0);
    mol.add_atom(9, Element::O, 0);
    mol.add_atom(10, Element::N, 0);
    mol.add_atom(11, Element::H, 0);
    mol.add_atom(12, Element::H, 0);
    mol.add_atom(13, Element::H, 0);
    mol.add_atom(14, Element::H, 0);
    mol.add_atom(15, Element::H, 0);

    mol.add_bond(0, 8, 3, BondOrder::Aromatic);
    mol.add_bond(1, 3, 2, BondOrder::Aromatic);
    mol.add_bond(2, 2, 1, BondOrder::Aromatic);
    mol.add_bond(3, 1, 0, BondOrder::Aromatic);
    mol.add_bond(4, 0, 5, BondOrder::Aromatic);
    mol.add_bond(5, 5, 4, BondOrder::Aromatic);
    mol.add_bond(6, 4, 3, BondOrder::Aromatic);
    mol.add_bond(7, 4, 6, BondOrder::Aromatic);
    mol.add_bond(8, 6, 7, BondOrder::Aromatic);
    mol.add_bond(9, 7, 8, BondOrder::Aromatic);
    mol.add_bond(10, 5, 9, BondOrder::Double);
    mol.add_bond(11, 1, 10, BondOrder::Single);
    mol.add_bond(12, 0, 11, BondOrder::Single);
    mol.add_bond(13, 7, 12, BondOrder::Single);
    mol.add_bond(14, 8, 13, BondOrder::Single);
    mol.add_bond(15, 10, 14, BondOrder::Single);
    mol.add_bond(16, 10, 15, BondOrder::Single);

    mol
}

pub fn build_guanine_kekule() -> TestMolecule {
    let mut mol = TestMolecule::new();

    mol.add_atom(0, Element::N, 0);
    mol.add_atom(1, Element::C, 0);
    mol.add_atom(2, Element::N, 0);
    mol.add_atom(3, Element::C, 0);
    mol.add_atom(4, Element::C, 0);
    mol.add_atom(5, Element::C, 0);
    mol.add_atom(6, Element::N, 0);
    mol.add_atom(7, Element::C, 0);
    mol.add_atom(8, Element::N, 0);
    mol.add_atom(9, Element::O, 0);
    mol.add_atom(10, Element::N, 0);
    mol.add_atom(11, Element::H, 0);
    mol.add_atom(12, Element::H, 0);
    mol.add_atom(13, Element::H, 0);
    mol.add_atom(14, Element::H, 0);
    mol.add_atom(15, Element::H, 0);

    mol.add_bond(0, 8, 3, BondOrder::Single);
    mol.add_bond(1, 3, 2, BondOrder::Single);
    mol.add_bond(2, 2, 1, BondOrder::Double);
    mol.add_bond(3, 1, 0, BondOrder::Single);
    mol.add_bond(4, 0, 5, BondOrder::Single);
    mol.add_bond(5, 5, 4, BondOrder::Double);
    mol.add_bond(6, 4, 3, BondOrder::Double);
    mol.add_bond(7, 4, 6, BondOrder::Single);
    mol.add_bond(8, 6, 7, BondOrder::Double);
    mol.add_bond(9, 7, 8, BondOrder::Single);
    mol.add_bond(10, 5, 9, BondOrder::Double);
    mol.add_bond(11, 1, 10, BondOrder::Single);
    mol.add_bond(12, 0, 11, BondOrder::Single);
    mol.add_bond(13, 7, 12, BondOrder::Single);
    mol.add_bond(14, 8, 13, BondOrder::Single);
    mol.add_bond(15, 10, 14, BondOrder::Single);
    mol.add_bond(16, 10, 15, BondOrder::Single);

    mol
}

pub fn build_deoxyadenosine_aromatic() -> TestMolecule {
    let mut mol = TestMolecule::new();

    mol.add_atom(0, Element::N, 0);
    mol.add_atom(1, Element::C, 0);
    mol.add_atom(2, Element::N, 0);
    mol.add_atom(3, Element::C, 0);
    mol.add_atom(4, Element::C, 0);
    mol.add_atom(5, Element::C, 0);
    mol.add_atom(6, Element::N, 0);
    mol.add_atom(7, Element::C, 0);
    mol.add_atom(8, Element::N, 0);
    mol.add_atom(9, Element::N, 0);
    mol.add_atom(10, Element::H, 0);
    mol.add_atom(11, Element::H, 0);
    mol.add_atom(12, Element::H, 0);
    mol.add_atom(13, Element::H, 0);
    mol.add_atom(14, Element::C, 0);
    mol.add_atom(15, Element::C, 0);
    mol.add_atom(16, Element::C, 0);
    mol.add_atom(17, Element::C, 0);
    mol.add_atom(18, Element::C, 0);
    mol.add_atom(19, Element::O, 0);
    mol.add_atom(20, Element::O, 0);
    mol.add_atom(21, Element::O, 0);
    mol.add_atom(22, Element::H, 0);
    mol.add_atom(23, Element::H, 0);
    mol.add_atom(24, Element::H, 0);
    mol.add_atom(25, Element::H, 0);
    mol.add_atom(26, Element::H, 0);
    mol.add_atom(27, Element::H, 0);
    mol.add_atom(28, Element::H, 0);
    mol.add_atom(29, Element::H, 0);
    mol.add_atom(30, Element::H, 0);

    mol.add_bond(0, 8, 3, BondOrder::Aromatic);
    mol.add_bond(1, 3, 2, BondOrder::Aromatic);
    mol.add_bond(2, 2, 1, BondOrder::Aromatic);
    mol.add_bond(3, 1, 0, BondOrder::Aromatic);
    mol.add_bond(4, 0, 5, BondOrder::Aromatic);
    mol.add_bond(5, 5, 4, BondOrder::Aromatic);
    mol.add_bond(6, 4, 3, BondOrder::Aromatic);
    mol.add_bond(7, 4, 6, BondOrder::Aromatic);
    mol.add_bond(8, 6, 7, BondOrder::Aromatic);
    mol.add_bond(9, 7, 8, BondOrder::Aromatic);
    mol.add_bond(10, 5, 9, BondOrder::Single);
    mol.add_bond(11, 9, 12, BondOrder::Single);
    mol.add_bond(12, 9, 13, BondOrder::Single);
    mol.add_bond(13, 1, 10, BondOrder::Single);
    mol.add_bond(14, 7, 11, BondOrder::Single);
    mol.add_bond(15, 8, 14, BondOrder::Single);
    mol.add_bond(16, 14, 15, BondOrder::Single);
    mol.add_bond(17, 15, 16, BondOrder::Single);
    mol.add_bond(18, 16, 17, BondOrder::Single);
    mol.add_bond(19, 17, 19, BondOrder::Single);
    mol.add_bond(20, 19, 14, BondOrder::Single);
    mol.add_bond(21, 16, 20, BondOrder::Single);
    mol.add_bond(22, 17, 18, BondOrder::Single);
    mol.add_bond(23, 18, 21, BondOrder::Single);
    mol.add_bond(24, 14, 22, BondOrder::Single);
    mol.add_bond(25, 15, 23, BondOrder::Single);
    mol.add_bond(26, 15, 24, BondOrder::Single);
    mol.add_bond(27, 16, 25, BondOrder::Single);
    mol.add_bond(28, 17, 26, BondOrder::Single);
    mol.add_bond(29, 18, 27, BondOrder::Single);
    mol.add_bond(30, 18, 28, BondOrder::Single);
    mol.add_bond(31, 20, 29, BondOrder::Single);
    mol.add_bond(32, 21, 30, BondOrder::Single);

    mol
}

pub fn build_deoxyadenosine_kekule() -> TestMolecule {
    let mut mol = TestMolecule::new();

    mol.add_atom(0, Element::N, 0);
    mol.add_atom(1, Element::C, 0);
    mol.add_atom(2, Element::N, 0);
    mol.add_atom(3, Element::C, 0);
    mol.add_atom(4, Element::C, 0);
    mol.add_atom(5, Element::C, 0);
    mol.add_atom(6, Element::N, 0);
    mol.add_atom(7, Element::C, 0);
    mol.add_atom(8, Element::N, 0);
    mol.add_atom(9, Element::N, 0);
    mol.add_atom(10, Element::H, 0);
    mol.add_atom(11, Element::H, 0);
    mol.add_atom(12, Element::H, 0);
    mol.add_atom(13, Element::H, 0);
    mol.add_atom(14, Element::C, 0);
    mol.add_atom(15, Element::C, 0);
    mol.add_atom(16, Element::C, 0);
    mol.add_atom(17, Element::C, 0);
    mol.add_atom(18, Element::C, 0);
    mol.add_atom(19, Element::O, 0);
    mol.add_atom(20, Element::O, 0);
    mol.add_atom(21, Element::O, 0);
    mol.add_atom(22, Element::H, 0);
    mol.add_atom(23, Element::H, 0);
    mol.add_atom(24, Element::H, 0);
    mol.add_atom(25, Element::H, 0);
    mol.add_atom(26, Element::H, 0);
    mol.add_atom(27, Element::H, 0);
    mol.add_atom(28, Element::H, 0);
    mol.add_atom(29, Element::H, 0);
    mol.add_atom(30, Element::H, 0);

    mol.add_bond(0, 8, 3, BondOrder::Single);
    mol.add_bond(1, 3, 2, BondOrder::Double);
    mol.add_bond(2, 2, 1, BondOrder::Single);
    mol.add_bond(3, 1, 0, BondOrder::Double);
    mol.add_bond(4, 0, 5, BondOrder::Single);
    mol.add_bond(5, 5, 4, BondOrder::Double);
    mol.add_bond(6, 4, 3, BondOrder::Single);
    mol.add_bond(7, 4, 6, BondOrder::Single);
    mol.add_bond(8, 6, 7, BondOrder::Double);
    mol.add_bond(9, 7, 8, BondOrder::Single);
    mol.add_bond(10, 5, 9, BondOrder::Single);
    mol.add_bond(11, 9, 12, BondOrder::Single);
    mol.add_bond(12, 9, 13, BondOrder::Single);
    mol.add_bond(13, 1, 10, BondOrder::Single);
    mol.add_bond(14, 7, 11, BondOrder::Single);
    mol.add_bond(15, 8, 14, BondOrder::Single);
    mol.add_bond(16, 14, 15, BondOrder::Single);
    mol.add_bond(17, 15, 16, BondOrder::Single);
    mol.add_bond(18, 16, 17, BondOrder::Single);
    mol.add_bond(19, 17, 19, BondOrder::Single);
    mol.add_bond(20, 19, 14, BondOrder::Single);
    mol.add_bond(21, 16, 20, BondOrder::Single);
    mol.add_bond(22, 17, 18, BondOrder::Single);
    mol.add_bond(23, 18, 21, BondOrder::Single);
    mol.add_bond(24, 14, 22, BondOrder::Single);
    mol.add_bond(25, 15, 23, BondOrder::Single);
    mol.add_bond(26, 15, 24, BondOrder::Single);
    mol.add_bond(27, 16, 25, BondOrder::Single);
    mol.add_bond(28, 17, 26, BondOrder::Single);
    mol.add_bond(29, 18, 27, BondOrder::Single);
    mol.add_bond(30, 18, 28, BondOrder::Single);
    mol.add_bond(31, 20, 29, BondOrder::Single);
    mol.add_bond(32, 21, 30, BondOrder::Single);

    mol
}

pub fn build_dinucleotide_backbone() -> TestMolecule {
    let mut mol = TestMolecule::new();

    mol.add_atom(0, Element::C, 0);
    mol.add_atom(1, Element::O, 0);
    mol.add_atom(2, Element::P, 0);
    mol.add_atom(3, Element::O, 0);
    mol.add_atom(4, Element::O, -1);
    mol.add_atom(5, Element::O, 0);
    mol.add_atom(6, Element::C, 0);
    mol.add_atom(7, Element::H, 0);
    mol.add_atom(8, Element::H, 0);
    mol.add_atom(9, Element::H, 0);
    mol.add_atom(10, Element::H, 0);
    mol.add_atom(11, Element::H, 0);
    mol.add_atom(12, Element::H, 0);

    mol.add_bond(0, 0, 1, BondOrder::Single);
    mol.add_bond(1, 1, 2, BondOrder::Single);
    mol.add_bond(2, 2, 5, BondOrder::Single);
    mol.add_bond(3, 5, 6, BondOrder::Single);
    mol.add_bond(4, 2, 3, BondOrder::Double);
    mol.add_bond(5, 2, 4, BondOrder::Single);
    mol.add_bond(6, 0, 7, BondOrder::Single);
    mol.add_bond(7, 0, 8, BondOrder::Single);
    mol.add_bond(8, 0, 9, BondOrder::Single);
    mol.add_bond(9, 6, 10, BondOrder::Single);
    mol.add_bond(10, 6, 11, BondOrder::Single);
    mol.add_bond(11, 6, 12, BondOrder::Single);

    mol
}

pub fn build_acbuol() -> TestMolecule {
    let mut mol = TestMolecule::new();

    mol.add_atom(0, Element::C, 0);
    mol.add_atom(1, Element::C, 0);
    mol.add_atom(2, Element::C, 0);
    mol.add_atom(3, Element::C, 0);
    mol.add_atom(4, Element::C, 0);
    mol.add_atom(5, Element::O, 0);
    mol.add_atom(6, Element::O, 0);
    mol.add_atom(7, Element::C, 0);
    mol.add_atom(8, Element::O, 0);
    mol.add_atom(9, Element::C, 0);
    mol.add_atom(10, Element::O, 0);
    mol.add_atom(11, Element::H, 0);
    mol.add_atom(12, Element::H, 0);
    mol.add_atom(13, Element::H, 0);
    mol.add_atom(14, Element::H, 0);
    mol.add_atom(15, Element::H, 0);
    mol.add_atom(16, Element::H, 0);
    mol.add_atom(17, Element::H, 0);
    mol.add_atom(18, Element::H, 0);
    mol.add_atom(19, Element::H, 0);
    mol.add_atom(20, Element::H, 0);

    mol.add_bond(0, 0, 1, BondOrder::Single);
    mol.add_bond(1, 1, 2, BondOrder::Single);
    mol.add_bond(2, 2, 3, BondOrder::Single);
    mol.add_bond(3, 3, 0, BondOrder::Single);
    mol.add_bond(4, 0, 4, BondOrder::Single);
    mol.add_bond(5, 1, 8, BondOrder::Single);
    mol.add_bond(6, 2, 10, BondOrder::Single);
    mol.add_bond(7, 3, 7, BondOrder::Single);
    mol.add_bond(8, 4, 5, BondOrder::Double);
    mol.add_bond(9, 4, 6, BondOrder::Single);
    mol.add_bond(10, 6, 7, BondOrder::Single);
    mol.add_bond(11, 8, 9, BondOrder::Single);
    mol.add_bond(12, 0, 11, BondOrder::Single);
    mol.add_bond(13, 1, 12, BondOrder::Single);
    mol.add_bond(14, 2, 13, BondOrder::Single);
    mol.add_bond(15, 3, 14, BondOrder::Single);
    mol.add_bond(16, 7, 15, BondOrder::Single);
    mol.add_bond(17, 7, 16, BondOrder::Single);
    mol.add_bond(18, 9, 17, BondOrder::Single);
    mol.add_bond(19, 9, 18, BondOrder::Single);
    mol.add_bond(20, 9, 19, BondOrder::Single);
    mol.add_bond(21, 10, 20, BondOrder::Single);

    mol
}

pub fn build_afurpo10() -> TestMolecule {
    let mut mol = TestMolecule::new();

    mol.add_atom(0, Element::C, 0);
    mol.add_atom(1, Element::C, 0);
    mol.add_atom(2, Element::C, 0);
    mol.add_atom(3, Element::C, 0);
    mol.add_atom(4, Element::C, 0);
    mol.add_atom(5, Element::C, 0);
    mol.add_atom(6, Element::O, 0);
    mol.add_atom(7, Element::O, 0);
    mol.add_atom(8, Element::H, 0);
    mol.add_atom(9, Element::H, 0);
    mol.add_atom(10, Element::H, 0);
    mol.add_atom(11, Element::H, 0);
    mol.add_atom(12, Element::H, 0);
    mol.add_atom(13, Element::H, 0);
    mol.add_atom(14, Element::H, 0);
    mol.add_atom(15, Element::H, 0);

    mol.add_bond(0, 0, 1, BondOrder::Single);
    mol.add_bond(1, 0, 5, BondOrder::Single);
    mol.add_bond(2, 0, 6, BondOrder::Single);
    mol.add_bond(3, 1, 2, BondOrder::Single);
    mol.add_bond(4, 2, 3, BondOrder::Single);
    mol.add_bond(5, 2, 7, BondOrder::Single);
    mol.add_bond(6, 3, 4, BondOrder::Single);
    mol.add_bond(7, 3, 6, BondOrder::Single);
    mol.add_bond(8, 4, 5, BondOrder::Single);
    mol.add_bond(9, 4, 7, BondOrder::Single);
    mol.add_bond(10, 0, 8, BondOrder::Single);
    mol.add_bond(11, 1, 9, BondOrder::Single);
    mol.add_bond(12, 1, 10, BondOrder::Single);
    mol.add_bond(13, 2, 11, BondOrder::Single);
    mol.add_bond(14, 3, 12, BondOrder::Single);
    mol.add_bond(15, 4, 13, BondOrder::Single);
    mol.add_bond(16, 5, 14, BondOrder::Single);
    mol.add_bond(17, 5, 15, BondOrder::Single);

    mol
}

pub fn build_adamantane() -> TestMolecule {
    let mut mol = TestMolecule::new();

    mol.add_atom(0, Element::C, 0);
    mol.add_atom(1, Element::C, 0);
    mol.add_atom(2, Element::C, 0);
    mol.add_atom(3, Element::C, 0);
    mol.add_atom(4, Element::C, 0);
    mol.add_atom(5, Element::C, 0);
    mol.add_atom(6, Element::C, 0);
    mol.add_atom(7, Element::C, 0);
    mol.add_atom(8, Element::C, 0);
    mol.add_atom(9, Element::C, 0);
    mol.add_atom(10, Element::H, 0);
    mol.add_atom(11, Element::H, 0);
    mol.add_atom(12, Element::H, 0);
    mol.add_atom(13, Element::H, 0);
    mol.add_atom(14, Element::H, 0);
    mol.add_atom(15, Element::H, 0);
    mol.add_atom(16, Element::H, 0);
    mol.add_atom(17, Element::H, 0);
    mol.add_atom(18, Element::H, 0);
    mol.add_atom(19, Element::H, 0);
    mol.add_atom(20, Element::H, 0);
    mol.add_atom(21, Element::H, 0);
    mol.add_atom(22, Element::H, 0);
    mol.add_atom(23, Element::H, 0);
    mol.add_atom(24, Element::H, 0);
    mol.add_atom(25, Element::H, 0);

    mol.add_bond(0, 0, 1, BondOrder::Single);
    mol.add_bond(1, 0, 3, BondOrder::Single);
    mol.add_bond(2, 0, 5, BondOrder::Single);
    mol.add_bond(3, 2, 1, BondOrder::Single);
    mol.add_bond(4, 2, 7, BondOrder::Single);
    mol.add_bond(5, 2, 9, BondOrder::Single);
    mol.add_bond(6, 4, 3, BondOrder::Single);
    mol.add_bond(7, 4, 7, BondOrder::Single);
    mol.add_bond(8, 4, 8, BondOrder::Single);
    mol.add_bond(9, 6, 5, BondOrder::Single);
    mol.add_bond(10, 6, 8, BondOrder::Single);
    mol.add_bond(11, 6, 9, BondOrder::Single);
    mol.add_bond(12, 0, 10, BondOrder::Single);
    mol.add_bond(13, 2, 11, BondOrder::Single);
    mol.add_bond(14, 4, 12, BondOrder::Single);
    mol.add_bond(15, 6, 13, BondOrder::Single);
    mol.add_bond(16, 1, 14, BondOrder::Single);
    mol.add_bond(17, 1, 15, BondOrder::Single);
    mol.add_bond(18, 3, 16, BondOrder::Single);
    mol.add_bond(19, 3, 17, BondOrder::Single);
    mol.add_bond(20, 5, 18, BondOrder::Single);
    mol.add_bond(21, 5, 19, BondOrder::Single);
    mol.add_bond(22, 7, 20, BondOrder::Single);
    mol.add_bond(23, 7, 21, BondOrder::Single);
    mol.add_bond(24, 8, 22, BondOrder::Single);
    mol.add_bond(25, 8, 23, BondOrder::Single);
    mol.add_bond(26, 9, 24, BondOrder::Single);
    mol.add_bond(27, 9, 25, BondOrder::Single);

    mol
}

pub fn build_decalin() -> TestMolecule {
    let mut mol = TestMolecule::new();

    mol.add_atom(0, Element::C, 0);
    mol.add_atom(1, Element::C, 0);
    mol.add_atom(2, Element::C, 0);
    mol.add_atom(3, Element::C, 0);
    mol.add_atom(4, Element::C, 0);
    mol.add_atom(5, Element::C, 0);
    mol.add_atom(6, Element::C, 0);
    mol.add_atom(7, Element::C, 0);
    mol.add_atom(8, Element::C, 0);
    mol.add_atom(9, Element::C, 0);
    mol.add_atom(10, Element::H, 0);
    mol.add_atom(11, Element::H, 0);
    mol.add_atom(12, Element::H, 0);
    mol.add_atom(13, Element::H, 0);
    mol.add_atom(14, Element::H, 0);
    mol.add_atom(15, Element::H, 0);
    mol.add_atom(16, Element::H, 0);
    mol.add_atom(17, Element::H, 0);
    mol.add_atom(18, Element::H, 0);
    mol.add_atom(19, Element::H, 0);
    mol.add_atom(20, Element::H, 0);
    mol.add_atom(21, Element::H, 0);
    mol.add_atom(22, Element::H, 0);
    mol.add_atom(23, Element::H, 0);
    mol.add_atom(24, Element::H, 0);
    mol.add_atom(25, Element::H, 0);
    mol.add_atom(26, Element::H, 0);
    mol.add_atom(27, Element::H, 0);

    mol.add_bond(0, 0, 1, BondOrder::Single);
    mol.add_bond(1, 1, 2, BondOrder::Single);
    mol.add_bond(2, 2, 3, BondOrder::Single);
    mol.add_bond(3, 3, 4, BondOrder::Single);
    mol.add_bond(4, 4, 5, BondOrder::Single);
    mol.add_bond(5, 5, 0, BondOrder::Single);
    mol.add_bond(6, 5, 6, BondOrder::Single);
    mol.add_bond(7, 6, 7, BondOrder::Single);
    mol.add_bond(8, 7, 8, BondOrder::Single);
    mol.add_bond(9, 8, 9, BondOrder::Single);
    mol.add_bond(10, 9, 0, BondOrder::Single);
    mol.add_bond(11, 0, 10, BondOrder::Single);
    mol.add_bond(12, 1, 11, BondOrder::Single);
    mol.add_bond(13, 1, 12, BondOrder::Single);
    mol.add_bond(14, 2, 13, BondOrder::Single);
    mol.add_bond(15, 2, 14, BondOrder::Single);
    mol.add_bond(16, 3, 15, BondOrder::Single);
    mol.add_bond(17, 3, 16, BondOrder::Single);
    mol.add_bond(18, 4, 17, BondOrder::Single);
    mol.add_bond(19, 4, 18, BondOrder::Single);
    mol.add_bond(20, 5, 19, BondOrder::Single);
    mol.add_bond(21, 6, 20, BondOrder::Single);
    mol.add_bond(22, 6, 21, BondOrder::Single);
    mol.add_bond(23, 7, 22, BondOrder::Single);
    mol.add_bond(24, 7, 23, BondOrder::Single);
    mol.add_bond(25, 8, 24, BondOrder::Single);
    mol.add_bond(26, 8, 25, BondOrder::Single);
    mol.add_bond(27, 9, 26, BondOrder::Single);
    mol.add_bond(28, 9, 27, BondOrder::Single);

    mol
}

pub fn build_abaxes() -> TestMolecule {
    let mut mol = TestMolecule::new();

    mol.add_atom(0, Element::S, 0);
    mol.add_atom(1, Element::C, 0);
    mol.add_atom(2, Element::C, 0);
    mol.add_atom(3, Element::S, 0);
    mol.add_atom(4, Element::C, 0);
    mol.add_atom(5, Element::C, 0);
    mol.add_atom(6, Element::H, 0);
    mol.add_atom(7, Element::H, 0);
    mol.add_atom(8, Element::H, 0);
    mol.add_atom(9, Element::H, 0);
    mol.add_atom(10, Element::H, 0);
    mol.add_atom(11, Element::H, 0);

    mol.add_bond(0, 0, 1, BondOrder::Single);
    mol.add_bond(1, 1, 2, BondOrder::Double);
    mol.add_bond(2, 2, 3, BondOrder::Single);
    mol.add_bond(3, 3, 4, BondOrder::Single);
    mol.add_bond(4, 4, 5, BondOrder::Single);
    mol.add_bond(5, 5, 0, BondOrder::Single);
    mol.add_bond(6, 1, 6, BondOrder::Single);
    mol.add_bond(7, 2, 7, BondOrder::Single);
    mol.add_bond(8, 4, 8, BondOrder::Single);
    mol.add_bond(9, 4, 9, BondOrder::Single);
    mol.add_bond(10, 5, 10, BondOrder::Single);
    mol.add_bond(11, 5, 11, BondOrder::Single);

    mol
}

pub fn build_4_methylpyridine_aromatic() -> TestMolecule {
    let mut mol = TestMolecule::new();

    mol.add_atom(0, Element::N, 0);
    mol.add_atom(1, Element::C, 0);
    mol.add_atom(2, Element::C, 0);
    mol.add_atom(3, Element::C, 0);
    mol.add_atom(4, Element::C, 0);
    mol.add_atom(5, Element::C, 0);
    mol.add_atom(6, Element::C, 0);
    mol.add_atom(7, Element::H, 0);
    mol.add_atom(8, Element::H, 0);
    mol.add_atom(9, Element::H, 0);
    mol.add_atom(10, Element::H, 0);
    mol.add_atom(11, Element::H, 0);
    mol.add_atom(12, Element::H, 0);
    mol.add_atom(13, Element::H, 0);

    mol.add_bond(0, 0, 1, BondOrder::Aromatic);
    mol.add_bond(1, 1, 2, BondOrder::Aromatic);
    mol.add_bond(2, 2, 3, BondOrder::Aromatic);
    mol.add_bond(3, 3, 4, BondOrder::Aromatic);
    mol.add_bond(4, 4, 5, BondOrder::Aromatic);
    mol.add_bond(5, 5, 0, BondOrder::Aromatic);
    mol.add_bond(6, 3, 6, BondOrder::Single);
    mol.add_bond(7, 1, 7, BondOrder::Single);
    mol.add_bond(8, 2, 8, BondOrder::Single);
    mol.add_bond(9, 4, 9, BondOrder::Single);
    mol.add_bond(10, 5, 10, BondOrder::Single);
    mol.add_bond(11, 6, 11, BondOrder::Single);
    mol.add_bond(12, 6, 12, BondOrder::Single);
    mol.add_bond(13, 6, 13, BondOrder::Single);

    mol
}

pub fn build_4_methylpyridine_kekule() -> TestMolecule {
    let mut mol = TestMolecule::new();

    mol.add_atom(0, Element::N, 0);
    mol.add_atom(1, Element::C, 0);
    mol.add_atom(2, Element::C, 0);
    mol.add_atom(3, Element::C, 0);
    mol.add_atom(4, Element::C, 0);
    mol.add_atom(5, Element::C, 0);
    mol.add_atom(6, Element::C, 0);
    mol.add_atom(7, Element::H, 0);
    mol.add_atom(8, Element::H, 0);
    mol.add_atom(9, Element::H, 0);
    mol.add_atom(10, Element::H, 0);
    mol.add_atom(11, Element::H, 0);
    mol.add_atom(12, Element::H, 0);
    mol.add_atom(13, Element::H, 0);

    mol.add_bond(0, 0, 1, BondOrder::Double);
    mol.add_bond(1, 1, 2, BondOrder::Single);
    mol.add_bond(2, 2, 3, BondOrder::Double);
    mol.add_bond(3, 3, 4, BondOrder::Single);
    mol.add_bond(4, 4, 5, BondOrder::Double);
    mol.add_bond(5, 5, 0, BondOrder::Single);
    mol.add_bond(6, 3, 6, BondOrder::Single);
    mol.add_bond(7, 1, 7, BondOrder::Single);
    mol.add_bond(8, 2, 8, BondOrder::Single);
    mol.add_bond(9, 4, 9, BondOrder::Single);
    mol.add_bond(10, 5, 10, BondOrder::Single);
    mol.add_bond(11, 6, 11, BondOrder::Single);
    mol.add_bond(12, 6, 12, BondOrder::Single);
    mol.add_bond(13, 6, 13, BondOrder::Single);

    mol
}

pub fn build_acridine_aromatic() -> TestMolecule {
    let mut mol = TestMolecule::new();

    mol.add_atom(0, Element::C, 0);
    mol.add_atom(1, Element::C, 0);
    mol.add_atom(2, Element::C, 0);
    mol.add_atom(3, Element::C, 0);
    mol.add_atom(4, Element::C, 0);
    mol.add_atom(5, Element::C, 0);
    mol.add_atom(6, Element::C, 0);
    mol.add_atom(7, Element::C, 0);
    mol.add_atom(8, Element::C, 0);
    mol.add_atom(9, Element::N, 0);
    mol.add_atom(10, Element::C, 0);
    mol.add_atom(11, Element::C, 0);
    mol.add_atom(12, Element::C, 0);
    mol.add_atom(13, Element::C, 0);
    mol.add_atom(14, Element::H, 0);
    mol.add_atom(15, Element::H, 0);
    mol.add_atom(16, Element::H, 0);
    mol.add_atom(17, Element::H, 0);
    mol.add_atom(18, Element::H, 0);
    mol.add_atom(19, Element::H, 0);
    mol.add_atom(20, Element::H, 0);
    mol.add_atom(21, Element::H, 0);
    mol.add_atom(22, Element::H, 0);

    mol.add_bond(0, 0, 1, BondOrder::Aromatic);
    mol.add_bond(1, 1, 2, BondOrder::Aromatic);
    mol.add_bond(2, 2, 3, BondOrder::Aromatic);
    mol.add_bond(3, 3, 4, BondOrder::Aromatic);
    mol.add_bond(4, 4, 5, BondOrder::Aromatic);
    mol.add_bond(5, 5, 0, BondOrder::Aromatic);
    mol.add_bond(6, 4, 6, BondOrder::Aromatic);
    mol.add_bond(7, 6, 7, BondOrder::Aromatic);
    mol.add_bond(8, 7, 8, BondOrder::Aromatic);
    mol.add_bond(9, 8, 9, BondOrder::Aromatic);
    mol.add_bond(10, 9, 5, BondOrder::Aromatic);
    mol.add_bond(11, 9, 10, BondOrder::Aromatic);
    mol.add_bond(12, 10, 11, BondOrder::Aromatic);
    mol.add_bond(13, 11, 12, BondOrder::Aromatic);
    mol.add_bond(14, 12, 13, BondOrder::Aromatic);
    mol.add_bond(15, 13, 8, BondOrder::Aromatic);
    mol.add_bond(16, 0, 14, BondOrder::Single);
    mol.add_bond(17, 1, 15, BondOrder::Single);
    mol.add_bond(18, 2, 16, BondOrder::Single);
    mol.add_bond(19, 3, 17, BondOrder::Single);
    mol.add_bond(20, 6, 18, BondOrder::Single);
    mol.add_bond(21, 7, 19, BondOrder::Single);
    mol.add_bond(22, 10, 20, BondOrder::Single);
    mol.add_bond(23, 11, 21, BondOrder::Single);
    mol.add_bond(24, 12, 22, BondOrder::Single);

    mol
}

pub fn build_acridine_kekule() -> TestMolecule {
    let mut mol = TestMolecule::new();

    mol.add_atom(0, Element::C, 0);
    mol.add_atom(1, Element::C, 0);
    mol.add_atom(2, Element::C, 0);
    mol.add_atom(3, Element::C, 0);
    mol.add_atom(4, Element::C, 0);
    mol.add_atom(5, Element::C, 0);
    mol.add_atom(6, Element::C, 0);
    mol.add_atom(7, Element::C, 0);
    mol.add_atom(8, Element::C, 0);
    mol.add_atom(9, Element::N, 0);
    mol.add_atom(10, Element::C, 0);
    mol.add_atom(11, Element::C, 0);
    mol.add_atom(12, Element::C, 0);
    mol.add_atom(13, Element::C, 0);
    mol.add_atom(14, Element::H, 0);
    mol.add_atom(15, Element::H, 0);
    mol.add_atom(16, Element::H, 0);
    mol.add_atom(17, Element::H, 0);
    mol.add_atom(18, Element::H, 0);
    mol.add_atom(19, Element::H, 0);
    mol.add_atom(20, Element::H, 0);
    mol.add_atom(21, Element::H, 0);
    mol.add_atom(22, Element::H, 0);

    mol.add_bond(0, 0, 1, BondOrder::Double);
    mol.add_bond(1, 1, 2, BondOrder::Single);
    mol.add_bond(2, 2, 3, BondOrder::Double);
    mol.add_bond(3, 3, 4, BondOrder::Single);
    mol.add_bond(4, 4, 5, BondOrder::Double);
    mol.add_bond(5, 5, 0, BondOrder::Single);
    mol.add_bond(6, 4, 6, BondOrder::Single);
    mol.add_bond(7, 6, 7, BondOrder::Double);
    mol.add_bond(8, 7, 8, BondOrder::Single);
    mol.add_bond(9, 8, 9, BondOrder::Double);
    mol.add_bond(10, 9, 5, BondOrder::Single);
    mol.add_bond(11, 9, 10, BondOrder::Single);
    mol.add_bond(12, 10, 11, BondOrder::Double);
    mol.add_bond(13, 11, 12, BondOrder::Single);
    mol.add_bond(14, 12, 13, BondOrder::Double);
    mol.add_bond(15, 13, 8, BondOrder::Single);
    mol.add_bond(16, 0, 14, BondOrder::Single);
    mol.add_bond(17, 1, 15, BondOrder::Single);
    mol.add_bond(18, 2, 16, BondOrder::Single);
    mol.add_bond(19, 3, 17, BondOrder::Single);
    mol.add_bond(20, 6, 18, BondOrder::Single);
    mol.add_bond(21, 7, 19, BondOrder::Single);
    mol.add_bond(22, 10, 20, BondOrder::Single);
    mol.add_bond(23, 11, 21, BondOrder::Single);
    mol.add_bond(24, 12, 22, BondOrder::Single);

    mol
}

pub fn build_phosphinane() -> TestMolecule {
    let mut mol = TestMolecule::new();

    mol.add_atom(0, Element::P, 0);
    mol.add_atom(1, Element::C, 0);
    mol.add_atom(2, Element::C, 0);
    mol.add_atom(3, Element::C, 0);
    mol.add_atom(4, Element::C, 0);
    mol.add_atom(5, Element::C, 0);
    mol.add_atom(6, Element::H, 0);
    mol.add_atom(7, Element::H, 0);
    mol.add_atom(8, Element::H, 0);
    mol.add_atom(9, Element::H, 0);
    mol.add_atom(10, Element::H, 0);
    mol.add_atom(11, Element::H, 0);
    mol.add_atom(12, Element::H, 0);
    mol.add_atom(13, Element::H, 0);
    mol.add_atom(14, Element::H, 0);
    mol.add_atom(15, Element::H, 0);
    mol.add_atom(16, Element::H, 0);

    mol.add_bond(0, 0, 1, BondOrder::Single);
    mol.add_bond(1, 1, 2, BondOrder::Single);
    mol.add_bond(2, 2, 3, BondOrder::Single);
    mol.add_bond(3, 3, 4, BondOrder::Single);
    mol.add_bond(4, 4, 5, BondOrder::Single);
    mol.add_bond(5, 5, 0, BondOrder::Single);
    mol.add_bond(6, 0, 6, BondOrder::Single);
    mol.add_bond(7, 1, 7, BondOrder::Single);
    mol.add_bond(8, 1, 8, BondOrder::Single);
    mol.add_bond(9, 2, 9, BondOrder::Single);
    mol.add_bond(10, 2, 10, BondOrder::Single);
    mol.add_bond(11, 3, 11, BondOrder::Single);
    mol.add_bond(12, 3, 12, BondOrder::Single);
    mol.add_bond(13, 4, 13, BondOrder::Single);
    mol.add_bond(14, 4, 14, BondOrder::Single);
    mol.add_bond(15, 5, 15, BondOrder::Single);
    mol.add_bond(16, 5, 16, BondOrder::Single);

    mol
}

pub fn build_trinitrobenzene_aromatic() -> TestMolecule {
    let mut mol = TestMolecule::new();

    mol.add_atom(0, Element::C, 0);
    mol.add_atom(1, Element::C, 0);
    mol.add_atom(2, Element::C, 0);
    mol.add_atom(3, Element::C, 0);
    mol.add_atom(4, Element::C, 0);
    mol.add_atom(5, Element::C, 0);
    mol.add_atom(6, Element::N, 1);
    mol.add_atom(7, Element::O, -1);
    mol.add_atom(8, Element::O, 0);
    mol.add_atom(9, Element::N, 1);
    mol.add_atom(10, Element::O, -1);
    mol.add_atom(11, Element::O, 0);
    mol.add_atom(12, Element::N, 1);
    mol.add_atom(13, Element::O, -1);
    mol.add_atom(14, Element::O, 0);
    mol.add_atom(15, Element::H, 0);
    mol.add_atom(16, Element::H, 0);
    mol.add_atom(17, Element::H, 0);

    mol.add_bond(0, 0, 1, BondOrder::Aromatic);
    mol.add_bond(1, 1, 2, BondOrder::Aromatic);
    mol.add_bond(2, 2, 3, BondOrder::Aromatic);
    mol.add_bond(3, 3, 4, BondOrder::Aromatic);
    mol.add_bond(4, 4, 5, BondOrder::Aromatic);
    mol.add_bond(5, 5, 0, BondOrder::Aromatic);
    mol.add_bond(6, 0, 6, BondOrder::Single);
    mol.add_bond(7, 6, 7, BondOrder::Single);
    mol.add_bond(8, 6, 8, BondOrder::Double);
    mol.add_bond(9, 2, 9, BondOrder::Single);
    mol.add_bond(10, 9, 10, BondOrder::Single);
    mol.add_bond(11, 9, 11, BondOrder::Double);
    mol.add_bond(12, 4, 12, BondOrder::Single);
    mol.add_bond(13, 12, 13, BondOrder::Single);
    mol.add_bond(14, 12, 14, BondOrder::Double);
    mol.add_bond(15, 1, 15, BondOrder::Single);
    mol.add_bond(16, 3, 16, BondOrder::Single);
    mol.add_bond(17, 5, 17, BondOrder::Single);

    mol
}

pub fn build_trinitrobenzene_kekule() -> TestMolecule {
    let mut mol = TestMolecule::new();

    mol.add_atom(0, Element::C, 0);
    mol.add_atom(1, Element::C, 0);
    mol.add_atom(2, Element::C, 0);
    mol.add_atom(3, Element::C, 0);
    mol.add_atom(4, Element::C, 0);
    mol.add_atom(5, Element::C, 0);
    mol.add_atom(6, Element::N, 1);
    mol.add_atom(7, Element::O, -1);
    mol.add_atom(8, Element::O, 0);
    mol.add_atom(9, Element::N, 1);
    mol.add_atom(10, Element::O, -1);
    mol.add_atom(11, Element::O, 0);
    mol.add_atom(12, Element::N, 1);
    mol.add_atom(13, Element::O, -1);
    mol.add_atom(14, Element::O, 0);
    mol.add_atom(15, Element::H, 0);
    mol.add_atom(16, Element::H, 0);
    mol.add_atom(17, Element::H, 0);

    mol.add_bond(0, 0, 1, BondOrder::Double);
    mol.add_bond(1, 1, 2, BondOrder::Single);
    mol.add_bond(2, 2, 3, BondOrder::Double);
    mol.add_bond(3, 3, 4, BondOrder::Single);
    mol.add_bond(4, 4, 5, BondOrder::Double);
    mol.add_bond(5, 5, 0, BondOrder::Single);
    mol.add_bond(6, 0, 6, BondOrder::Single);
    mol.add_bond(7, 6, 7, BondOrder::Single);
    mol.add_bond(8, 6, 8, BondOrder::Double);
    mol.add_bond(9, 2, 9, BondOrder::Single);
    mol.add_bond(10, 9, 10, BondOrder::Single);
    mol.add_bond(11, 9, 11, BondOrder::Double);
    mol.add_bond(12, 4, 12, BondOrder::Single);
    mol.add_bond(13, 12, 13, BondOrder::Single);
    mol.add_bond(14, 12, 14, BondOrder::Double);
    mol.add_bond(15, 1, 15, BondOrder::Single);
    mol.add_bond(16, 3, 16, BondOrder::Single);
    mol.add_bond(17, 5, 17, BondOrder::Single);

    mol
}

pub fn build_tetramethylthiourea() -> TestMolecule {
    let mut mol = TestMolecule::new();

    mol.add_atom(0, Element::C, 0);
    mol.add_atom(1, Element::S, 0);
    mol.add_atom(2, Element::N, 0);
    mol.add_atom(3, Element::N, 0);
    mol.add_atom(4, Element::C, 0);
    mol.add_atom(5, Element::C, 0);
    mol.add_atom(6, Element::C, 0);
    mol.add_atom(7, Element::C, 0);
    mol.add_atom(8, Element::H, 0);
    mol.add_atom(9, Element::H, 0);
    mol.add_atom(10, Element::H, 0);
    mol.add_atom(11, Element::H, 0);
    mol.add_atom(12, Element::H, 0);
    mol.add_atom(13, Element::H, 0);
    mol.add_atom(14, Element::H, 0);
    mol.add_atom(15, Element::H, 0);
    mol.add_atom(16, Element::H, 0);
    mol.add_atom(17, Element::H, 0);
    mol.add_atom(18, Element::H, 0);
    mol.add_atom(19, Element::H, 0);

    mol.add_bond(0, 0, 1, BondOrder::Double);
    mol.add_bond(1, 0, 2, BondOrder::Single);
    mol.add_bond(2, 0, 3, BondOrder::Single);
    mol.add_bond(3, 2, 4, BondOrder::Single);
    mol.add_bond(4, 2, 5, BondOrder::Single);
    mol.add_bond(5, 3, 6, BondOrder::Single);
    mol.add_bond(6, 3, 7, BondOrder::Single);
    mol.add_bond(7, 4, 8, BondOrder::Single);
    mol.add_bond(8, 4, 9, BondOrder::Single);
    mol.add_bond(9, 4, 10, BondOrder::Single);
    mol.add_bond(10, 5, 11, BondOrder::Single);
    mol.add_bond(11, 5, 12, BondOrder::Single);
    mol.add_bond(12, 5, 13, BondOrder::Single);
    mol.add_bond(13, 6, 14, BondOrder::Single);
    mol.add_bond(14, 6, 15, BondOrder::Single);
    mol.add_bond(15, 6, 16, BondOrder::Single);
    mol.add_bond(16, 7, 17, BondOrder::Single);
    mol.add_bond(17, 7, 18, BondOrder::Single);
    mol.add_bond(18, 7, 19, BondOrder::Single);

    mol
}

pub fn build_dimethyl_sulfoxide() -> TestMolecule {
    let mut mol = TestMolecule::new();

    mol.add_atom(0, Element::S, 1);
    mol.add_atom(1, Element::O, -1);
    mol.add_atom(2, Element::C, 0);
    mol.add_atom(3, Element::C, 0);
    mol.add_atom(4, Element::H, 0);
    mol.add_atom(5, Element::H, 0);
    mol.add_atom(6, Element::H, 0);
    mol.add_atom(7, Element::H, 0);
    mol.add_atom(8, Element::H, 0);
    mol.add_atom(9, Element::H, 0);

    mol.add_bond(0, 0, 1, BondOrder::Single);
    mol.add_bond(1, 0, 2, BondOrder::Single);
    mol.add_bond(2, 0, 3, BondOrder::Single);
    mol.add_bond(3, 2, 4, BondOrder::Single);
    mol.add_bond(4, 2, 5, BondOrder::Single);
    mol.add_bond(5, 2, 6, BondOrder::Single);
    mol.add_bond(6, 3, 7, BondOrder::Single);
    mol.add_bond(7, 3, 8, BondOrder::Single);
    mol.add_bond(8, 3, 9, BondOrder::Single);

    mol
}

pub fn build_methanesulfonamide() -> TestMolecule {
    let mut mol = TestMolecule::new();

    mol.add_atom(0, Element::S, 0);
    mol.add_atom(1, Element::O, 0);
    mol.add_atom(2, Element::O, 0);
    mol.add_atom(3, Element::N, 0);
    mol.add_atom(4, Element::C, 0);
    mol.add_atom(5, Element::H, 0);
    mol.add_atom(6, Element::H, 0);
    mol.add_atom(7, Element::H, 0);
    mol.add_atom(8, Element::H, 0);
    mol.add_atom(9, Element::H, 0);

    mol.add_bond(0, 0, 1, BondOrder::Double);
    mol.add_bond(1, 0, 2, BondOrder::Double);
    mol.add_bond(2, 0, 3, BondOrder::Single);
    mol.add_bond(3, 0, 4, BondOrder::Single);
    mol.add_bond(4, 3, 5, BondOrder::Single);
    mol.add_bond(5, 3, 6, BondOrder::Single);
    mol.add_bond(6, 4, 7, BondOrder::Single);
    mol.add_bond(7, 4, 8, BondOrder::Single);
    mol.add_bond(8, 4, 9, BondOrder::Single);

    mol
}
