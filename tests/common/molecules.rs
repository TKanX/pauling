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
