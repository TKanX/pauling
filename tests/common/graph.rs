use pauling::traits::{AtomView, BondView, MoleculeGraph};
use pauling::{AtomId, BondId, BondOrder, Element};

#[derive(Clone)]
pub struct TestAtom {
    pub id: AtomId,
    pub element: Element,
    pub formal_charge: i8,
}

impl AtomView for TestAtom {
    fn id(&self) -> AtomId {
        self.id
    }

    fn element(&self) -> Element {
        self.element
    }

    fn formal_charge(&self) -> i8 {
        self.formal_charge
    }
}

#[derive(Clone)]
pub struct TestBond {
    pub id: BondId,
    pub order: BondOrder,
    pub start: AtomId,
    pub end: AtomId,
}

impl BondView for TestBond {
    fn id(&self) -> BondId {
        self.id
    }

    fn order(&self) -> BondOrder {
        self.order
    }

    fn start_atom_id(&self) -> AtomId {
        self.start
    }

    fn end_atom_id(&self) -> AtomId {
        self.end
    }
}

#[derive(Clone, Default)]
pub struct TestMolecule {
    pub atoms: Vec<TestAtom>,
    pub bonds: Vec<TestBond>,
}

impl TestMolecule {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add_atom(&mut self, id: AtomId, element: Element, formal_charge: i8) {
        self.atoms.push(TestAtom {
            id,
            element,
            formal_charge,
        });
    }

    pub fn add_bond(&mut self, id: BondId, start: AtomId, end: AtomId, order: BondOrder) {
        self.bonds.push(TestBond {
            id,
            order,
            start,
            end,
        });
    }
}

impl MoleculeGraph for TestMolecule {
    type Atom = TestAtom;
    type Bond = TestBond;

    fn atoms(&self) -> impl Iterator<Item = &Self::Atom> {
        self.atoms.iter()
    }

    fn bonds(&self) -> impl Iterator<Item = &Self::Bond> {
        self.bonds.iter()
    }
}
