#![allow(non_upper_case_globals)]

mod graph;
mod molecules;

pub use graph::TestMolecule;
use pauling::{AtomId, BondId, ResonanceSystem};

pub struct ResonanceSystemExpectation {
    pub atoms: fn() -> Vec<AtomId>,
    pub bonds: fn() -> Vec<BondId>,
}

impl ResonanceSystemExpectation {
    pub fn to_system(&self) -> ResonanceSystem {
        let atoms = (self.atoms)();
        let bonds = (self.bonds)();
        ResonanceSystem::new(atoms, bonds)
    }
}

pub struct ResonanceCase {
    pub slug: &'static str,
    pub title: &'static str,
    pub tags: &'static [&'static str],
    pub build: fn() -> TestMolecule,
    pub expected: &'static [ResonanceSystemExpectation],
}

impl ResonanceCase {
    pub fn expected_systems(&self) -> Vec<ResonanceSystem> {
        self.expected
            .iter()
            .map(ResonanceSystemExpectation::to_system)
            .collect()
    }
}

macro_rules! resonance_cases {
    ( $( $name:ident => {
            title: $title:expr,
            tags: $tags:expr,
            molecule: $builder:path,
            systems: [ $( { atoms: $atoms:expr, bonds: $bonds:expr } ),* $(,)? ],
        }
    ),* $(,)? ) => {
        $(
            pub const $name: ResonanceCase = ResonanceCase {
                slug: stringify!($name),
                title: $title,
                tags: $tags,
                build: $builder,
                expected: &[
                    $(
                        ResonanceSystemExpectation {
                            atoms: $atoms,
                            bonds: $bonds,
                        }
                    ),*
                ],
            };
        )*

        pub const ALL_CASES: &[ResonanceCase] = &[
            $( $name ),*
        ];

        macro_rules! for_each_resonance_case {
            ($macro:ident) => {
                $(
                    $macro!($name);
                )*
            };
        }
    };
}

resonance_cases! {
    glycine_zwitterion => {
        title: "Glycine Zwitterion",
        tags: &["biomolecule", "zwitterion", "carboxylate"],
        molecule: molecules::build_glycine_zwitterion,
        systems: [
            { atoms: || vec![2, 3, 4], bonds: || vec![2, 3] },
        ],
    },
    alanine_zwitterion => {
        title: "Alanine Zwitterion",
        tags: &["biomolecule", "zwitterion", "carboxylate"],
        molecule: molecules::build_alanine_zwitterion,
        systems: [
            { atoms: || vec![2, 3, 4], bonds: || vec![2, 3] },
        ],
    },
    valine_zwitterion => {
        title: "Valine Zwitterion",
        tags: &["biomolecule", "zwitterion", "carboxylate"],
        molecule: molecules::build_valine_zwitterion,
        systems: [
            { atoms: || vec![2, 3, 4], bonds: || vec![2, 3] },
        ],
    },
    leucine_zwitterion => {
        title: "Leucine Zwitterion",
        tags: &["biomolecule", "zwitterion", "carboxylate"],
        molecule: molecules::build_leucine_zwitterion,
        systems: [
            { atoms: || vec![2, 3, 4], bonds: || vec![2, 3] },
        ],
    },
    isoleucine_zwitterion => {
        title: "Isoleucine Zwitterion",
        tags: &["biomolecule", "zwitterion", "carboxylate"],
        molecule: molecules::build_isoleucine_zwitterion,
        systems: [
            { atoms: || vec![2, 3, 4], bonds: || vec![2, 3] },
        ],
    },
    proline_zwitterion => {
        title: "Proline Zwitterion",
        tags: &["biomolecule", "zwitterion", "carboxylate", "cyclic"],
        molecule: molecules::build_proline_zwitterion,
        systems: [
            { atoms: || vec![5, 6, 7], bonds: || vec![6, 7] },
        ],
    },
    serine_zwitterion => {
        title: "Serine Zwitterion",
        tags: &["biomolecule", "zwitterion", "carboxylate"],
        molecule: molecules::build_serine_zwitterion,
        systems: [
            { atoms: || vec![2, 3, 4], bonds: || vec![2, 3] },
        ],
    },
    threonine_zwitterion => {
        title: "Threonine Zwitterion",
        tags: &["biomolecule", "zwitterion", "carboxylate"],
        molecule: molecules::build_threonine_zwitterion,
        systems: [
            { atoms: || vec![2, 3, 4], bonds: || vec![2, 3] },
        ],
    },
    cysteine_zwitterion => {
        title: "Cysteine Zwitterion",
        tags: &["biomolecule", "zwitterion", "carboxylate", "sulfur"],
        molecule: molecules::build_cysteine_zwitterion,
        systems: [
            { atoms: || vec![2, 3, 4], bonds: || vec![2, 3] },
        ],
    },
    methionine_zwitterion => {
        title: "Methionine Zwitterion",
        tags: &["biomolecule", "zwitterion", "sulfur"],
        molecule: molecules::build_methionine_zwitterion,
        systems: [
            { atoms: || vec![2, 3, 4], bonds: || vec![2, 3] },
        ],
    },
    aspartate_zwitterion => {
        title: "Aspartate Zwitterion",
        tags: &["biomolecule", "zwitterion", "carboxylate"],
        molecule: molecules::build_aspartate_zwitterion,
        systems: [
            { atoms: || vec![2, 3, 4], bonds: || vec![2, 3] },
            { atoms: || vec![6, 7, 8], bonds: || vec![6, 7] },
        ],
    },
    asparagine_zwitterion => {
        title: "Asparagine Zwitterion",
        tags: &["biomolecule", "zwitterion", "amide", "carboxylate", "disconnected"],
        molecule: molecules::build_asparagine_zwitterion,
        systems: [
            { atoms: || vec![2, 3, 4], bonds: || vec![2, 3] },
            { atoms: || vec![6, 7, 8], bonds: || vec![6, 7] },
        ],
    },
    glutamate_zwitterion => {
        title: "Glutamate Zwitterion",
        tags: &["biomolecule", "zwitterion", "carboxylate"],
        molecule: molecules::build_glutamate_zwitterion,
        systems: [
            { atoms: || vec![2, 3, 4], bonds: || vec![2, 3] },
            { atoms: || vec![7, 8, 9], bonds: || vec![7, 8] },
        ],
    },
    glutamine_zwitterion => {
        title: "Glutamine Zwitterion",
        tags: &["biomolecule", "zwitterion", "carboxylate", "amide"],
        molecule: molecules::build_glutamine_zwitterion,
        systems: [
            { atoms: || vec![2, 3, 4], bonds: || vec![2, 3] },
            { atoms: || vec![7, 8, 9], bonds: || vec![7, 8] },
        ],
    },
    lysine_zwitterion => {
        title: "Lysine Zwitterion",
        tags: &["biomolecule", "zwitterion", "carboxylate", "amine"],
        molecule: molecules::build_lysine_zwitterion,
        systems: [
            { atoms: || vec![2, 3, 4], bonds: || vec![2, 3] },
        ],
    },
}
