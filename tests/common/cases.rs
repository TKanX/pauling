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
    phenylalanine_zwitterion_aromatic => {
        title: "Phenylalanine Zwitterion (Aromatic Input)",
        tags: &["biomolecule", "zwitterion", "aromatic"],
        molecule: molecules::build_phenylalanine_zwitterion_aromatic,
        systems: [
            { atoms: || vec![2, 3, 4], bonds: || vec![2, 3] },
            { atoms: || vec![6, 7, 8, 9, 10, 11], bonds: || vec![6, 7, 8, 9, 10, 11] },
        ],
    },
    phenylalanine_zwitterion_kekule => {
        title: "Phenylalanine Zwitterion (Kekulé Input)",
        tags: &["biomolecule", "zwitterion", "aromatic", "kekule"],
        molecule: molecules::build_phenylalanine_zwitterion_kekule,
        systems: [
            { atoms: || vec![2, 3, 4], bonds: || vec![2, 3] },
            { atoms: || vec![6, 7, 8, 9, 10, 11], bonds: || vec![6, 7, 8, 9, 10, 11] },
        ],
    },
    tyrosine_zwitterion_aromatic => {
        title: "Tyrosine Zwitterion (Aromatic Input)",
        tags: &["biomolecule", "zwitterion", "aromatic", "carboxylate", "phenol"],
        molecule: molecules::build_tyrosine_zwitterion_aromatic,
        systems: [
            { atoms: || vec![2, 3, 4], bonds: || vec![2, 3] },
            { atoms: || vec![6, 7, 8, 9, 10, 11, 12], bonds: || vec![6, 7, 8, 9, 10, 11, 12] },
        ],
    },
    tyrosine_zwitterion_kekule => {
        title: "Tyrosine Zwitterion (Kekulé Input)",
        tags: &["biomolecule", "zwitterion", "aromatic", "carboxylate", "phenol"],
        molecule: molecules::build_tyrosine_zwitterion_kekule,
        systems: [
            { atoms: || vec![2, 3, 4], bonds: || vec![2, 3] },
            { atoms: || vec![6, 7, 8, 9, 10, 11, 12], bonds: || vec![6, 7, 8, 9, 10, 11, 12] },
        ],
    },
    histidine_zwitterion_aromatic => {
        title: "Histidine Zwitterion (Aromatic Input)",
        tags: &["biomolecule", "zwitterion", "aromatic", "heteroaromatic"],
        molecule: molecules::build_histidine_zwitterion_aromatic,
        systems: [
            { atoms: || vec![2, 3, 4], bonds: || vec![2, 3] },
            { atoms: || vec![6, 7, 8, 9, 10], bonds: || vec![6, 7, 8, 9, 10] },
        ],
    },
    histidine_zwitterion_kekule => {
        title: "Histidine Zwitterion (Kekulé Input)",
        tags: &["biomolecule", "zwitterion", "aromatic", "heteroaromatic", "kekule"],
        molecule: molecules::build_histidine_zwitterion_kekule,
        systems: [
            { atoms: || vec![2, 3, 4], bonds: || vec![2, 3] },
            { atoms: || vec![6, 7, 8, 9, 10], bonds: || vec![6, 7, 8, 9, 10] },
        ],
    },
    tryptophan_zwitterion_aromatic => {
        title: "Tryptophan Zwitterion (Aromatic Input)",
        tags: &["biomolecule", "zwitterion", "aromatic", "indole"],
        molecule: molecules::build_tryptophan_zwitterion_aromatic,
        systems: [
            { atoms: || vec![2, 3, 4], bonds: || vec![2, 3] },
            { atoms: || vec![6, 7, 8, 9, 10, 11, 12, 13, 14], bonds: || vec![6, 7, 8, 9, 10, 11, 12, 13, 14, 15] },
        ],
    },
    tryptophan_zwitterion_kekule => {
        title: "Tryptophan Zwitterion (Kekulé Input)",
        tags: &["biomolecule", "zwitterion", "kekule", "indole"],
        molecule: molecules::build_tryptophan_zwitterion_kekule,
        systems: [
            { atoms: || vec![2, 3, 4], bonds: || vec![2, 3] },
            { atoms: || vec![6, 7, 8, 9, 10, 11, 12, 13, 14], bonds: || vec![6, 7, 8, 9, 10, 11, 12, 13, 14, 15] },
        ],
    },
    uracil_aromatic_input => {
        title: "Uracil (Aromatic Input)",
        tags: &["biomolecule", "nucleobase", "heteroaromatic", "aromatic-input"],
        molecule: molecules::build_uracil_aromatic,
        systems: [
            { atoms: || vec![0, 1, 2, 3, 4, 5, 6, 7], bonds: || vec![0, 1, 2, 3, 4, 5, 6, 7] },
        ],
    },
    uracil_kekule_input => {
        title: "Uracil (Kekulé Input)",
        tags: &["biomolecule", "nucleobase", "heteroaromatic", "kekule-input"],
        molecule: molecules::build_uracil_kekule,
        systems: [
            { atoms: || vec![0, 1, 2, 3, 4, 5, 6, 7], bonds: || vec![0, 1, 2, 3, 4, 5, 6, 7] },
        ],
    },
    thymine_aromatic_input => {
        title: "Thymine (Aromatic Input)",
        tags: &["biomolecule", "nucleobase", "heterocycle"],
        molecule: molecules::build_thymine_aromatic,
        systems: [
            { atoms: || vec![0, 1, 2, 3, 4, 5, 6, 7], bonds: || vec![0, 1, 2, 3, 4, 5, 6, 7] },
        ],
    },
    thymine_kekule_input => {
        title: "Thymine (Kekulé Input)",
        tags: &["biomolecule", "nucleobase", "heterocycle", "kekule"],
        molecule: molecules::build_thymine_kekule,
        systems: [
            { atoms: || vec![0, 1, 2, 3, 4, 5, 6, 7], bonds: || vec![0, 1, 2, 3, 4, 5, 6, 7] },
        ],
    },
    cytosine_aromatic_input => {
        title: "Cytosine (Aromatic Input)",
        tags: &["biomolecule", "nucleobase", "heterocycle"],
        molecule: molecules::build_cytosine_aromatic,
        systems: [
            { atoms: || vec![0, 1, 2, 3, 4, 5, 6, 7], bonds: || vec![0, 1, 2, 3, 4, 5, 6, 7] },
        ],
    },
    cytosine_kekule_input => {
        title: "Cytosine (Kekulé Input)",
        tags: &["biomolecule", "nucleobase", "heterocycle", "kekule"],
        molecule: molecules::build_cytosine_kekule,
        systems: [
            { atoms: || vec![0, 1, 2, 3, 4, 5, 6, 7], bonds: || vec![0, 1, 2, 3, 4, 5, 6, 7] },
        ],
    },
    adenine_aromatic_input => {
        title: "Adenine (Aromatic Input)",
        tags: &["biomolecule", "nucleobase", "heterocycle", "aromatic"],
        molecule: molecules::build_adenine_aromatic,
        systems: [
            { atoms: || vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9], bonds: || vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10] },
        ],
    },
    adenine_kekule_input => {
        title: "Adenine (Kekulé Input)",
        tags: &["biomolecule", "nucleobase", "heterocycle", "kekule"],
        molecule: molecules::build_adenine_kekule,
        systems: [
            { atoms: || vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9], bonds: || vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10] },
        ],
    },
    guanine_aromatic_input => {
        title: "Guanine (Aromatic Input)",
        tags: &["biomolecule", "nucleobase", "heterocycle", "fused"],
        molecule: molecules::build_guanine_aromatic,
        systems: [
            { atoms: || vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10], bonds: || vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11] },
        ],
    },
    guanine_kekule_input => {
        title: "Guanine (Kekulé Input)",
        tags: &["biomolecule", "nucleobase", "heterocycle", "fused", "kekule"],
        molecule: molecules::build_guanine_kekule,
        systems: [
            { atoms: || vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10], bonds: || vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11] },
        ],
    },
    deoxyadenosine_aromatic_input => {
        title: "Deoxyadenosine (Aromatic Input)",
        tags: &["biomolecule", "nucleoside", "aromatic", "purine"],
        molecule: molecules::build_deoxyadenosine_aromatic,
        systems: [
            { atoms: || vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9], bonds: || vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10] },
        ],
    },
    deoxyadenosine_kekule_input => {
        title: "Deoxyadenosine (Kekulé Input)",
        tags: &["biomolecule", "nucleoside", "kekule", "purine"],
        molecule: molecules::build_deoxyadenosine_kekule,
        systems: [
            { atoms: || vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9], bonds: || vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10] },
        ],
    },
    dinucleotide_backbone => {
        title: "Dinucleotide Backbone Fragment",
        tags: &["biomolecule", "phosphate", "backbone"],
        molecule: molecules::build_dinucleotide_backbone,
        systems: [
            { atoms: || vec![2, 3, 4], bonds: || vec![4, 5] },
        ],
    },
    acbuol_isolated_lactone => {
        title: "ACBUOL - Isolated Lactone Resonance",
        tags: &["lactone", "polycyclic", "isolated-resonance"],
        molecule: molecules::build_acbuol,
        systems: [
            { atoms: || vec![4, 5, 6], bonds: || vec![8, 9] },
        ],
    },
    afurpo10_cage_ether => {
        title: "AFURPO10 - Saturated Cage Ether",
        tags: &["saturated", "polycyclic", "negative-control"],
        molecule: molecules::build_afurpo10,
        systems: [],
    },
    adamantane_saturated_polycycle => {
        title: "Adamantane (Saturated Polycyclic Alkane)",
        tags: &["saturated", "polycyclic", "negative-case", "alkane"],
        molecule: molecules::build_adamantane,
        systems: [],
    },
    decalin_saturated_fused_bicycle => {
        title: "Decalin (Saturated Fused Bicyclic Alkane)",
        tags: &["saturated", "fused", "negative-case", "alkane"],
        molecule: molecules::build_decalin,
        systems: [],
    },
    abaxes_dithiacyclohexene => {
        title: "ABAXES - 1,4-dithiacyclohex-2-ene",
        tags: &["heterocycle", "sulfur", "partial-conjugation"],
        molecule: molecules::build_abaxes,
        systems: [
            { atoms: || vec![0, 1, 2, 3], bonds: || vec![0, 1, 2] },
        ],
    },
    methyl_pyridine_aromatic_input => {
        title: "4-Methylpyridine (Aromatic Input)",
        tags: &["aromatic", "heterocycle"],
        molecule: molecules::build_4_methylpyridine_aromatic,
        systems: [
            { atoms: || vec![0, 1, 2, 3, 4, 5], bonds: || vec![0, 1, 2, 3, 4, 5] },
        ],
    },
    methyl_pyridine_kekule_input => {
        title: "4-Methylpyridine (Kekulé Input)",
        tags: &["aromatic", "heterocycle", "kekule"],
        molecule: molecules::build_4_methylpyridine_kekule,
        systems: [
            { atoms: || vec![0, 1, 2, 3, 4, 5], bonds: || vec![0, 1, 2, 3, 4, 5] },
        ],
    },
    acridine_aromatic_input => {
        title: "Acridine (Aromatic Input)",
        tags: &["aromatic", "heterocycle", "fused"],
        molecule: molecules::build_acridine_aromatic,
        systems: [
            { atoms: || (0..14).collect(), bonds: || (0..16).collect() },
        ],
    },
    acridine_kekule_input => {
        title: "Acridine (Kekulé Input)",
        tags: &["aromatic", "heterocycle", "fused", "kekule"],
        molecule: molecules::build_acridine_kekule,
        systems: [
            { atoms: || (0..14).collect(), bonds: || (0..16).collect() },
        ],
    },
    phosphinane_saturated_heterocycle => {
        title: "Phosphinane (Saturated Heterocycle)",
        tags: &["saturated", "heterocycle", "phosphorus", "negative-case"],
        molecule: molecules::build_phosphinane,
        systems: [],
    },
    trinitrobenzene_aromatic_input => {
        title: "1,3,5-Trinitrobenzene (Aromatic Input)",
        tags: &["aromatic", "nitro", "cross-conjugation"],
        molecule: molecules::build_trinitrobenzene_aromatic,
        systems: [
            { atoms: || (0..15).collect(), bonds: || (0..15).collect() },
        ],
    },
    trinitrobenzene_kekule_input => {
        title: "1,3,5-Trinitrobenzene (Kekulé Input)",
        tags: &["aromatic", "nitro", "cross-conjugation", "kekule"],
        molecule: molecules::build_trinitrobenzene_kekule,
        systems: [
            { atoms: || (0..15).collect(), bonds: || (0..15).collect() },
        ],
    },
}
