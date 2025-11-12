# Pauling

**Pauling** is a molecular perception library focused on high-fidelity resonance detection. It transforms generic molecular graphs into chemically rich descriptions by identifying rings, aromatic systems, Kekulé patterns, atomic states, and ultimately distinct resonance systems. The library is written in **Rust** to deliver predictable performance, strong correctness guarantees, and seamless integration into larger cheminformatics pipelines.

The goal of Pauling is to serve as a dependable building block for applications that require accurate conjugation and aromaticity perception, such as reactivity prediction, force-field assignment, and advanced molecular visualization.

> _The theory of chemical resonance was first developed by Linus **Pauling** in the 1930s._

## Features

- **Multi-Stage Perception Pipeline**: Ring detection, aromaticity analysis, Kekulé assignment, and atomic state inference are performed in sequence before resonance systems are reported.
- **Graph Agnostic**: Operates on any data structure that implements the `MoleculeGraph` trait, allowing zero-copy integration with existing tooling.
- **High-Quality Chemical Heuristics**: Implements Hückel aromaticity tests, conjugation heuristics for heteroatoms, and lone-pair promotion rules for amides and similar motifs.
- **Robust Error Reporting**: Provides descriptive errors when graph integrity or perception steps fail, simplifying debugging in downstream applications.

## Getting Started

Pauling is currently distributed as a library crate. Add it to your `Cargo.toml` and bring the traits and helpers you need into scope.

```toml
[dependencies]
pauling = "0.1.0"
```

You can then compile the documentation locally:

```sh
cargo doc --open
```

## Example: Benzene Resonance Systems

The example below constructs benzene using the bundled `Molecule` type and runs the full perception pipeline to recover its aromatic resonance system.

```rust
use pauling::{find_resonance_systems, BondOrder, Element, Molecule};

fn main() -> Result<(), pauling::PerceptionError> {
    let mut molecule = Molecule::new();

    let carbons: Vec<_> = (0..6).map(|_| molecule.add_atom(Element::C, 0)).collect();
    let hydrogens: Vec<_> = (0..6).map(|_| molecule.add_atom(Element::H, 0)).collect();

    let orders = [
        BondOrder::Double,
        BondOrder::Single,
        BondOrder::Double,
        BondOrder::Single,
        BondOrder::Double,
        BondOrder::Single,
    ];

    for i in 0..6 {
        let next = (i + 1) % 6;
        molecule.add_bond(carbons[i], carbons[next], orders[i]).unwrap();
        molecule.add_bond(carbons[i], hydrogens[i], BondOrder::Single).unwrap();
    }

    let systems = find_resonance_systems(&molecule)?;

    assert_eq!(systems.len(), 1);
    assert_eq!(systems[0].atoms, carbons);
    assert_eq!(systems[0].bonds.len(), 6);

    Ok(())
}
```

> **Note**: This is a simplified example. For more advanced usage, including custom graph implementations and error handling, please refer to the [API Documentation](https://docs.rs/pauling).

## Documentation

- [API Documentation](https://docs.rs/pauling) – comprehensive reference for public types and functions.
- [Architecture Overview](ARCHITECTURE.md) – detailed explanation of the internal design and algorithms used in Pauling.

## Tech Stack

- **Core Language**: Rust
- **Testing & Docs**: `cargo test`, `cargo doc`
- **Error Handling**: `thiserror`

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
