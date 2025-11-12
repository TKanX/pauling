//! Core primitives reused across every perception stage.
//!
//! The submodules provide canonical definitions for atoms and bonds, including
//! the identifiers and enumerations that appear throughout the public API.

/// Atom-centric primitives such as [`AtomId`](crate::AtomId) and [`Element`](crate::Element).
pub mod atom;
/// Bond-centric primitives such as [`BondId`](crate::BondId) and [`BondOrder`](crate::BondOrder).
pub mod bond;
