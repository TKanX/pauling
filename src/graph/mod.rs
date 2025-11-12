//! Graph abstraction layer used by the perception pipeline.
//!
//! Consumers can implement the traits in this module to expose their molecular
//! data structures without copying into the built-in `Molecule` container.

/// Traits that describe atoms, bonds, and the containing graph.
pub mod traits;
