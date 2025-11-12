//! Error types emitted while running the perception pipeline.

use crate::core::atom::AtomId;
use thiserror::Error;

/// Error returned when a perception stage cannot complete successfully.
#[derive(Debug, Error)]
pub enum PerceptionError {
    /// The input graph references an atom identifier that was not provided.
    #[error("graph integrity error: bond references non-existent atom ID {0}")]
    InconsistentGraph(AtomId),

    /// Multiple edges connect the same pair of atoms in the source graph.
    #[error("invalid graph topology: duplicate bond detected between atoms {start} and {end}")]
    DuplicateBond { start: AtomId, end: AtomId },

    /// Kekulization exhausted its attempt budget without finding a valid pattern.
    #[error("kekulization failed for an aromatic component after {0} attempts")]
    KekulizationFailed(usize),

    /// The ring perception stage reported a failure.
    #[error("ring perception failed: {0}")]
    RingPerceptionFailed(String),
}
