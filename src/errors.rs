use crate::core::atom::AtomId;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum PerceptionError {
    #[error("graph integrity error: bond references non-existent atom ID {0}")]
    InconsistentGraph(AtomId),

    #[error("kekulization failed for an aromatic component after {0} attempts")]
    KekulizationFailed(usize),

    #[error("ring perception failed: {0}")]
    RingPerceptionFailed(String),
}
