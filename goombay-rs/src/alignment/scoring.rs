pub trait Scoring {
    fn get_match_score(&self) -> usize;
    fn get_mismatch_score(&self) -> usize;
    fn get_gap_score(&self) -> usize;
    fn get_extended_gap_score(&self) -> usize;
    fn get_transpose_score(&self) -> usize;
}

#[derive(Clone)]
pub struct LevenshteinScoring {
    pub substitution: usize,
    pub gap: usize,
}

impl Scoring for LevenshteinScoring {
    fn get_match_score(&self) -> usize {
        0_usize
    }
    fn get_mismatch_score(&self) -> usize {
        self.substitution
    }
    fn get_gap_score(&self) -> usize {
        self.gap
    }
    fn get_extended_gap_score(&self) -> usize {
        0_usize
    }
    fn get_transpose_score(&self) -> usize {
        0_usize
    }
}

#[derive(Clone)]
pub struct GeneralScoring {
    pub identity: usize,
    pub mismatch: usize,
    pub gap: usize,
}

impl Scoring for GeneralScoring {
    fn get_match_score(&self) -> usize {
        self.identity
    }
    fn get_mismatch_score(&self) -> usize {
        self.mismatch
    }
    fn get_gap_score(&self) -> usize {
        self.gap
    }
    fn get_extended_gap_score(&self) -> usize {
        0_usize
    }
    fn get_transpose_score(&self) -> usize {
        0_usize
    }
}

#[derive(Clone)]
pub struct TransposeScoring {
    pub identity: usize,
    pub mismatch: usize,
    pub gap: usize,
    pub transpose: usize,
}

impl Scoring for TransposeScoring {
    fn get_match_score(&self) -> usize {
        self.identity
    }
    fn get_mismatch_score(&self) -> usize {
        self.mismatch
    }
    fn get_gap_score(&self) -> usize {
        self.gap
    }
    fn get_extended_gap_score(&self) -> usize {
        0_usize
    }
    fn get_transpose_score(&self) -> usize {
        self.transpose
    }
}

#[derive(Clone)]
pub struct ExtendedGapScoring {
    pub identity: usize,
    pub mismatch: usize,
    pub gap: usize,
    pub extended_gap: usize,
}

impl Scoring for ExtendedGapScoring {
    fn get_match_score(&self) -> usize {
        self.identity
    }
    fn get_mismatch_score(&self) -> usize {
        self.mismatch
    }
    fn get_gap_score(&self) -> usize {
        self.gap
    }
    fn get_extended_gap_score(&self) -> usize {
        self.extended_gap
    }
    fn get_transpose_score(&self) -> usize {
        0_usize
    }
}
