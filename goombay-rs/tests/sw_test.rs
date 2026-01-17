use goombay_rs::align::{LocalAlignmentMatrices, SmithWaterman};
use goombay_rs::scoring::GeneralScoring;

#[test]
fn test_identical_sequences() {
    let sw = SmithWaterman::compute("ACTG", "ACTG");
    let aligned = sw.align();
    let sim = sw.similarity();

    assert_eq!(aligned[0], "ACTG\nACTG");
    assert_eq!(sim, (4 * sw.identity) as i32); 
}

#[test]
fn test_local_match() {
    // "ACTG" is a local match within "AAACTGAA" and "TTTACTGTT"
    let sw = SmithWaterman::compute("AAACTGAA", "TTTACTGTT");
    let aligned = sw.align();
    let sim = sw.similarity();

    // The best local alignment should be ACTG vs ACTG
    assert_eq!(aligned[0], "ACTG\nACTG");
    assert_eq!(sim, (4 * sw.identity) as i32);
}

#[test]
fn test_no_similarity() {
    let sw = SmithWaterman::compute("AAAA", "TTTT");
    let aligned = sw.align();
    let sim = sw.similarity();

    // With mismatch score 1 and identity 2, A vs T gives 0 if negative.
    // In our implementation, a single mismatch would be -1, which is maxed at 0.
    assert_eq!(sim, 0);
    assert!(aligned.is_empty() || aligned[0].trim().is_empty() || aligned[0] == "\n");
}

#[test]
fn test_different_length_local() {
    let sw = SmithWaterman::compute("GCATGCG", "GATTACA");
    let aligned = sw.align();
    // G C A T G C G
    // G A T T A C A
    // G vs G = 2
    // A vs A = 2
    // T vs T = 2
    // T vs T = 2
    // G-A Match? 
    // Let's see what a standard SW would give
    // G C A T G C G
    // G A T T A C A
    // G..C..A..T..G..C..G
    // G..A..T..T..A..C..A
    // G   -   A   T
    // G   A   -   T  ?
    
    assert!(sw.similarity() > 0);
    assert!(!aligned[0].is_empty());
}

#[test]
fn test_scoring_parameters_sw() {
    let custom_scores = GeneralScoring {
        identity: 3,
        mismatch: 3,
        gap: 2,
    };
    let custom_sw = SmithWaterman::set_scores(&custom_scores);
    let sw_alignment = custom_sw.calculate_matrix("ACGT", "ACG");
    assert_eq!(sw_alignment.similarity(), 9); // ACG vs ACG
    assert_eq!(sw_alignment.align()[0], "ACG\nACG");
}

#[test]
fn test_all_alignments_sw() {
    // Two equal local matches
    let sw = SmithWaterman::compute("ACTGNNNACTG", "ACTG");
    let all_aligned = sw.all_alignments(true).align();
    
    // It should find both "ACTG" matches if they have the same max score
    assert_eq!(all_aligned.len(), 2);
    for alignment in all_aligned {
        assert_eq!(alignment, "ACTG\nACTG");
    }
}
