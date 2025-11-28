use goombay_rs::align::{AlignmentData, NeedlemanWunsch, Scoring};

fn main() {
    // Set scoring parameters for Needleman Wunsch
    let scores = Scoring {
        identity: 2,
        mismatch: 1,
        gap: 1,
        transpose: None,    // This is an Optional parameter
        extended_gap: None, // This is an Optional parameter
    };
    let nw = NeedlemanWunsch::new(&scores);

    // Set sequences to be aligned
    let query = "attain";
    let subject = "atin";
    // Create struct for above sequences so matrices don't need to be
    // recalculated when scoring alignments
    let mut nw_alignment = AlignmentData::new(query, subject);

    // Build the scoring matrix and pointer matrix
    nw.build(&mut nw_alignment);

    // Align the sequences based on the pointer matrix
    let aligned = nw.align(&nw_alignment, true);
    println!("{}", nw_alignment.score_matrix());
    println!("{}", nw_alignment.pointer_matrix());
    for (i, alignment) in aligned.iter().enumerate() {
        println!("{}.", i + 1);
        println!("{alignment}");
    }

    // Calculate alignment scores for aligned sequences
    let sim = nw.similarity(&nw_alignment);
    let dist = nw.distance(&nw_alignment);
    let norm_sim = nw.normalized_similarity(&nw_alignment);
    let norm_dist = nw.normalized_distance(&nw_alignment);
    println!(
        "Similarity: {sim}\nDistance: {dist}\nNormalized Similarity: {norm_sim}\nNormalized Distance: {norm_dist}"
    );
}
