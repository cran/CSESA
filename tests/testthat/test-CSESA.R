context("Test CSESA functions")

test_that( "Test GetNewSpacerCode function", {
    spacer <- GetNewSpacer("GCGCCGGGAACACCAACGTCGGTTTATCCCCGCTGGCGCCGGGAACACAGGCGGACCGAAAAACCGTTTTCAGCCAACGTCGGTTTATCCCCGCTGGCGCCGGGAACACCAACGTCGGTTT")
    expect_equal(spacer, "AGGCGGACCGAAAAACCGTTTTCAGCCAACGT")

    spacerCode <- GetNewSpacerCode("GCGCCGGGAACACCAACGTCGGTTTATCCCCGCTGGCGCCGGGAACACAGGCGGACCGAAAAACCGTTTTCAGCCAACGTCGGTTTATCCCCGCTGGCGCCGGGAACACCAACGTCGGTTT")
    expect_equal(spacerCode, "Ent8")
})

test_that( "Test CSESA function with PCR method", {
    outString <- "The newly incorporated spacer in the first CRISPR sequence: Ent8\\nThe newly incorporated spacer in the second CRISPR sequence: EntB9\\nPredicted serotype\\(s\\): \\[Enteritidis 95%\\] \\[Nitra 2%\\] \\[Rosenberg 2%\\] \\[Blegdam 1%\\]"
    expect_output(CSESA(method = "PCR"), "ERROR : No such file\\(s\\)! ")
    expect_output(CSESA(in.file1 = "1.txt", method = "PCR"), "ERROR : 1\\.txt is not existed! ")
    expect_output(CSESA(system.file("extdata", "sequence_CRIPSR1.fasta", package = "CSESA"), system.file("extdata", "sequence_CRIPSR2.fasta", package = "CSESA"), method = "PCR"), outString)
})
