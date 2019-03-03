#' CSESA (CRISPR-based Salmonella enterica Serotype Analyzer).
#'
#' @description The main function in CSESA package.
#'
#' @param in.file1 The first input file, the default value is NULL.
#' @param in.file2 The second input file (optional), the default value is NULL.
#' @param out.file Into which results will be saved if this value is set. Otherwise results will be displayed on the screen.
#' @param method  The method to handle the input file(s), which can be set as "PCR" or "WGS". Choose "PCR" if the CRISPR sequence(s) from PCR amplification is entered, and choose "WGS" when entering the whole genome assembly of a Salmonella isolate.
#'
#' @note If you use the "WGS" method, please make sure you have installed the BLAST software and included it within the working path.
#'
#' @examples
#'   CSESA(system.file("extdata", "sequence_CRIPSR1.fasta", package = "CSESA"),
#'   system.file("extdata", "sequence_CRIPSR2.fasta", package = "CSESA"), method = "PCR")
#'   CSESA(system.file("extdata", "sequence_CRIPSR1.fasta", package = "CSESA"), method = "PCR")
#'   CSESA(system.file("extdata", "Salmonella_whole_genome_assembly.fasta",
#'   package = "CSESA"), method = "WGS")
#' @importFrom utils read.table
#' @import Biostrings
#'
#' @export
#'
CSESA <- function(in.file1 = NULL, in.file2 = NULL, out.file = NULL, method = c("PCR", "WGS")) {
    tryCatch({
        if (is.null(in.file1) && is.null(in.file2)) {
            stop("No such file(s)!")
        }
        method <- toupper(method)
        method <- match.arg(method)
        if (method == "PCR") {
            seq1 <- ReadInFile(in.file1)
            seq2 <- ReadInFile(in.file2)
            PCR(seq1, seq2, out.file)
        }
        else {
            if (is.null(in.file1) == FALSE && file.exists(in.file1)) {
                file <- in.file1
                if (is.null(in.file2) == FALSE)
                    print("Warning: under the WGS mode, CSESA would ignore the second file when receiving two input files.")
            }
            else if (is.null(in.file2) == FALSE && file.exists(in.file2)){
                file <- in.file2
            }
            else {
                stop("The input file(s) does not exist!")
            }
            wgs.list <- WGS(file)
            seq1 <- wgs.list$seq1
            seq2 <- wgs.list$seq2
            PCR(seq1, seq2, out.file)
        }
    }, error = function(e) {
        cat("ERROR :",conditionMessage(e),"\n")
    })
}


#' Get the CSESA obeject through the two sequence.
#'
#' @param seq1 The first DNA sequence.
#' @param seq2 The second DNA sequence.
#' @param out.file Into which results will be saved if this value is set. Otherwise results will be displayed on the screen.
#'
PCR <- function(seq1, seq2, out.file) {
    csesa.result <- list()

    csesa.result$spacer1 <- GetAllNewSpacers(seq1)
    csesa.result$spacer2 <- GetAllNewSpacers(seq2)

    csesa.result$serotype <- FindSerotype(csesa.result$spacer1, csesa.result$spacer2)
    class(csesa.result) <- "CSESA"

    if (is.null(out.file)) {
        cat(GetStr(csesa.result))
    }
    else {
        write(GetStr(csesa.result), file = out.file)
    }
}


#' Find the serotype based on the analysis of the new spacers.
#'
#' @param file The input fasta file.
#'
#' @return The two DNA molecular sequence.
#'
WGS <- function(file) {
    path <- Sys.which("blastn")
    if (all(path == "")) {
        stop("Blast does not exist!")
    }

    # loading the database
    db <- system.file("primerDB", package = "CSESA")
    db <- file.path(db, "primers")
    db2 <- system.file("primerB2DB", package = "CSESA")
    db2 <- file.path(db2, "primerB2.fasta")

    blastn <- Sys.which("blastn")

    if (all(blastn == "")) {
        stop("The BLAST software has not been employed. Please install it first and check it within the working path.")
    }
    blastn = blastn[which(blastn != "")[1]]

    tmpwd <- tempdir()
    curwd <- getwd()
    tmp.prefix <- basename(tempfile(tmpdir = tmpwd))
    on.exit({
        file.remove(Sys.glob(paste(tmp.prefix, "*")))
        setwd(curwd)
    })
    setwd(tmpwd)

    outfile <- paste(tmp.prefix, ".out", sep = "")
    infile <- paste(tmp.prefix, ".in", sep = "")

    text <- readLines(file)
    text <- gsub(pattern = " ", replacement = "_", x = text)
    text <- gsub(pattern = ",|#", replacement = "", x = text)
    writeLines(text, con = infile)

    data <- readDNAStringSet(infile)

    system(paste(blastn, "-db", db, "-query", infile, "-out", outfile, "-outfmt 10", "-task blastn-short"), ignore.stdout = FALSE, ignore.stderr = FALSE)

    result.table <- read.table(outfile, sep=",", quote = "")
    colnames(result.table) <- c("Query_id",  "Subject_id", "Perc_ident",
                                "Align_length", "Mismatches", "Gap_openings", "Query_start", "Query_end",
                                "S_start", "S_end", "E", "Bits" )

    config.primerA2 = result.table[which(result.table$Subject_id == "locusA_primer_R" & result.table$Align_length == 23), ]
    config.primerA1 = result.table[which(result.table$Subject_id == "locusA_primer_F" & result.table$Align_length == 20), ]
    config.primerB1 = result.table[which(result.table$Subject_id == "locusB_primer_F" & result.table$Align_length == 25), ]
    seq1 <- NA
    seq2 <- NA
    # primer A
    if (nrow(config.primerA1) == 1 && nrow(config.primerA2) == 1 && config.primerA1$Query_id == config.primerA2$Query_id) {
        id <- as.character(config.primerA1$Query_id)
        locus.primerA1 = config.primerA1$Query_start
        locus.primerA2 = config.primerA2$Query_start
        seq <- ""
        idx <- 1
        for (tn in names(data)) {
            if (startsWith(tn, id)) {
                seq <- as.character(unlist(data[idx, ]))
                break
            }
            idx <- idx + 1
        }
        seq1 <- substr(seq, min(locus.primerA1, locus.primerA2), max(locus.primerA1, locus.primerA2))
    }
    # primer B
    if (nrow(config.primerB1) == 1) {
        id = as.character(config.primerB1$Query_id)
        locus.primerB1 = config.primerB1$Query_start

        # find the primerB2
        system(paste(blastn, "-db", db2, "-query", infile, "-out", outfile, "-outfmt 10", "-task blastn-short"), ignore.stdout = FALSE, ignore.stderr = FALSE)
        result.table <- read.table(outfile, sep=",", quote = "")
        colnames(result.table) <- c("Query_id",  "Subject_id", "Perc_ident",
                                    "Align_length", "Mismatches", "Gap_openings", "Query_start", "Query_end",
                                    "S_start", "S_end", "E", "Bits")

        seq <- ""
        idx <- 1
        for (tn in names(data)) {
            if (startsWith(tn, id)) {
                seq <- as.character(unlist(data[idx, ]))
                break
            }
            idx <- idx + 1
        }
        if (nrow(result.table) >= 1) {
            config.primerB2 = result.table[1, ]
            locus.primerB2 = config.primerB2$Query_start
            # If the position of primerB2 is in the range [locus.primerB1 - 3000, locus.primerB1 + 3000]
            if (locus.primerB2 >= locus.primerB1 - 3000 && locus.primerB2 <= locus.primerB1 + 3000) {
                seq2 <- substr(seq, min(locus.primerB1, locus.primerB2), max(locus.primerB1, locus.primerB2))
            }
            else {
                seq2 <- substr(seq, locus.primerB1 - 3000, locus.primerB1 + 3000)
            }
        }
        else {
            seq2 <- substr(seq, locus.primerB1 - 3000, locus.primerB1 + 3000)
        }
    }
    return (list(seq1 = seq1, seq2 = seq2))
}


#' Get the new spacers from the molecular sequence and its reverse complement.
#'
#' @param molecular.seq The molecular sequence.
#' @return The vector of the new spacers, which is extracted from the molecular sequence and its reverse complement.
#'
#' @note If there doesn't exist any new spacer, the function would return NA.
#'
GetAllNewSpacers <- function(molecular.seq = NULL) {
    if (is.null(molecular.seq) || is.na(molecular.seq) || molecular.seq == "")
        return (NA)
    molecular.seq.rev <- GetReverseComplement(molecular.seq)

    # handle the cases specific to Typhi
    typhi <- "ACGGCTATCCTTGTTGACGTGGGGAATACTGCTACACGCAAAAATTCCAGTCGTTGGCGCACGGTTTATCCCCGCTGGCGCGGGGAACAC"
    if(grepl(typhi,molecular.seq) || grepl(typhi,molecular.seq.rev))
        return (c("EntB0var1"))

    new.spacer <- GetNewSpacerCode(molecular.seq)
    new.spacer.rev <- GetNewSpacerCode(molecular.seq.rev)
    new.spacer.arr <- character()
    if (is.null(new.spacer) == FALSE && is.na(new.spacer) == FALSE)
        new.spacer.arr <- c(new.spacer.arr, new.spacer)
    if (is.null(new.spacer.rev) == FALSE && is.na(new.spacer.rev) == FALSE)
        new.spacer.arr <- c(new.spacer.arr, new.spacer.rev)
    if (length(new.spacer.arr) == 0) {
        return (NA)
    }
    else return (new.spacer.arr)
}


#' Find the serotype based on the analysis of the new spacers.
#'
#' @param csesa1 The new spacer of the first sequence.
#' @param csesa2 The new spacer of the second sequence.
#'
#' @return The data frame which represents the serotype.
#'
FindSerotype <- function(csesa1 = NA, csesa2 = NA) {
    if (is.na(csesa1) == TRUE && is.na(csesa2) == TRUE) {
        stop("Sorry. We did not find any corresponding serotype in the lib!")
    }
    V1 = V2 = V3 = V4 = NULL
    mapping.table <- read.table(system.file("packageData", "mapping_tbl.txt", package = "CSESA"), sep = "\t")
    if (is.na(csesa1) == TRUE || is.na(csesa2) == TRUE) {
        csesa <- csesa1
        if (is.na(csesa1))
            csesa <- csesa2
        serotype <- subset(mapping.table, is.element(V1, csesa) | is.element(V2, csesa), select = V3)
        if (nrow(serotype) > 1)
            serotype <- unique(serotype)
    }
    else {
        serotype <- subset(mapping.table, is.element(V1, csesa1) & is.element(V2, csesa2) |
                               is.element(V1, csesa2) & is.element(V2, csesa1), select = c(V3, V4))
    }
    return (serotype)
}


#' Get the new spacer from the molecular sequence and map it to the code.
#'
#' @param molecular.seq The molecular sequence.
#' @return The new spacer code as a string.
#'
GetNewSpacerCode <- function(molecular.seq = NULL) {
    if (is.null(molecular.seq))
        return (NULL)
    new.spacer <- GetNewSpacer(molecular.seq)
    if (is.null(new.spacer))
        return (NULL)
    V1 = V3 = NULL
    spacers.table <- read.table(system.file("packageData", "spacer_tbl.txt", package = "CSESA"))
    spacer.char <- as.character(subset(spacers.table, V3 == new.spacer, select = V1)[1, 1])
}


#' Get the new spacer from the molecular sequence.
#'
#' @param molecular.seq The molecular sequence.
#' @return The new spacer sequence as a string.
#'
#' @examples
#' GetNewSpacer("AGAGGCGGACCGAAAAACCGTTTTCAGCCAACGTAT")
#'
#' @export
GetNewSpacer <- function(molecular.seq = NULL) {
    if (is.null(molecular.seq))
        return (NULL)
    max.match <- "-"
    max.count <- -1
    dr.table <- read.table(system.file("packageData", "dr_tbl.txt", package = "CSESA"))
    for (x in dr.table$V3[-1]) {
        t <- gregexpr(pattern = x, text = molecular.seq)
        if (t[[1]][1] != -1 && length(t[[1]]) > max.count) {
            max.match = x
            max.count = length(t[[1]])
        }
    }
    if (max.count < 1)
        return (NULL)
    spacers <- strsplit(molecular.seq, max.match)[[1]]
    if (substr(molecular.seq, nchar(molecular.seq) - nchar(max.match) + 1, nchar(molecular.seq)) == max.match) {
        return (spacers[length(spacers)])
    }
    spacers[length(spacers) - 1]
}


#' Get the information string from the CSESA s3 object.
#'
#' @param csesa The S3 object CSESA.
#' @return The string record the newly spacers and serotype information.
#'
GetStr <- function(csesa) {
    if (is.null(csesa)) {
        stop("The csesa object should be set!")
    }

    if (is.na(csesa$spacer1))
        str <- "The newly incorporated spacer of the first CRISPR sequence is not available for prediction.\n"
    else
        str <- paste("The newly incorporated spacer in the first CRISPR sequence: ", csesa$spacer1, "\n", sep = "")
    if (is.na(csesa$spacer2))
        str <- paste0(str, "The newly incorporated spacer of the second CRISPR sequence is not available for prediction.\n")
    else
        str <- paste(str, "The newly incorporated spacer in the second CRISPR sequence: ", csesa$spacer2, "\n", sep = "")

    if (is.na(csesa$spacer1) || is.na(csesa$spacer2)) {
        result <- ""
        if (is.atomic(csesa$serotype))
            result <- csesa$serotype
        else
            result <- paste(csesa$serotype[, 1], collapse = "] [")
        result <- paste("Predicted serotype(s): [", result, sep = "")
        str <- paste(str, result, "]", "\n", sep = "")
    }
    else {
        result <- paste(csesa$serotype[, 1], csesa$serotype[, 2])
        result <- paste("Predicted serotype(s): [", paste(result, collapse = "] ["), sep = "")
        str <- paste(str, result, "]", "\n", sep = "")
    }
    str
}


#' Read the three types of input file.
#'
#' @param file.name The input file name.
#' @return The molecular sequence as a string.
#'
ReadInFile <- function(file.name) {
    if (is.null(file.name) || is.na(file.name))
        return (NA)
    if (file.exists(file.name) == FALSE) {
        stop(paste0(file.name, " does not existed!"))
    }
    data <- scan(file.name, what = "", quiet = TRUE)
    if (substring(data[1], 1, 1) == '>')
        data <- data[-1]
    toupper(paste(data, collapse = ""))
}


#' Return the reverse complement of the sequence.
#'
#' @param x The input sequence.
#' @return The reverse complement sequence as a string.
#'
GetReverseComplement <- function(x) {
    a <- chartr("ATGC","TACG",x)
    paste(rev(substring(a, 1 : nchar(a), 1 : nchar(a))), collapse = "")
}
