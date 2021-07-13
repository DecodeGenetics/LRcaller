#pragma once

#include <filesystem>
#include <string>

#include <seqan/arg_parse.h>

enum class genotyping_model
{
    multi,
    ad,
    va,
    presence,
    va_old,
    joint
};

/*  Option handling */
struct LRCOptions
{
    std::string       bam;
    seqan::CharString bam2;
    std::string       faFile = "genome.fa";
    std::string       vcfInFile;
    std::string       vcfOutFile;
    seqan::CharString genotypingModel = "multi";

    bool   dynamicWSize            = false;
    size_t wSize                   = 500; // Window Size,
    size_t maxBARcount             = 200; // maximum number of reads to use for a variant
    bool   verbose                 = false;
    bool   genotypeRightBreakpoint = false;
    double logScaleFactor          = 2.0;
    size_t bandedAlignmentPercent  = 100;

    int match     = 1;
    int mismatch  = -1;
    int gapOpen   = -1;
    int gapExtend = -1;

    bool   mask                  = false;
    bool   cropRead              = true;
    bool   outputRefAlt          = false;  // For debug purposes
    size_t maxSoftClipped        = 500;    // Read cannot be softclipped at invariant end
    double maxAlignBits          = 10.0;   // needs to be added as option
    size_t overlapBits           = 10;     // needs to be added as option
    size_t varWindow             = 100;    // How far
    size_t maxReadLenForCropping = 200000; // Should not be needed
    size_t minMapQ               = 30;     // minimum std::mapQ for a read to be used
    size_t minDelIns             = 6;      // minum number of consecutive bp to be deleted or inserted to count it
    size_t minPresent            = 30;     // minum number of del/ins for a variant to be considered present

    // Fraction of deletion/insertion length that can be deleted/inserted but still called ref
    double refThreshFraction    = 1.0 / 10.0;
    // Fraction of deletion/insertion length that needs to be deleted/inserted to call alt
    double altThreshFraction    = 1.0 / 2.0;
    // Fraction of deletion/insertion length can be at most be observed to be deleted/inserted to call alt
    double altThreshFractionMax = 100.0;

    genotyping_model gtModel  = genotyping_model::multi;
    size_t           nThreads = std::thread::hardware_concurrency();

    bool                  cacheDataInTmp = false; // whether to copy BAM and BAI to /tmp at start
    std::filesystem::path cacheDir;               // where to store cache
};

inline int parseLRCArguments(int argc, char const ** argv, LRCOptions & O)
{
    seqan::ArgumentParser parser("LRcaller");
    addUsageLine(parser, "[\\fIOPTIONS\\fP]  \"\\fIBAMFILE\\fP\"  \"\\fIVCF_FILE_IN\\fP\" \"\\fIVCF_FILE_OUT\\fP\" ");
    addDescription(parser, "Genotypes variants using long reads\n in bam file/file of bam files ");

    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "BAMFILE bam file/file of bam files"));
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "VCF_FILE_IN - input vcf file"));
    addArgument(parser,
                seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "VCF_FILE_OUT - output genotyped vcf file"));
    addOption(parser, seqan::ArgParseOption("fa", "fa", "Fastafile", seqan::ArgParseArgument::STRING, "FA"));
    addOption(
      parser,
      seqan::ArgParseOption(
        "b2",
        "bam2",
        "Second bam file/file of bam files, use if you want to use multiple basecallings/mappings for the same reads",
        seqan::ArgParseArgument::STRING,
        "b2"));

    addOption(parser,
              seqan::ArgParseOption("gtm",
                                    "genotyper",
                                    "Genotyper used joint/ad/va/multi (default joint).",
                                    seqan::ArgParseArgument::STRING,
                                    "GENOTYPER"));

    addOption(parser,
              seqan::ArgParseOption("faA",
                                    "faA",
                                    "Fasta file alt, used when VCF file refers to an alternate fasta file",
                                    seqan::ArgParseArgument::STRING,
                                    "FAA"));
    addOption(parser, seqan::ArgParseOption("v", "verbose", "Verbose output"));
    addOption(parser,
              seqan::ArgParseOption("rb",
                                    "right_breakpoint",
                                    "Genotype right breakpoint, left breakpoint genotyped by default"));
    addOption(parser, seqan::ArgParseOption("ora", "get_ref_alt", "Output ref and alt allele, for debugging purposes"));

    addOption(
      parser,
      seqan::ArgParseOption("A", "match", "Match score for alignment", seqan::ArgParseArgument::INTEGER, "INT"));
    addOption(
      parser,
      seqan::ArgParseOption("B", "mismatch", "Error penalty for alignment.", seqan::ArgParseArgument::INTEGER, "INT"));
    addOption(
      parser,
      seqan::ArgParseOption("O", "gap_open", "Gap open for alignment", seqan::ArgParseArgument::INTEGER, "INT"));
    addOption(
      parser,
      seqan::ArgParseOption("E", "gap_extend", "Gap extend for alignment", seqan::ArgParseArgument::INTEGER, "INT"));
    addOption(parser,
              seqan::ArgParseOption("w", "window_size", "Window size", seqan::ArgParseArgument::INTEGER, "INT"));

    addOption(parser, seqan::ArgParseOption("", "dyn-w-size", "Dynamically adjust window size to allele length."));

    addOption(parser,
              seqan::ArgParseOption("", "cacheDataInTmp", "Copy reads and index to (local) tmp directory before run."));

    addOption(
      parser,
      seqan::ArgParseOption("", "mask", "Reduce stretches of the same base to a single base before alignment."));

    seqan::hideOption(parser, "mask"); // THIS IS CURRENTLY NOT EXPOSED TO THE USER

    setDefaultValue(parser, "A", O.match);
    setDefaultValue(parser, "B", O.mismatch);
    setDefaultValue(parser, "O", O.gapOpen);
    setDefaultValue(parser, "E", O.gapExtend);
    setDefaultValue(parser, "w", O.wSize);

    addOption(
      parser,
      seqan::ArgParseOption("ab", "align_bits", "Max bits for alignment", seqan::ArgParseArgument::DOUBLE, "DOUBLE"));
    setDefaultValue(parser, "ab", O.maxAlignBits);

    addOption(parser,
              seqan::ArgParseOption("ob",
                                    "overlap_bits",
                                    "Score for overlap (log_2)",
                                    seqan::ArgParseArgument::INTEGER,
                                    "INT"));
    setDefaultValue(parser, "ob", O.overlapBits);

    addOption(
      parser,
      seqan::ArgParseOption("nt", "number_of_threads", "Number of threads", seqan::ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "nt", O.nThreads);

    addOption(parser,
              seqan::ArgParseOption("vw",
                                    "var_window",
                                    "Var window size, look for del/ins inside this window",
                                    seqan::ArgParseArgument::INTEGER,
                                    "INT"));
    setDefaultValue(parser, "vw", O.varWindow);

    addOption(
      parser,
      seqan::ArgParseOption("mapQ", "read_mapq", "Minimum read std::mapQ", seqan::ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "mapQ", O.minMapQ);

    addOption(parser,
              seqan::ArgParseOption("di",
                                    "min_del_ins",
                                    "Minimum del/ins in cigar to consider",
                                    seqan::ArgParseArgument::INTEGER,
                                    "INT"));
    setDefaultValue(parser, "di", O.minDelIns);

    addOption(parser,
              seqan::ArgParseOption("mp",
                                    "min_present",
                                    "Minimum del/ins var a variant to be considered present",
                                    seqan::ArgParseArgument::INTEGER,
                                    "INT"));
    setDefaultValue(parser, "mp", O.minPresent);

    addOption(parser,
              seqan::ArgParseOption("",
                                    "max-alignments",
                                    "Maximum alignments per variant to consider",
                                    seqan::ArgParseArgument::INTEGER,
                                    "INT"));
    setDefaultValue(parser, "max-alignments", O.maxBARcount);

    addOption(parser,
              seqan::ArgParseOption("rt",
                                    "ref_thresh_fraction",
                                    "Threshold for fraction of del/ins bp to be considered ref",
                                    seqan::ArgParseArgument::DOUBLE,
                                    "DOUBLE"));
    setDefaultValue(parser, "rt", O.refThreshFraction);

    addOption(parser,
              seqan::ArgParseOption("at",
                                    "alt_thresh_fraction",
                                    "Threshold for fraction of del/ins bp to be considered alt",
                                    seqan::ArgParseArgument::DOUBLE,
                                    "DOUBLE"));
    setDefaultValue(parser, "at", O.altThreshFraction);

    addOption(parser,
              seqan::ArgParseOption("am",
                                    "alt_thresh_fraction_max",
                                    "Max threshold for fraction of del/ins bp to be considered alt",
                                    seqan::ArgParseArgument::DOUBLE,
                                    "DOUBLE"));
    setDefaultValue(parser, "am", O.altThreshFractionMax);

    addOption(parser,
              seqan::ArgParseOption("lsf",
                                    "logScaleFactor",
                                    "Log scale factor for comparing alignment scores",
                                    seqan::ArgParseArgument::DOUBLE,
                                    "DOUBLE"));

    addOption(parser,
              seqan::ArgParseOption("",
                                    "band",
                                    "Percentage of window size to use as band in alignment (100 == no band).",
                                    seqan::ArgParseArgument::INTEGER,
                                    "INTEGER"));
    setDefaultValue(parser, "band", O.bandedAlignmentPercent);

    setDefaultValue(parser, "fa", O.faFile);
    setDefaultValue(parser, "b2", O.bam2);
    setValidValues(parser, "genotyper", "joint ad va multi");
    setDefaultValue(parser, "genotyper", "joint");
    setDefaultValue(parser, "logScaleFactor", O.logScaleFactor);

    seqan::ArgumentParser::ParseResult res = parse(parser, argc, argv);
    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;
    getArgumentValue(O.bam, parser, 0);
    getArgumentValue(O.vcfInFile, parser, 1);
    getArgumentValue(O.vcfOutFile, parser, 2);

    if (isSet(parser, "align_bits"))
        getOptionValue(O.maxAlignBits, parser, "align_bits");

    if (isSet(parser, "overlap_bits"))
        getOptionValue(O.overlapBits, parser, "overlap_bits");

    if (isSet(parser, "fa"))
        getOptionValue(O.faFile, parser, "fa");

    if (isSet(parser, "number_of_threads"))
        getOptionValue(O.nThreads, parser, "number_of_threads");

    if (isSet(parser, "var_window"))
        getOptionValue(O.varWindow, parser, "var_window");

    if (isSet(parser, "read_mapq"))
        getOptionValue(O.minMapQ, parser, "read_mapq");

    if (isSet(parser, "min_del_ins"))
        getOptionValue(O.minDelIns, parser, "min_del_ins");

    if (isSet(parser, "min_present"))
        getOptionValue(O.minPresent, parser, "min_present");

    if (isSet(parser, "max-alignments"))
        getOptionValue(O.maxBARcount, parser, "max-alignments");

    if (isSet(parser, "ref_thresh_fraction"))
        getOptionValue(O.refThreshFraction, parser, "ref_thresh_fraction");

    if (isSet(parser, "alt_thresh_fraction"))
        getOptionValue(O.altThreshFraction, parser, "alt_thresh_fraction");
    if (isSet(parser, "alt_thresh_fraction_max"))
        getOptionValue(O.altThreshFractionMax, parser, "alt_thresh_fraction_max");

    if (isSet(parser, "match"))
        getOptionValue(O.match, parser, "match");
    if (isSet(parser, "mismatch"))
        getOptionValue(O.mismatch, parser, "mismatch");
    if (isSet(parser, "gap_open"))
        getOptionValue(O.gapOpen, parser, "gap_open");
    if (isSet(parser, "gap_extend"))
        getOptionValue(O.gapExtend, parser, "gap_extend");
    if (isSet(parser, "window_size"))
        getOptionValue(O.wSize, parser, "window_size");
    if (isSet(parser, "bam2"))
        getOptionValue(O.bam2, parser, "bam2");
    if (isSet(parser, "band"))
        getOptionValue(O.bandedAlignmentPercent, parser, "band");

    seqan::CharString gtModelName;
    if (isSet(parser, "genotyper"))
        getOptionValue(gtModelName, parser, "genotyper");
    if (gtModelName == "multi")
        O.gtModel = genotyping_model::multi;
    else if (gtModelName == "ad")
        O.gtModel = genotyping_model::ad;
    else if (gtModelName == "va")
        O.gtModel = genotyping_model::va;
    else if (gtModelName == "va2")
        O.gtModel = genotyping_model::va_old;
    else if (gtModelName == "p")
        O.gtModel = genotyping_model::presence;
    else if (gtModelName == "joint")
        O.gtModel = genotyping_model::joint;
    if (isSet(parser, "logScaleFactor"))
        getOptionValue(O.logScaleFactor, parser, "logScaleFactor");
    // O.cropRead = isSet(parser, "cropread");
    O.verbose                 = isSet(parser, "verbose");
    O.genotypeRightBreakpoint = isSet(parser, "right_breakpoint");
    O.outputRefAlt            = isSet(parser, "get_ref_alt");
    O.cacheDataInTmp          = isSet(parser, "cacheDataInTmp");
    O.dynamicWSize            = isSet(parser, "dyn-w-size");
    O.mask                    = isSet(parser, "mask");
    // get options
    return res;
}
