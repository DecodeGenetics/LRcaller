#include <cstddef>
// BEFORE EVERYTHING
inline size_t lrcaller_bgzf_threads = 1;
#define SEQAN_BGZF_NUM_THREADS lrcaller_bgzf_threads

#include "algo.hpp"
#include "misc.hpp"
#include "options.hpp"

void mainProgram(LRCOptions & O)
{
    // THIS NEEDS TO BE SET BEFORE ANY BAM OBJECTS ARE DECLARED; parallelism happens on higher level, so this is 1
    lrcaller_bgzf_threads = 1;

    if (O.verbose)
        std::cerr << "Number of threads requested: " << O.nThreads << ". Got: " << omp_get_max_threads() << ".\n";

    seqan::VcfFileIn              vcfIn(O.vcfInFile.c_str());
    seqan::VcfHeader              header;
    std::vector<seqan::VcfRecord> vcfRecords;

    readHeader(header, vcfIn);
    while (!atEnd(vcfIn))
    {
        seqan::VcfRecord r;
        readRecord(r, vcfIn);

        vcfRecords.push_back(std::move(r));
    }

    // Open the input VCF file and prepare output VCF stream.
    appendValue(header, seqan::VcfHeaderRecord("FORMAT", "<ID=GT,Number=1,Type=String,Description=\"Genotype\">"));
    appendValue(header,
                seqan::VcfHeaderRecord("FORMAT",
                                       "<ID=AD,Number=3,Type=Integer,Description=\"Allelic depths from alignment "
                                       "supporting ref and alt allele and total number of reads\">"));
    appendValue(header,
                seqan::VcfHeaderRecord("FORMAT",
                                       "<ID=VA,Number=3,Type=Integer,Description=\"Allelic depths from bam file "
                                       "supporting ref and alt allele and total number of reads\">"));
    appendValue(
      header,
      seqan::VcfHeaderRecord("FORMAT",
                             "<ID=PL,Number=G,Type=Integer,Description=\"PHRED-scaled genotype likelihoods\">"));
    std::ofstream     vcfStream(O.vcfOutFile.c_str());
    seqan::VcfFileOut vcfOut(vcfIn);
    open(vcfOut, vcfStream, seqan::Vcf());
    writeHeader(vcfOut, header);

//     bool useBam2 = false;
//     if (O.bam2 != "")
//     {
//         useBam2 = true;
//         if (O.verbose)
//             std::cerr << "Using bam2 " << O.bam2 << '\n';
//     }

    if (O.cacheDataInTmp)
    {
        O.cacheDir         = std::filesystem::temp_directory_path() / "tmp-lrcaller-XXXXXX";
        std::string buffer = O.cacheDir.string();
        char * r = mkdtemp(buffer.data());
        if (r == nullptr)
        {
            std::cerr << "ERROR: Could not create temporary directory.\n";
            std::exit(1);
        }
        O.cacheDir = buffer;
        if (O.verbose)
            std::cerr << "Using local cachedir: " << O.cacheDir << '\n';
    }

    // Build an index of the insertion sequences' fasta file, if such a file is used.
    if (O.verbose)
        std::cerr << "Genotyping VCF records in \'" << O.vcfInFile << std::endl;

    omp_set_num_threads(O.nThreads);

    struct thread_cache_t
    {
        std::vector<seqan::BamAlignmentRecord> bars;
        seqan::CharString                      chrom;

        std::vector<seqan::BamFileIn>            bamFiles;
        std::vector<seqan::BamIndex<seqan::Bai>> bamIndexes;

        seqan::FaiIndex faIndex;
    };

    std::vector<thread_cache_t> per_thread;
    per_thread.resize(O.nThreads);

    for (thread_cache_t & c : per_thread)
    {
        parseBamFileName(O.bam, c.bamFiles, c.bamIndexes, O);

        if (!open(c.faIndex, O.faFile.c_str()))
            if (!build(c.faIndex, O.faFile.c_str()))
                throw error{"Could neither find nor build the index of ", O.faFile};
    }

    //     if (useBam2)
    //         parseBamFileName(seqan::toCString(O.bam2), bamIndex2Handles, bam2Handles, O);

    // split input into chunks of adjacent variants so that reads are only read once
    std::vector<std::span<seqan::VcfRecord>> chunks;
    size_t                                   chunk_first = 0;
    for (size_t i = 1; i < length(vcfRecords); ++i)
    {
        seqan::VcfRecord & var = vcfRecords[i];
        if (var.rID == -1)
            throw error{"Invalid ID in VCF record number: ", i};

        if (seqan::VcfRecord & lastVar = vcfRecords[i - 1];
            (var.rID != lastVar.rID) || (var.beginPos > lastVar.beginPos + (ssize_t)O.wSize)) // new chunk
        {
            chunks.emplace_back(vcfRecords.begin() + chunk_first, vcfRecords.begin() + i);
            chunk_first = i;
        }
    }
    // last chunk
    chunks.emplace_back(vcfRecords.begin() + chunk_first, vcfRecords.begin() + vcfRecords.size());

    #pragma omp parallel for
    for (size_t i = 0; i < chunks.size(); ++i)
    {
        std::span<seqan::VcfRecord> chunk = chunks[i];
        thread_cache_t & thread_cache = per_thread[omp_get_thread_num()];

        thread_cache.bars.clear();
        thread_cache.chrom = seqan::contigNames(seqan::context(vcfIn))[chunk.begin()->rID];

        processChunk(thread_cache.bamFiles,
                     thread_cache.bamIndexes,
                     thread_cache.faIndex,
                     thread_cache.chrom,
                     thread_cache.bars,
                     chunk,
                     O);
    }

    for (seqan::VcfRecord /*const ? */ & var : vcfRecords)
    {
        writeRecord(vcfOut, var);
    }

    if (O.cacheDataInTmp)
    {
        if (O.verbose)
            std::cerr << "Removing cache directory...";

        std::filesystem::remove_all(O.cacheDir);

        if (O.verbose)
            std::cerr << "done.\n";
    }
}

int main(int argc, char const ** argv)
{
    LRCOptions O;

    try
    {
        auto res = parseLRCArguments(argc, argv, O);

        if (res == seqan::ArgumentParser::PARSE_ERROR)
            throw error{"Could not parse command line arguments."};
        else if (res == seqan::ArgumentParser::PARSE_OK)
            mainProgram(O);
        // else the help page was shown
    }
    catch (error const & e)
    {
        std::cerr << "ERROR: " << e.what() << '\n';
        return 1;
    }

    return 0;
}
