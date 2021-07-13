#pragma once

#include <cassert>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <map>
#include <omp.h>
#include <set>
#include <span>
#include <string>
#include <time.h>
#include <vector>

#include <seqan/align.h>
#include <seqan/arg_parse.h>
#include <seqan/bam_io.h>
#include <seqan/basic.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/store.h>
#include <seqan/stream.h>
#include <seqan/vcf_io.h>

#include "misc.hpp"
#include "options.hpp"

// Sequence, alignment, and alignment row.
typedef seqan::String<seqan::Dna5>                TSequence;
typedef seqan::Align<TSequence, seqan::ArrayGaps> TAlign;
typedef seqan::Row<TAlign>::Type                  TRow;
using ssize_t = std::make_signed_t<size_t>;

#define LL_THRESHOLD -25.5
#define LG10         3.322
#define NO_ALIGNMENT -10000

inline constexpr size_t NO_BEST = size_t(-1ull);

/* Stores information of how a read aligns across a variant */
class varAlignInfo
{
public:
    size_t              nD;
    size_t              nI;
    size_t              nAlleles;
    std::vector<double> alignS;
    bool                softClipped;
    bool                alignsLeft;
    bool                alignsRight;
    // int refAlignLen;

    // Supports and rejects logic is not correct for very long variants
    // Alternate allele is supported as judged from bam alignment record
    bool supports(float const refLen, float const altLen, LRCOptions const & O) const
    {
        if (altLen > refLen)
        { // insertion (these are simplistic for insertion/deletion type variants)
            // doesn't work properly if alt and ref are of similar size

            return ((alignsLeft && alignsRight && ((float)nI > altLen * O.altThreshFraction) &&
                     ((float)nI < altLen * O.altThreshFractionMax)) ||
                    softClipped);
        }
        else
        {
            return ((alignsLeft && alignsRight && ((float)nD > refLen * O.altThreshFraction) &&
                     ((float)nD < refLen * O.altThreshFractionMax)) ||
                    softClipped);
        }
    }

    // Alternate allele is rejected as judged from bam alignment record
    bool rejects(float const refLen, float const altLen, LRCOptions const & O) const
    {
        if (altLen > refLen)
        { // insertion
            return ((alignsLeft && alignsRight && ((float)nI < altLen * O.refThreshFraction)) && (!softClipped));
        }
        else
        {
            return ((alignsLeft && alignsRight && ((float)nD < refLen * O.refThreshFraction)) && (!softClipped));
        }
    }

    bool present(LRCOptions const & O) const
    {
        return (nI >= O.minPresent || nD >= O.minPresent);
    }

    bool aligns() const
    {
        return (alignsLeft && alignsRight);
    }

    void reset()
    {
        nD = 0;
        nI = 0;

        for (size_t i = 0; i < alignS.size(); i++)
            alignS[i] = NO_ALIGNMENT;

        softClipped = false;
        alignsLeft  = false;
        alignsRight = false;
    }

    // Likelihood of variant relative to the most likely,
    // Updates pref values. A value x, represents that the allele is 2^-x less likely than the most likely allele
    // Returns an index to the most likely allele, if one exists
    size_t alignmentPreference(size_t const wSizeActual, LRCOptions const & O, std::vector<double> & pref) const
    {
        int    maxScore      = alignS[0];
        size_t maxI          = 0;
        int    minAlignScore = static_cast<double>(wSizeActual) * 1.2;

        for (size_t i = 0; i < nAlleles; i++)
        {
            // pref[i] = 0;
            if (alignS[i] > maxScore)
            {
                maxI     = i;
                maxScore = alignS[i];
            }
        }

        if (maxScore == NO_ALIGNMENT || maxScore <= minAlignScore)
        {
            return NO_BEST;
        }
        else
        {
            for (size_t i = 0; i < nAlleles; i++)
            {
                float d = (maxScore - alignS[i]) / O.logScaleFactor;
                if (alignS[i] == NO_ALIGNMENT || alignS[i] <= minAlignScore)
                {
                    d = (maxScore - minAlignScore) / O.logScaleFactor;
                }

                if (d > O.maxAlignBits)
                    d = O.maxAlignBits;
                if (d < 0 && O.verbose)
                    std::cerr << "WTF negative alignment score";
                pref[i] += d;
            }

            return maxI;
        }
    }

    // Likelihood of variant relative to the most likely,
    // Increments pref values. A value x, represents that the allele is 2^-x less likely than the most likely allele
    // Returns an index to the most likely allele, if one exists
    size_t vaPreference(LRCOptions const &          O,
                        size_t const                refLen,
                        std::vector<size_t> const & altLens,
                        std::vector<double> &       pref) const
    {
        if (softClipped)
        { // Softclipped, doesn't support the reference, all other alleles are equally likely
            pref[0] += O.overlapBits;
            return NO_BEST;
        }

        if (!alignsLeft || !alignsRight)
        {
            return NO_BEST;
        }

        // Count number of deletions && insertions, find the variant that is closest in size

        ssize_t insDel = nI - nD;
        ssize_t minD   = std::abs(insDel);
        size_t  minDi  = 0;

        for (size_t i = 1; i < nAlleles; i++)
        {
            ssize_t cD = altLens[i - 1] - refLen;
            if (std::abs(cD - insDel) < minD)
            {
                minDi = i;
                minD  = std::abs(cD - insDel);
            }
        }

        for (size_t i = 0; i < nAlleles; i++)
        {
            if (i != minDi)
                pref[i] += O.overlapBits;
        }

        return minDi;
    }

    //  varAlignInfo(): nD(0),nI(0),
    //  refS(NO_ALIGNMENT),altS(NO_ALIGNMENT),softClipped(false),alignsLeft(false),alignsRight(false){};
    varAlignInfo(size_t const nAllelesIn)
    {
        nAlleles = nAllelesIn;
        nD       = 0;
        nI       = 0;
        alignS.resize(nAlleles);

        for (size_t i = 0; i < nAlleles; i++)
        {
            alignS[i] = NO_ALIGNMENT;
        }

        softClipped = false;
        alignsLeft  = false;
        alignsRight = false;
    }

    varAlignInfo() : varAlignInfo{2}
    {}
};

/* Turns genotyping into std::string */
void getGtString(std::vector<double> &       lls,
                 std::vector<size_t> const & ads,
                 std::vector<size_t> const & vas,
                 std::string &               gtString)
{
    size_t gtLen = lls.size();
    for (size_t i = 0; i < gtLen; i++)
        lls[i] = -lls[i]; // Really silly hack for historical reasons, would confuse the hell out of me to fix it

    //  int max = 0;
    double maxP  = lls[0];
    size_t a1    = 0;
    size_t a2    = 0;
    size_t maxA1 = 0;
    size_t maxA2 = 0;

    for (size_t i = 0; i < gtLen; i++)
    {
        if (lls[i] > maxP)
        {
            maxP  = lls[i];
            //      max = i;
            maxA1 = a1;
            maxA2 = a2;
        }

        if (a2 < a1)
        {
            a2++;
        }
        else
        {
            a1++;
            a2 = 0;
        }
    }
    std::ostringstream buff;
    buff << maxA2 << "/" << maxA1 << ":";
    for (size_t i = 0; i < ads.size() - 1; i++)
        buff << ads[i] << ",";
    buff << ads[ads.size() - 1] << ":";

    for (size_t i = 0; i < vas.size() - 1; i++)
        buff << vas[i] << ",";
    buff << vas[vas.size() - 1] << ":";

    for (size_t i = 0; i < gtLen; i++)
    {
        double lp = (lls[i] - maxP) / LG10;
        if (lp < LL_THRESHOLD)
            lp = LL_THRESHOLD;
        buff << int(-10 * lp);
        if (i != gtLen - 1)
            buff << ",";
    }
    gtString = buff.str();
}

// Input: variant and seqan::VarAlignInfo records for each read overlapping variant
// Output: Relative genotype likelihoods in log_2 scale and read info counts
inline void multiUpdateVC(seqan::VcfRecord const &          var,
                          std::vector<varAlignInfo> const & vais,
                          std::vector<double> &             vC,
                          std::vector<size_t> &             rI,
                          std::vector<size_t> &             VAs,
                          size_t const                      wSizeActual,
                          LRCOptions const &                O,
                          genotyping_model const            gtm)
{
    seqan::StringSet<seqan::CharString> altSet;
    strSplit(altSet, var.alt, seqan::EqualsChar<','>());

    size_t nAlts = length(altSet); // TODO: = var.nAlts

    std::vector<size_t> altLens;
    altLens.resize(nAlts);
    for (size_t i = 0; i < nAlts; i++)
    {
        altLens[i] = length(altSet[i]);
    }

    // Loop over all reads in bam file(s) 1 (deletion biased calls)
    for (auto & vai : vais)
    {
        std::vector<double> prefs;
        prefs.resize(nAlts + 1);
        // read does not occur in bam file(s) 2 (insertion biased calls)
        if (gtm == genotyping_model::ad || gtm == genotyping_model::joint)
        {
            size_t bestI = vai.alignmentPreference(wSizeActual, O, prefs);
            if (bestI != NO_BEST)
                rI[bestI]++;
            rI[rI.size() - 1]++;
        }

        if (gtm == genotyping_model::va || gtm == genotyping_model::joint)
        {
            size_t bestI = vai.vaPreference(O, length(var.ref), altLens, prefs);
            if (bestI != NO_BEST)
                VAs[bestI]++;
            VAs[VAs.size() - 1]++;
            if (O.verbose)
                std::cerr << "va " << /*TODO vai.first <<*/ " " << vai.nD << " " << vai.nI << " " << prefs[0] << " "
                          << prefs[1] << " " << bestI << '\n';
        }

        if (gtm == genotyping_model::va_old)
        {
            size_t bestI     = 0;
            double bestScore = 0; // std::numeric_limits< float>::max();
            for (size_t iP = 0; iP < nAlts; iP++)
            {
                double cScore =
                  O.overlapBits * (-vai.supports(length(var.ref), altLens[iP], O) * 1.0 +
                                   vai.rejects(length(var.ref), altLens[iP], O) * 1.0); // ACHTUNG: needs fixing
                prefs[iP + 1] += cScore;
                if (cScore < bestScore)
                {
                    bestScore = cScore;
                    bestI     = iP + 1;
                }
            }
            VAs[bestI]++;
            VAs[VAs.size() - 1]++;
            if (O.verbose)
                std::cerr << "va_old " << /*TODO vai.first <<*/ " " << vai.nD << " " << vai.nI << " " << prefs[0] << " "
                          << prefs[1] << " " << bestI << '\n';
        }

        if (gtm == genotyping_model::presence)
        {
            if (vai.present(O))
            {
                prefs[0] += O.overlapBits;
            }
            else
            {
                prefs[1] += O.overlapBits;
            }
            for (size_t iP = 2; iP <= nAlts; iP++)
            {
                prefs[iP] += O.overlapBits;
            }
        }

        //    if( O.verbose ) std::cerr << "Update VC bam1 " << vai.first << " " << prefs[0] << " " << prefs[1] << " "
        //    <<
        //    vai.supports( (int) length( var.ref ), (int) length( altSet[0] ), O ) << " " <<
        //    vai.rejects( (int) length( var.ref), (int) length( altSet[0]), O ) << '\n';

        float minPref = std::numeric_limits<float>::max();
        float maxPref = std::numeric_limits<float>::lowest();
        for (size_t iP = 0; iP < nAlts + 1; iP++)
        {
            if (prefs[iP] < minPref)
                minPref = prefs[iP];
            if (prefs[iP] > maxPref)
                maxPref = prefs[iP];
        }
        for (size_t iP = 0; iP < nAlts + 1; iP++)
            prefs[iP] -= minPref;

#define MINIMUM_PREF_DIFF 2.0

        if (maxPref - minPref > MINIMUM_PREF_DIFF)
        {
            size_t vCI = 0;
            for (size_t a1 = 0; a1 < nAlts + 1; a1++)
            {
                for (size_t a2 = 0; a2 <= a1; a2++)
                {
                    if (a1 != a2)
                    {
                        if (prefs[a1] == prefs[a2])
                        {
                            vC[vCI] += prefs[a1];
                        }
                        else if (prefs[a1] > prefs[a2] + 2)
                        {
                            vC[vCI] += prefs[a2] + 1;
                        }
                        else if (prefs[a2] > prefs[a1] + 2)
                        {
                            vC[vCI] += prefs[a1] + 1;
                        }
                        else if (prefs[a1] > prefs[a2])
                        {
                            vC[vCI] += (prefs[a1] + prefs[a2]) / 2.0;
                        }
                    }
                    else
                    {
                        vC[vCI] += prefs[a1];
                    }
                    vCI++;
                }
            }
        }
    }

    if (O.verbose)
        std::cerr << "multiUpdateVC " << vC[0] << " " << vC[1] << " " << vC[2] << '\n';
}

// Input: seqan::VcfRecord, reference and alternate reference
// Output: Chrom and position of variant, ref sequence and alt sequence
// TODO: Change this for a library that does this
inline void getLocRefAlt(seqan::VcfRecord const &  variant,
                         seqan::FaiIndex const &   faiI,
                         seqan::CharString const & chrom,
                         TSequence &               refSeq,
                         std::vector<TSequence> &  altSeqs,
                         size_t const              wSizeActual,
                         LRCOptions const &        O)
{
    TSequence                           ref = variant.ref;
    seqan::StringSet<seqan::CharString> altSet;
    strSplit(altSet, variant.alt, seqan::EqualsChar<','>());

    size_t nAlts = length(altSet);

    int32_t  beginPos = variant.beginPos;
    unsigned idx      = 0;

    if (!getIdByName(idx, faiI, chrom))
    {
        if (O.verbose)
            std::cerr << "rID " << chrom << " " << beginPos
                      << " WARNING: reference FAI index has no entry for rID in Ref std::mapped.\n";
    }

    if (O.genotypeRightBreakpoint)
        readRegion(refSeq, faiI, idx, beginPos - wSizeActual + length(ref), beginPos + length(ref) + wSizeActual);
    else
        readRegion(refSeq, faiI, idx, beginPos - wSizeActual, beginPos + wSizeActual);

    if (O.verbose)
        std::cerr << "refSeq " << refSeq << " " << chrom << " " << beginPos << std::endl;

    std::vector<size_t> altLens;
    altLens.resize(nAlts);
    for (size_t i = 0; i < nAlts; i++)
        altLens[i] = length(altSet[i]);

    size_t refLen = length(variant.ref);
    for (size_t i = 0; i < nAlts; i++)
    {
        if (!O.genotypeRightBreakpoint)
        {
            readRegion(altSeqs[i],
                       faiI,
                       idx,
                       beginPos - wSizeActual,
                       beginPos); // beginPos is included in the alt sequence
            TSequence post;
            if (altLens[i] < (size_t)wSizeActual)
            {
                post = altSet[i];
                TSequence post2;
                readRegion(post2, faiI, idx, beginPos + refLen, beginPos + refLen + wSizeActual - altLens[i]);
                append(post, post2);
            }
            else
            {
                post = infixWithLength(altSet[i], 0, wSizeActual);
            }
            append(altSeqs[i], post);
        }
        else
        {
            if (altLens[i] < (size_t)wSizeActual)
            {
                readRegion(altSeqs[i], faiI, idx, beginPos - wSizeActual + altLens[i], beginPos);
                append(altSeqs[i], altSet[i]);
            }
            else
            {
                altSeqs[i] = infixWithLength(altSet[i], altLens[i] - wSizeActual, wSizeActual);
            }
            TSequence post;
            readRegion(post, faiI, idx, beginPos + refLen, beginPos + refLen + wSizeActual);
            append(altSeqs[i], post);
        }
    }

    if (O.verbose)
    {
        std::cerr << "Printing altSeq " << '\n';
        for (size_t i = 0; i < nAlts; i++)
            std::cerr << "altSeq " << i << " " << altSeqs[i] << '\n';
        std::cerr << "Done printing altSeq " << '\n';
    }
}

// Crops a subsequence in a seqan::BamAlignmentRecord and rewrites into the seqan::BamAlignmentRecord
inline void cropSeq(seqan::BamAlignmentRecord const & bar,
                    seqan::VcfRecord const &          var,
                    ssize_t const                     wSizeActual, // <- this is signed here!
                    LRCOptions const &                O,
                    TSequence &                       croppedSeq)
{
    auto &  cigarString    = bar.cigar;
    ssize_t alignPos       = bar.beginPos;
    ssize_t readPos        = 0;
    ssize_t lReadPos       = 0;
    size_t  cigarI         = 0;
    char    cigarOperation = cigarString[cigarI].operation;

    // Searches for the first position overlapping our window (right insert) or last position overlapping window (left
    // insert)
    // TODO the following can underflow :o
    ssize_t searchPos = var.beginPos - wSizeActual; // Want to change the search intervals for TRs
    if (O.genotypeRightBreakpoint)
        searchPos = var.beginPos + length(var.ref) + wSizeActual;

    searchPos = std::max<ssize_t>(searchPos, 0);

    while (alignPos < searchPos && cigarI < length(cigarString))
    {
        // lAlignPos = alignPos;
        lReadPos       = readPos;
        cigarOperation = cigarString[cigarI].operation;

        switch (cigarOperation)
        {
            case 'D':
                alignPos += cigarString[cigarI].count;
                break;
            case 'M':
                alignPos += cigarString[cigarI].count;
                [[fallthrough]];
            case 'S':
            case 'H': // TODO THIS IS PROBABLY WRONG
            case 'I':
                readPos += cigarString[cigarI].count;
                break;
            default:
                std::cerr << "WARNING: cigar string case not accounted for \n";
        }

        if (O.verbose)
            std::cerr << bar.qName << " readpos " << alignPos << " " << var.beginPos << " " << searchPos << " "
                      << cigarI << " " << readPos << " " << cigarString[cigarI].count << " " << cigarOperation << '\n';
        cigarI++;
    }
    if (alignPos < searchPos && O.verbose)
    {
        std::cerr << "Read clipped " << alignPos << " " << var.beginPos << " " << searchPos << " " << cigarI << " "
                  << length(cigarString) << '\n';
    }

    if (cigarOperation == 'S' || cigarOperation == 'H')
        readPos = lReadPos;

    ssize_t rBeg = 0;
    ssize_t rEnd = 0;
    if (O.genotypeRightBreakpoint)
    {
        //    if( alignPos >= searchPos ){
        //   rBeg = readPos - 2*wSizeActual;
        // rEnd = readPos;
        //}else
        if (alignPos >= searchPos - 2 * wSizeActual)
        {
            ssize_t rShift = searchPos - alignPos;
            rBeg           = readPos - 2 * wSizeActual + rShift;
            rEnd           = readPos + rShift;
        }
        else
        {
            rBeg = readPos;
            rEnd = readPos + wSizeActual;
            if (O.verbose)
                std::cerr << "Insensible case for read " << bar.qName << " " << alignPos << " " << var.beginPos << " "
                          << searchPos << " " << cigarI << " " << length(cigarString) << '\n';
        }
    }
    else
    {
        ssize_t rShift = alignPos - searchPos;
        rBeg           = readPos - rShift;
        rEnd           = readPos + 2 * wSizeActual - rShift;
        if (rShift < 0 && O.verbose)
            std::cerr << "Poorly formatted read, case not accounted for " << bar.qName << " " << alignPos << " "
                      << var.beginPos << " " << searchPos << " " << cigarI << " " << length(cigarString) << '\n';
    }
    if (O.verbose)
        std::cerr << "Cropped read " << bar.qName << " " << alignPos << " " << var.beginPos << " " << searchPos << " "
                  << cigarI << " " << length(cigarString) << " " << rBeg << " " << rEnd << '\n';
    if (rBeg < 0)
        rBeg = 0;
    if (rEnd < 2 * wSizeActual)
        rEnd = 2 * wSizeActual;
    if (rEnd > (ssize_t)length(bar.seq))
        rEnd = (ssize_t)length(bar.seq);
    if (O.verbose)
        std::cerr << "ToInfix " << rBeg << " " << rEnd << " " << length(bar.seq) << '\n';

    croppedSeq = infixWithLength(bar.seq, rBeg, rEnd - rBeg);

    if (O.verbose)
        std::cerr << "Successful crop " << bar.qName << " " << croppedSeq << '\n';
}

template <typename TSeq>
inline TSeq mask(TSeq const & in)
{
    TSeq ret;
    seqan::appendValue(ret, in[0]);

    for (size_t i = 1; i < seqan::length(in); ++i)
        if (in[i] != in[i - 1])
            seqan::appendValue(ret, in[i]);

    return ret;
}

/** Input: bamStream, a VCF entry, reference fasta and alt fasta file if required by VCF entry
    Output: variant alignment info for each read near the VCF entry
 */
inline void LRprocessReads(seqan::VcfRecord const &                               variant,
                           seqan::CharString const &                              chrom,
                           seqan::FaiIndex const &                                faiI,
                           std::vector<seqan::BamAlignmentRecord const *> const & overlappingBars,
                           std::vector<varAlignInfo> &                            vais,
                           size_t const                                           wSizeActual,
                           LRCOptions const &                                     O)
{
    TSequence              refSeq;
    std::vector<TSequence> altSeqs;

    seqan::StringSet<seqan::CharString> altSet;
    strSplit(altSet, variant.alt, seqan::EqualsChar<','>());

    size_t nAlts = length(altSet);
    altSeqs.resize(nAlts);
    if (O.verbose)
        std::cerr << "nAlts " << nAlts << '\n';

    // move outside of this function
    getLocRefAlt(variant, faiI, chrom, refSeq, altSeqs, wSizeActual, O);

    if (O.outputRefAlt)
    {
        std::cerr << chrom << " " << /*beginPos +*/ 1 << " " << variant.info << " " << refSeq;
        for (size_t i = 0; i < nAlts; i++)
            std::cerr << " " << altSeqs[i];
        std::cerr << '\n';
        return;
    }

    if (O.mask)
        refSeq = mask(refSeq);

    seqan::Score<int, seqan::Simple> scoringScheme(O.match, O.mismatch, O.gapExtend, O.gapOpen);
    double                           band_fac = std::min<double>(O.bandedAlignmentPercent, 100.0) / 100.0;

    TSequence seqToAlign;
    for (size_t i = 0; i < overlappingBars.size(); ++i)
    {
        seqan::BamAlignmentRecord const & b   = *overlappingBars[i];
        varAlignInfo &                    vai = vais[i];

        clear(seqToAlign);

        if (O.cropRead)
        {
            cropSeq(b, variant, wSizeActual, O, seqToAlign);
        }
        else
        {
            seqToAlign = b.seq; // TODO: we shouldn't copy here
        }

        if (O.mask)
            seqToAlign = mask(seqToAlign);

        int lband = static_cast<double>(seqan::length(refSeq)) * band_fac;
        int rband = static_cast<double>(seqan::length(seqToAlign)) * band_fac;

        vai.alignS[0] = (int)localAlignmentScore(refSeq, seqToAlign, scoringScheme, -lband, +rband);

        for (size_t iAlt = 0; iAlt < nAlts; iAlt++)
        {
            if (O.mask)
            {
                auto al              = mask(altSeqs[iAlt]);
                lband                = static_cast<double>(seqan::length(al)) * band_fac;
                vai.alignS[iAlt + 1] = (int)localAlignmentScore(al, seqToAlign, scoringScheme, -lband, +rband);
            }
            else
            {
                lband = static_cast<double>(seqan::length(altSeqs[iAlt])) * band_fac;
                vai.alignS[iAlt + 1] =
                  (int)localAlignmentScore(altSeqs[iAlt], seqToAlign, scoringScheme, -lband, +rband);
            }
            if (O.verbose)
                std::cerr << "Aligned " << b.qName << " " << vai.alignS[0] << " " << vai.alignS[iAlt + 1] << " "
                          << altSeqs[iAlt] << '\n';
        }
    }
}

inline void initializeBam(std::string fileName, seqan::BamIndex<seqan::Bai> & bamIndex, seqan::BamFileIn & bamStream)
{
    if (!seqan::open(bamStream, fileName.data()))
        throw error{"Could not open ", fileName, " for reading."};

    fileName += ".bai";
    if (!seqan::open(bamIndex, fileName.data()))
        throw error{"Could not read BAI index file ", fileName};

    seqan::BamHeader header;
    readHeader(header, bamStream);
}

// Open a bam file or a set of bam files if the filename does not end with .bam
inline void parseBamFileName(std::filesystem::path const &              bfN,
                             std::vector<seqan::BamFileIn> &            bamStreamV,
                             std::vector<seqan::BamIndex<seqan::Bai>> & bamIndexV,
                             LRCOptions const &                         O)
{
    std::vector<std::filesystem::path> paths;

    if (bfN.native().ends_with(".bam") || bfN.native().ends_with(".sam.gz"))
    {
        paths.push_back(bfN);
    }
    else
    {
        std::ifstream fS;
        fS.open(bfN);
        std::string bf;
        while (fS >> bf)
            paths.push_back(bf);
    }

    if (O.verbose)
        std::cerr << "Checking input files" << (O.cacheDataInTmp ? " and copying to cache dir..." : "...");

    for (std::filesystem::path & p : paths)
    {
        if (!p.native().ends_with(".bam") && !p.native().ends_with(".sam.gz"))
            throw error{"Input file '", p, "' has unrecognized extension."};

        if (!std::filesystem::exists(p))
            throw error{"Input file '", p, "' does not exist."};

        std::filesystem::path p_bai = p;
        p_bai += ".bai";
        if (!std::filesystem::exists(p_bai))
            throw error{"Input file '", p, "' has no corresponding '.bai' index."};

        if (O.cacheDataInTmp) // copy to tmp
        {
            std::filesystem::path new_p     = O.cacheDir / p.filename();
            std::filesystem::path new_p_bai = O.cacheDir / p_bai.filename();

            if (std::filesystem::exists(new_p) || std::filesystem::exists(new_p_bai))
                throw error{"Cache file already exists. Does a filename appear twice in input?"};

            std::filesystem::copy(p, new_p);
            std::filesystem::copy(p_bai, new_p_bai);

            p = new_p; // update path in-place
        }
    }

    if (O.verbose)
        std::cerr << " done.";

    bamIndexV.resize(paths.size());
    bamStreamV.resize(paths.size());

    for (size_t i = 0; i < paths.size(); ++i)
        initializeBam(paths[i], bamIndexV[i], bamStreamV[i]);
}

// Examines a seqan::BamAlignmentRecord for evidence of supporting a variant and writes evidence into varAlignInfo
// record
inline void examineBamAlignment(seqan::BamAlignmentRecord const & bar,
                                seqan::VcfRecord const &          var,
                                varAlignInfo &                    vai,
                                LRCOptions const &                O)
{
    vai.reset();

    auto &  cigarString    = bar.cigar;
    int32_t alignPos       = bar.beginPos;
    size_t  cigarI         = 0;
    char    cigarOperation = cigarString[cigarI].operation;
    int32_t regionBeg      = var.beginPos - O.varWindow;
    int32_t regionEnd      = var.beginPos + (int32_t)length(var.ref) + O.varWindow;

    seqan::StringSet<seqan::CharString> infos;
    strSplit(infos, var.info, seqan::EqualsChar<';'>());

    for (size_t i = 0; i < length(infos); i++)
    {
        seqan::StringSet<seqan::CharString> info2;
        strSplit(info2, infos[i], seqan::EqualsChar<'='>());
        int32_t cVal;
        if (info2[0] == "TRRBEGIN")
        {
            if (info2[1] != ".")
            {
                cVal = std::atoi(toCString(info2[1])) - O.varWindow;
                if (cVal < regionBeg)
                    regionBeg = cVal;
            }
        }
        if (info2[0] == "TRREND")
        {
            if (info2[1] != ".")
            {
                cVal = std::atoi(toCString(info2[1])) + O.varWindow;
                if (cVal > regionEnd)
                    regionEnd = cVal;
            }
        }
        if (info2[0] == "REGBEGIN")
        {
            if (info2[1] != ".")
            {
                cVal = std::atoi(toCString(info2[1])) - O.varWindow;
                if (cVal < regionBeg)
                    regionBeg = cVal;
            }
        }
        if (info2[0] == "REGEND")
        {
            if (info2[1] != ".")
            {
                cVal = std::atoi(toCString(info2[1])) + O.varWindow;
                if (cVal > regionEnd)
                    regionEnd = cVal;
            }
        }
    }
    /*  if( infos.count( "TRRBEGIN" ) == 1 and infos.count( "TRREND" ) == 1 ){
      regionBeg = std::atoi(var.infos["TRRBEGIN"] ) - O.varWindow;
      seqan::CharString endStr = getValueById( var.genotypeInfos, trrEndID);
      regionEnd = std::atoi( var.infos["TRREND"] ) + O.varWindow;
      if( O.verbose ) std::cerr << "Found TRR " << regionBeg << " " << regionEnd << '\n';
  }else{*/
    if (O.verbose)
        std::cerr << "TRR " << regionBeg << " " << regionEnd << '\n';
    //  }

    if (alignPos < regionBeg)
        vai.alignsLeft = true;

    // Find the first position that overlaps the window we are interested in
    while (alignPos < regionBeg && cigarI < length(cigarString))
    {
        cigarOperation = cigarString[cigarI].operation;
        if (cigarOperation == 'M' || cigarOperation == 'D')
        {
            alignPos += cigarString[cigarI].count;
        }
        cigarI++;
    }

    // only counts the number of deleted bp in the window, probably fine to count a longer distance
    if (alignPos > regionBeg && cigarOperation == 'D' && (alignPos - (regionBeg) >= (int)O.minDelIns))
    {
        vai.nD = alignPos - regionBeg;
    }

    while (alignPos < regionEnd && cigarI < length(cigarString))
    {
        cigarOperation = cigarString[cigarI].operation;
        switch (cigarOperation)
        {
            case 'D':
                if (cigarString[cigarI].count >= O.minDelIns)
                    vai.nD += cigarString[cigarI].count;
                [[fallthrough]];
            case 'M':
                alignPos += cigarString[cigarI].count;
                break;
            case 'I':
                if (cigarString[cigarI].count >= O.minDelIns)
                    vai.nI += cigarString[cigarI].count;
                break;
            case 'S':
            case 'H':
                if (cigarString[cigarI].count > O.maxSoftClipped)
                {
                    if (!O.genotypeRightBreakpoint)
                    {
                        if (cigarI == length(cigarString) - 1)
                            vai.softClipped = true;
                    }
                    else
                    {
                        if (cigarI == 0)
                            vai.softClipped = true;
                    }
                }
                break;
            default:
                std::cerr << "WARNING: cigar string case not accounted for \n";
        }

        if (O.verbose)
            std::cerr << bar.qName << " " << cigarOperation << " " << cigarString[cigarI].count << " " << alignPos
                      << " " << cigarI << " " << vai.nD << " " << vai.nI << '\n';
        cigarI++;
    }
    if (alignPos > regionEnd)
        vai.alignsRight = true;

    if (O.verbose)
        std::cerr << "examinSeq " << bar.qName << " " << vai.nD << " " << vai.nI << " " << vai.softClipped << '\n';
}

// Gets reads in the region overlapping the variant
inline void parseReads(std::vector<seqan::BamAlignmentRecord> const &   bars,
                       seqan::VcfRecord const &                         var,
                       std::vector<seqan::BamAlignmentRecord const *> & overlappingBars,
                       std::vector<varAlignInfo> &                      align_infos,
                       size_t const                                     wSizeActual,
                       LRCOptions const &                               O)
{
    int32_t beg = var.beginPos - wSizeActual;
    int32_t end = var.beginPos + wSizeActual;

    if (O.genotypeRightBreakpoint)
    {
        beg = beg + length(var.ref);
        end = end + length(var.ref);
    }

    seqan::StringSet<seqan::CharString> altSet;
    strSplit(altSet, var.alt, seqan::EqualsChar<','>());

    size_t       nAlts = length(altSet);
    varAlignInfo vai(nAlts + 1);

    int32_t stopReading = beg;
    if (O.genotypeRightBreakpoint)
        stopReading = end;

    // TODO evaluate unordered_map here
    std::map<std::string_view, size_t> nameCache;

    for (seqan::BamAlignmentRecord const & record : bars)
    {
        if (overlappingBars.size() >= O.maxBARcount || record.beginPos > stopReading) // || record.rID != rID)
            return;

        // Ignore the read if it does not stretch to the region we are interested in
        if ((record.beginPos + (int32_t)length(record.seq) < beg) ||
            (record.beginPos + (int32_t)getAlignmentLengthInRef(record) < beg) || record.mapQ < O.minMapQ)
            continue;

        examineBamAlignment(record, var, vai, O);

        if (O.verbose)
            std::cerr << "Read record " << record.qName << '\n';

        // If we are on the next reference or at the end already then we stop.
        if (/*record.rID == -1 || record.rID > rID || */ record.beginPos >= end)
            break;

        // If we are left of the selected position then we skip this record.
        bool softClipRemove = false;
        if (not O.genotypeRightBreakpoint)
        {
            char cigarOperation = record.cigar[0].operation;
            if (cigarOperation == 'S' && (record.cigar[0].count > O.maxSoftClipped))
            {
                softClipRemove = true;
                if (O.verbose)
                    std::cerr << "SoftClip removed LeftBreakpoint " << record.qName << " " << O.genotypeRightBreakpoint
                              << " " << record.cigar[0].count << " " << O.maxSoftClipped << '\n';
            }
        }
        else if (O.genotypeRightBreakpoint)
        {
            char cigarOperation = record.cigar[length(record.cigar) - 1].operation;
            if (cigarOperation == 'S' && (record.cigar[length(record.cigar) - 1].count > O.maxSoftClipped))
            {
                softClipRemove = true;
                if (O.verbose)
                    std::cerr << "SoftClip removed RightBreakpoint " << record.qName << " " << O.genotypeRightBreakpoint
                              << " " << record.cigar[length(record.cigar) - 1].count << " " << O.maxSoftClipped << '\n';
            }
        }

        bool hardClipped     = false;
        char cigarOperationL = record.cigar[0].operation;
        char cigarOperationR = record.cigar[length(record.cigar) - 1].operation;

        if ((cigarOperationL == 'H') || (cigarOperationR == 'H'))
        {
            hardClipped = true;
            if (O.verbose)
                std::cerr << "Read " << record.qName << " is hardclipped at " << record.beginPos << '\n';
        }

        if ((!softClipRemove) && (!hasFlagDuplicate(record)) && (!hasFlagQCNoPass(record)) && (!hardClipped))
        {
            // prevent multiple alignments of the same read from being used
            std::string_view id = seqan::toCString(record.qName);
            if (nameCache.contains(id)) // replace existing
            {
                // TODO possibly check here which record is primary record and use that
                size_t index           = nameCache[id];
                overlappingBars[index] = &record;
                align_infos[index]     = vai;
            }
            else
            {
                nameCache[id] = overlappingBars.size();
                overlappingBars.push_back(&record);
                align_infos.push_back(vai); // this needs to be copied and not moved, because it is reused
            }
        }

        if (O.verbose)
            std::cerr << "Finished soft clipping " << '\n';
    }

    if (O.verbose)
        std::cerr << "Exiting readBamRegion " << '\n';
}

inline size_t getWSizeActual(std::span<seqan::VcfRecord> vcfRecords, LRCOptions const & O)
{
    if (O.dynamicWSize)
    {
        size_t maxAlleleLength = 0;
        for (seqan::VcfRecord & var : vcfRecords)
        {
            for (size_t i = 0, len = 1; i <= seqan::length(var.alt); ++i, ++len)
            {
                if (i == seqan::length(var.alt) || var.alt[i] == ',')
                {
                    maxAlleleLength = std::max(maxAlleleLength, len - 1);
                    len             = 0;
                }
            }
        }

        return maxAlleleLength + O.wSize;
    }
    else
    {
        return O.wSize;
    }
}

inline void processChunk(std::vector<seqan::BamFileIn> &            bamFiles,
                         std::vector<seqan::BamIndex<seqan::Bai>> & bamIndexes,
                         seqan::FaiIndex &                          faIndex,
                         seqan::CharString const &                  chrom,
                         std::vector<seqan::BamAlignmentRecord> &   bars,
                         std::span<seqan::VcfRecord>                vcfRecords,
                         LRCOptions const &                         O)
{
    size_t const wSizeActual = getWSizeActual(vcfRecords, O);

    /* Determine chromosome interval to fetch records for */
    size_t genome_begin = vcfRecords.front().beginPos;
    size_t genome_end   = vcfRecords.back().beginPos + 1;

    if (O.genotypeRightBreakpoint)
    {
        size_t minVarRef = std::numeric_limits<size_t>::max();
        size_t maxVarRef = std::numeric_limits<size_t>::min();

        for (seqan::VcfRecord & var : vcfRecords)
        {
            minVarRef = std::min<size_t>(minVarRef, seqan::length(var.ref));
            maxVarRef = std::max<size_t>(maxVarRef, seqan::length(var.ref));
        }

        genome_begin += minVarRef;
        genome_end += maxVarRef;
    }

    genome_begin = wSizeActual >= genome_begin ? 1 : genome_begin - wSizeActual;
    genome_end += wSizeActual;

    /* read BAM files for this chunk */
    for (size_t i = 0; i < bamFiles.size(); ++i)
    {
        size_t bamRID = 0;
        if (seqan::getIdByName(bamRID, seqan::contigNamesCache(seqan::context(bamFiles[i])), chrom))
            viewRecords(bars, bamFiles[i], bamIndexes[i], bamRID, genome_begin, genome_end);

        // else: BAM files that have no reads spanning the desired chromosome are quietly ignored
    }

    if (bamFiles.size() > 1)
    {
        std::ranges::sort(bars,
                          [](seqan::BamAlignmentRecord const & lhs, seqan::BamAlignmentRecord const & rhs)
                          { return lhs.beginPos < rhs.beginPos; });
    }

    /* process variants */
    for (seqan::VcfRecord & var : vcfRecords)
    {
        size_t nAlleles = 2;
        for (char c : var.alt)
            if (c == ',')
                ++nAlleles;

        // Implemented as a std::vector as we may have more than one output per marker
        std::vector<std::vector<double>> vC; // variant calls
        std::vector<std::vector<size_t>> AD; // Allele depth counts
        std::vector<std::vector<size_t>> VA; // Variant seqan::Alignment counts

        if (O.gtModel == genotyping_model::multi)
        {
            vC.resize(5);
            AD.resize(5);
            VA.resize(5);
        }
        else
        {
            vC.resize(1);
            AD.resize(1);
            VA.resize(1);
        }

        for (size_t mI = 0; mI < vC.size(); mI++)
        {
            vC[mI].resize(nAlleles * (nAlleles + 1) / 2);
            AD[mI].resize(nAlleles + 1);
            VA[mI].resize(nAlleles + 1);
        }

        std::vector<seqan::BamAlignmentRecord const *> overlappingBars;
        std::vector<varAlignInfo>                      alignInfos;
        parseReads(bars, var, overlappingBars, alignInfos, wSizeActual, O);

        //         for (seqan::BamAlignmentRecord const * bar : overlappingBars)
        //                 std::cerr << bar->qName << '\n';

        LRprocessReads(var, chrom, faIndex, overlappingBars, alignInfos, wSizeActual, O);

        //         if (useBam2)
        //         {
        //             for (size_t i = 0; i < bamIndex2Handles[omp_get_thread_num()].size(); i++)
        //             {
        //                 std::map<seqan::CharString, seqan::BamAlignmentRecord> bars;
        //                 readBamRegion(bamIndex2Handles[omp_get_thread_num()][i],
        //                               bam2Handles[omp_get_thread_num()][i],
        //                               varChrom,
        //                               var,
        //                               O,
        //                               bars,
        //
        //                               vaisR);
        //                 LRprocessReads(var, vcfIn, faIndexHandles[omp_get_thread_num()], bars, O, vaisR);
        //             }
        //         }

        if (O.gtModel == genotyping_model::multi)
        {
            multiUpdateVC(var, alignInfos, vC[0], AD[0], VA[0], wSizeActual, O, genotyping_model::ad);
            multiUpdateVC(var, alignInfos, vC[1], AD[1], VA[1], wSizeActual, O, genotyping_model::va);
            multiUpdateVC(var, alignInfos, vC[2], AD[2], VA[2], wSizeActual, O, genotyping_model::joint);
            multiUpdateVC(var, alignInfos, vC[3], AD[3], VA[3], wSizeActual, O, genotyping_model::presence);
            multiUpdateVC(var, alignInfos, vC[4], AD[4], VA[4], wSizeActual, O, genotyping_model::va_old);
        }
        else
        {
            multiUpdateVC(var, alignInfos, vC[0], AD[0], VA[0], wSizeActual, O, O.gtModel);
        }

        std::string gtString;
        for (size_t mI = 0; mI < vC.size(); mI++)
        {
            getGtString(vC[mI], AD[mI], VA[mI], gtString);
            appendValue(var.genotypeInfos, gtString);
            var.format = "GT:AD:VA:PL";
        }
    }
}
