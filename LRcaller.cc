#define SEQAN_HAS_ZLIB 1

#include <iostream>
#include <cassert>
#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <time.h>
#include <omp.h>
#include <seqan/align.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/bam_io.h>
#include <seqan/store.h>
#include <seqan/vcf_io.h>
#include <seqan/stream.h>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>


using namespace seqan;
using namespace std;


// Sequence, alignment, and alignment row.
typedef String< Dna5> TSequence;
typedef Align<TSequence, ArrayGaps> TAlign;
typedef Row<TAlign>::Type TRow;

#define LL_THRESHOLD -25.5
#define LG10 3.322
#define NO_ALIGNMENT -10000

enum aligner_type{ aligner_minimap, aligner_bwa, aligner_seqan };

enum genotyping_model{ genotype_multi, genotype_ad, genotype_va, genotype_joint};


/*  Option handling */
struct LRCOptions{
public:
  CharString bam;
  CharString bam2;
  CharString faFile;
  CharString vcfInFile;
  CharString vcfOutFile;
  int wSize;  // Window Size, 
  int maxBARcount;  // maximum number of reads to use for a variant
  bool verbose;
  bool genotypeRightBreakpoint;
  int minAlignScore;  // minimum alignment score to consider the alignment
  double logScaleFactor;
  CharString genotypingModel;

  int match;
  int mismatch;
  int gapOpen;
  int gapExtend;

  bool cropRead; 
  aligner_type aligner;
  int randID;  // To avoid file name collision when runing bwa/minimap
  bool outputRefAlt;  // For debug purposes
  int maxSoftClipped; // Read cannot be softclipped at invariant end
  double maxAlignBits;  //needs to be added as option
  int overlapBits;  // needs to be added as option
  int varWindow;  //How far 
  int maxReadLenForCropping;  // Should not be needed
  int minMapQ; // minimum mapQ for a read to be used
  int minDelIns; // minum number of consecutive bp to be deleted or inserted to count it
  double refThreshFraction;  // Fraction of deletion/insertion length that can be deleted/inserted but still called ref
  double altThreshFraction;  // Fraction of deletion/insertion length that needs to be deleted/inserted to call alt
  genotyping_model gtModel;
  int nThreads;
  LRCOptions():  bam2(""),faFile("genome.fa"),wSize(500), maxBARcount(200),verbose(false), genotypeRightBreakpoint(false), minAlignScore(400),logScaleFactor(2.0),genotypingModel("joint"),match(1),mismatch(-1),gapOpen(-1),gapExtend(-1),cropRead(true),aligner(aligner_seqan),outputRefAlt(false), maxSoftClipped(500), maxAlignBits( 10.0), overlapBits(10), varWindow(100), maxReadLenForCropping(50000), minMapQ(30), minDelIns(6), refThreshFraction(1.0/10.0), altThreshFraction(1.0/2.0), gtModel(genotype_joint), nThreads(1){};
  //regionWindowSize - how far we want to look in the bam file
};

template<typename TOptions>
int parseLRCArguments( int argc, char const ** argv, TOptions & O ){

  ArgumentParser parser("LRcaller");
  addUsageLine(parser, "[\\fIOPTIONS\\fP]  \"\\fIBAMFILE\\fP\"  \"\\fIVCF_FILE_IN\\fP\" \"\\fIVCF_OUT_FILE\\fP\" ");
  addDescription(parser, "Genotypes variants using long reads\n in bam file/file of bam files ");
 
  addArgument(parser, ArgParseArgument(ArgParseArgument::STRING, "BAMFILE bam file/file of bam files"));
  addArgument(parser, ArgParseArgument(ArgParseArgument::STRING, "VCF_FILE_IN - input vcf file"));
  addArgument(parser, ArgParseArgument(ArgParseArgument::STRING, "VCF_FILE_OUT - genotyped vcf file"));
  addOption(parser, ArgParseOption("fa", "fa", "Fastafile", ArgParseArgument::STRING, "FA"));
  addOption(parser, ArgParseOption("b2", "bam2", "Second bam file/file of bam files, use if you want to use multiple basecallings/mappings for the same reads", ArgParseArgument::STRING, "b2"));
 
  addOption(parser, ArgParseOption("a", "aligner", "Aligner used bwa/minimap/seqan (default seqan). (Using bwa or minimap requires programs to be in path)", ArgParseArgument::STRING, "ALIGNER"));
  addOption(parser, ArgParseOption("gtm", "genotyper", "Genotyper used joint/ad/va/multi (default joint).", ArgParseArgument::STRING, "GENOTYPER"));

  addOption(parser, ArgParseOption("faA", "faA", "Fasta file alt, used when VCF file refers to an alternate fasta file", ArgParseArgument::STRING, "FAA"));
  addOption(parser, ArgParseOption("v", "verbose", "Verbose output"));
  addOption(parser, ArgParseOption("rb", "right_breakpoint", "Genotype right breakpoint, left breakpoint genotyped by default"));
  addOption(parser, ArgParseOption("ora", "get_ref_alt", "Output ref and alt allele, for debugging purposes"));

  addOption(parser, ArgParseOption("A", "match", "Match score for alignment", ArgParseArgument::INTEGER, "INT"));
  addOption(parser, ArgParseOption("B", "mismatch", "Error penalty for alignment.", ArgParseArgument::INTEGER, "INT"));
  addOption(parser, ArgParseOption("O", "gap_open", "Gap open for alignment", ArgParseArgument::INTEGER, "INT"));
  addOption(parser, ArgParseOption("E", "gap_extend", "Gap extend for alignment", ArgParseArgument::INTEGER, "INT"));
  addOption(parser, ArgParseOption("w", "window_size", "Window size", ArgParseArgument::INTEGER, "INT"));
  setDefaultValue(parser, "A", O.match);    
  setDefaultValue(parser, "B", O.mismatch);
  setDefaultValue(parser, "O", O.gapOpen);
  setDefaultValue(parser, "E", O.gapExtend);
  setDefaultValue(parser, "w", O.wSize);

  addOption(parser, ArgParseOption("ab", "align_bits", "Max bits for alignment", ArgParseArgument::DOUBLE, "DOUBLE"));
  setDefaultValue(parser, "ab", O.maxAlignBits);    

  addOption(parser, ArgParseOption("ob", "overlap_bits", "Score for overlap (log_2)", ArgParseArgument::INTEGER, "INT"));
  setDefaultValue(parser, "ob", O.overlapBits);    

  addOption(parser, ArgParseOption("nt", "number_of_threads", "Number of threads", ArgParseArgument::INTEGER, "INT"));
  setDefaultValue(parser, "nt", O.nThreads);    

  addOption(parser, ArgParseOption("vw", "var_window", "Var window size, look for del/ins inside this window", ArgParseArgument::INTEGER, "INT"));
  setDefaultValue(parser, "vw", O.varWindow);    

  addOption(parser, ArgParseOption("mapQ", "read_mapq", "Minimum read mapQ", ArgParseArgument::INTEGER, "INT"));
  setDefaultValue(parser, "mapQ", O.minMapQ);    

  addOption(parser, ArgParseOption("di", "min_del_ins", "Minimum del/ins in cigar to consider", ArgParseArgument::INTEGER, "INT"));
  setDefaultValue(parser, "di", O.minDelIns);    

  addOption(parser, ArgParseOption("rt", "ref_thresh_fraction", "Threshold for fraction of del/ins bp to be considered ref", ArgParseArgument::DOUBLE, "DOUBLE"));
  setDefaultValue(parser, "rt", O.refThreshFraction);    

  addOption(parser, ArgParseOption("at", "alt_thresh_fraction", "Threshold for fraction of del/ins bp to be considered alt", ArgParseArgument::DOUBLE, "DOUBLE"));
  setDefaultValue(parser, "at", O.altThreshFraction);    

  addOption(parser, ArgParseOption("lsf", "logScaleFactor", "Log scale factor for comparing alignment scores", ArgParseArgument::DOUBLE, "DOUBLE"));

  setDefaultValue( parser, "fa", O.faFile );
  setDefaultValue( parser, "b2", O.bam2 );
  setValidValues(parser, "aligner", "bwa minimap seqan");
  setDefaultValue( parser, "aligner", "seqan" );
  setValidValues(parser, "genotyper", "joint ad va multi");
  setDefaultValue( parser, "genotyper", "joint" );
  setDefaultValue(parser, "logScaleFactor", O.logScaleFactor);
  ArgumentParser::ParseResult res = parse(parser, argc, argv);
  // Only extract  options if the program will continue after parseCommandLine()
  if (res != ArgumentParser::PARSE_OK)
    return res;
  getArgumentValue( O.bam, parser, 0);
  getArgumentValue( O.vcfInFile, parser, 1);
  getArgumentValue( O.vcfOutFile, parser, 2);

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

  if (isSet(parser, "ref_thresh_fraction"))
    getOptionValue(O.refThreshFraction, parser, "ref_thresh_fraction");

  if (isSet(parser, "alt_thresh_fraction"))
    getOptionValue(O.altThreshFraction, parser, "alt_thresh_fraction");

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
  CharString AlignerName;
  if (isSet(parser, "aligner"))
    getOptionValue(AlignerName, parser, "aligner");
  if( AlignerName == "bwa" )
    O.aligner = aligner_bwa;
  else if( AlignerName == "minimap" )
    O.aligner = aligner_minimap;
  else if( AlignerName == "seqan" )
    O.aligner = aligner_seqan;
  CharString GTmodelName;
  if (isSet(parser, "genotyper"))
    getOptionValue(GTmodelName, parser, "genotyper");
  if( GTmodelName == "multi" )
    O.gtModel = genotype_multi;
  else if( GTmodelName == "ad" )
    O.gtModel = genotype_ad;
  else if( GTmodelName == "va" )
    O.gtModel = genotype_va;
  else if( GTmodelName == "joint" )
    O.gtModel = genotype_joint;
  if (isSet(parser, "logScaleFactor"))
    getOptionValue(O.logScaleFactor, parser, "logScaleFactor");
  //O.cropRead = isSet(parser, "cropread");
  O.verbose = isSet(parser, "verbose");
  O.genotypeRightBreakpoint = isSet(parser, "right_breakpoint");
  O.outputRefAlt = isSet(parser, "get_ref_alt");
  //get options
  return res;
} 


/* Stores information of how a read aligns across a variant */
class varAlignInfo{
public:
  int nD;
  int nI;
  double refS;
  double altS;
  bool softClipped;
  bool alignsLeft;
  bool alignsRight;
  //int refAlignLen;

  // Supports and rejects logic is not correct for very long variants
  // Alternate allele is supported as judged from bam alignment record
  template<typename TOptions>
  bool supports( VcfRecord& var, TOptions& O ){
    if( length( var.alt ) > length(var.ref)){  // insertion (these are simplistic for insertion/deletion type variants)
      float threshSize = (float) length(var.alt);

      if( (alignsLeft and alignsRight and ((float) nI > (threshSize)*O.altThreshFraction)) or softClipped )
	return true;
      else
	return false;
    }else{
      if( (alignsLeft and alignsRight and ((float) nD > ((float) length(var.ref))*O.altThreshFraction)) or softClipped )
	return true;
      else
	return false;
    }
  }

  // Alternate allele is rejected as judged from bam alignment record
  template<typename TOptions>
  bool rejects( VcfRecord& var, TOptions& O  ){
    if( length( var.alt ) > length(var.ref)){  // insertion
      if( (alignsLeft and alignsRight and ((float) nI < ((float) length(var.alt))*O.refThreshFraction)) and (not softClipped) )
	return true;
      else
	return false;
    }else{
      if( (alignsLeft and alignsRight and ((float) nD < ((float) length(var.ref))*O.refThreshFraction)) and (not softClipped) )
	return true;
      else
	return false;
    }
  }
  bool aligns(){
    return alignsLeft and alignsRight;
  }
  void reset(){
    nD = 0;
    nI = 0;
    refS = NO_ALIGNMENT;
    altS = NO_ALIGNMENT;
    softClipped = false;
    alignsLeft = false;
    alignsRight = false;
  }

  // How much more likely the alternate is to the reference, reported on log_2 scale, i.e. a value 10, represents alt = 2^10=1024 times more likely than ref
  template<typename TOptions>
  double alignmentAltPreference( TOptions& O){
    double ret = 0;
    if( refS == NO_ALIGNMENT or refS <= O.minAlignScore ){
      if( altS == NO_ALIGNMENT or altS <= O.minAlignScore ){
	return 0;
      }else{
	ret = (altS - O.minAlignScore)/O.logScaleFactor;
      }
    }else{
      if( altS == NO_ALIGNMENT or altS <= O.minAlignScore ){
	ret = -(refS - O.minAlignScore)/O.logScaleFactor;
      }else{
	ret = (altS - refS)/O.logScaleFactor;
      }
    }
    if( ret < -O.maxAlignBits ) ret = -O.maxAlignBits;
    if( ret > O.maxAlignBits ) ret = O.maxAlignBits;
    return ret;
  }
  varAlignInfo(): nD(0),nI(0), refS(NO_ALIGNMENT),altS(NO_ALIGNMENT),softClipped(false),alignsLeft(false),alignsRight(false){};
};


/* Turns genotyping into string */
void
getGtString(std::vector<double> & lls, std::vector< int> & ads, std::vector< int>& vas, std::string & gtString )
{
  lls[0] = -lls[0];  // Really silly hack for historical reasons, would confuse the hell out of me to fix it
  lls[1] = -lls[1];
  lls[2] = -lls[2];
  int max = 0;
    double maxP = lls[0];
    if( lls[1] > maxP ){
        maxP = lls[1];
        max = 1;
    }
    if( lls[2] > maxP ){
        maxP = lls[2];
        max = 2;
    }
    std::ostringstream buff;
    if( max == 0 ){
        buff << "0/0:";
    }else if( max == 1 ){
        buff << "0/1:";
    }else{
        buff << "1/1:";
    }
    if( ads.size() == 2 )
      buff << ads[0] << "," << ads[1] << ":";
    else 
      buff << ads[0] << "," << ads[1] << "," << ads[2] << ":";
    if( vas.size() == 3 )
      buff << vas[0] << "," << vas[1] << "," << vas[2] << ":";

    double div = LG10;
    double lp0 = (lls[0]-maxP)/div;
    double lp1 = (lls[1]-maxP)/div;
    double lp2 = (lls[2]-maxP)/div;
    if( lp0 < LL_THRESHOLD ) lp0 = LL_THRESHOLD;
    if( lp1 < LL_THRESHOLD ) lp1 = LL_THRESHOLD;
    if( lp2 < LL_THRESHOLD ) lp2 = LL_THRESHOLD;

    int phr0 = int( -10*lp0 );
    int phr1 = int( -10*lp1 );
    int phr2 = int( -10*lp2 );
    buff  << phr0 << "," << phr1 << "," << phr2;
    gtString = buff.str();
}


//Input: VcfRecord, reference and alternate refernce
//Output: Chrom and position of variant, ref sequence and alt sequence
//TODO: Change this for a library that does this
template<typename TOptions>
int getLocRefAlt( VcfRecord & variant, VcfFileIn & vcfS,FaiIndex & faiI, 
		  CharString& chrom, int32_t& beginPos, TSequence& refSeq, TSequence& altSeq, TOptions & O )
{
    TSequence ref = variant.ref;
    CharString alt = variant.alt;

    int32_t rID = variant.rID;
    beginPos = variant.beginPos;
    unsigned idx = 0;
    chrom = contigNames(context(vcfS))[rID];
    if (!getIdByName(idx, faiI, chrom)){
        if( O.verbose ) cerr << "rID " << chrom << " " << beginPos << " ERROR: reference FAI index has no entry for rID in Ref mapped.\n";
    }


    if( O.genotypeRightBreakpoint ){
      readRegion( refSeq, faiI, idx, beginPos-O.wSize+length(ref), beginPos+length(ref)+O.wSize);
    }else{
      readRegion( refSeq, faiI, idx, beginPos-O.wSize, beginPos+O.wSize);
    }
    if( O.verbose ) cerr << "refSeq " << refSeq << " " << chrom << " " << beginPos << std::endl;

    int altLen = length( variant.alt );
    int refLen = length( variant.ref );
    if( not O.genotypeRightBreakpoint ){
      readRegion(altSeq, faiI, idx, beginPos-O.wSize, beginPos ); // beginPos is included in the alt sequence
      TSequence post;
      if( altLen < O.wSize ){
	post = variant.alt;
	TSequence post2;
	readRegion(post2, faiI, idx, beginPos+refLen, beginPos+refLen + O.wSize - altLen );
	append( post, post2);
      }else{
	post = infixWithLength( variant.alt, 0, O.wSize );
      }
      append( altSeq, post );
    }else{	
      if( altLen < O.wSize ){
	readRegion(altSeq, faiI, idx, beginPos-O.wSize+altLen, beginPos );
	append( altSeq, variant.alt );
      }else{
	altSeq = infixWithLength( variant.alt, altLen-O.wSize, O.wSize );
      }
      TSequence post;
      readRegion(post, faiI, idx, beginPos+refLen, beginPos + refLen+O.wSize );
      append( altSeq, post );
    }
    if( O.verbose ) cerr << "altSeq " << altSeq << endl;
    return 0;
}



// Crops a subsequence in a BamAlignmentRecord and rewrites into the BamAlignmentRecord
template<typename TOptions>
int cropSeq(  BamAlignmentRecord & bar, VcfRecord & var, TOptions& O, CharString& croppedSeq )
{
  String<CigarElement<> > cigarString = bar.cigar;
  int alignPos = bar.beginPos;
  // int lAlignPos = bar.beginPos;
  int readPos = 0;
  int lReadPos = 0;
  int cigarI = 0;
  CharString cigarOperation = cigarString[cigarI].operation;
  string cigarOperationStr = toCString(cigarOperation);
  
  // Searches for the first position overlapping our window (right insert) or last position overlapping window (left insert)
  int searchPos = var.beginPos-O.wSize;
  if( O.genotypeRightBreakpoint )
    searchPos = var.beginPos + (int) length( var.ref) + O.wSize;

  while( alignPos < searchPos and cigarI < (int) length( cigarString) ){
    //lAlignPos = alignPos;
    lReadPos = readPos;
    cigarOperation = cigarString[cigarI].operation;
    cigarOperationStr = toCString(cigarOperation);
    if( cigarOperationStr == "S" or cigarOperationStr == "H" ){
      readPos += cigarString[cigarI].count;
    }else if( cigarOperationStr == "M" ){
      readPos += cigarString[cigarI].count;
      alignPos += cigarString[cigarI].count;
    }else if( cigarOperationStr == "D" ){
      alignPos += cigarString[cigarI].count;
    }else if( cigarOperationStr == "I" ){
      readPos += cigarString[cigarI].count;
    }else{
      if( O.verbose ) cerr << "WTF cigar string case not accounted for " << endl;
    }
    cigarI++;
  }
  if( alignPos < searchPos and O.verbose ){
    cout  << "Read clipped " << alignPos << " " << var.beginPos << " " << searchPos << " " << cigarI << " " << length(cigarString) << endl;  
  }  

  if( cigarOperationStr == "S" or cigarOperationStr == "H") readPos = lReadPos;
  
  int rBeg,rEnd;
  if( O.genotypeRightBreakpoint ){
    if( alignPos >= searchPos ){
       rBeg = readPos - 2*O.wSize;
      rEnd = readPos;
    }else if( alignPos >= searchPos - 2*O.wSize ){
      int rShift = searchPos-alignPos;
      rBeg = readPos - 2*O.wSize+rShift;
      rEnd = readPos + rShift;
    }else{
      rBeg = readPos;
      rEnd = readPos+O.wSize;
      if( O.verbose ) cerr << "Insensible case for read " << bar.qName << " " << alignPos << " " << var.beginPos << " " << searchPos << " " << cigarI << " " << length(cigarString) << endl;
    }
  }else{
    int rShift = alignPos-searchPos;
    rBeg = readPos-rShift;
    rEnd = readPos + 2*O.wSize-rShift;
    if( rShift < 0 and O.verbose ) cout << "Poorly formatted read, case not accounted for " << bar.qName << " " << alignPos << " " << var.beginPos << " " << searchPos << " " << cigarI << " " << length(cigarString) << endl;
  }
  if( O.verbose ) cerr << "Cropped read " << bar.qName <<  " " << alignPos << " " << var.beginPos << " " << searchPos << " " << cigarI << " " << length(cigarString) << " " << rBeg << " " << rEnd << endl;
  if( rBeg < 0 ) rBeg = 0;
  if( rEnd < 2*O.wSize ) rEnd = 2*O.wSize;
  if( rEnd > (int) length( bar.seq ) ) rEnd = (int) length( bar.seq );
  if( O.verbose ) cerr << "ToInfix "<< rBeg << " " << rEnd << " " << length(bar.seq) << endl; 
  croppedSeq = infixWithLength( bar.seq, rBeg, rEnd-rBeg );
  if( O.verbose ) cerr << "Successful crop " << bar.qName << " " << croppedSeq << endl;
  return 0;
}


// Examines a BamAlignmentRecord for evidence of supporting a variant and writes evidence into varAlignInfo record 
template<typename TOptions>
int examineBamAlignment(  BamAlignmentRecord & bar, VcfRecord & var, varAlignInfo& vai, TOptions& O )
{
  String<CigarElement<> > cigarString = bar.cigar;
  int alignPos = bar.beginPos;
  int cigarI = 0;
  CharString cigarOperation = cigarString[cigarI].operation;
  string cigarOperationStr = toCString(cigarOperation);
  vai.reset();
  if( alignPos < var.beginPos - O.varWindow ) vai.alignsLeft = true;

  // Find the first position that overlaps the window we are interested in
  while( alignPos < var.beginPos-O.varWindow and cigarI < (int) length( cigarString) ){
    cigarOperation = cigarString[cigarI].operation;
    cigarOperationStr = toCString(cigarOperation);
    if( cigarOperationStr == "M" or cigarOperationStr == "D" ){
      alignPos += cigarString[cigarI].count;
    }
    cigarI++;
  }

  // only counts the number of deleted bp in the window, probably fine to count a longer distance
  if( alignPos > var.beginPos - O.varWindow and cigarOperationStr == "D" and (alignPos - (var.beginPos -O.varWindow) >= (int) O.minDelIns) ){
    vai.nD = alignPos - (var.beginPos - O.varWindow);
  }

  while( alignPos < var.beginPos+O.varWindow+ (int) length(var.ref) and cigarI < (int) length( cigarString) ){
    cigarOperation = cigarString[cigarI].operation;
    cigarOperationStr = toCString(cigarOperation);
    if( cigarOperationStr == "M" ){
      alignPos += cigarString[cigarI].count;
    }else if( cigarOperationStr == "D" ){
      if( (int) cigarString[cigarI].count >= O.minDelIns ) vai.nD += cigarString[cigarI].count;
      alignPos += cigarString[cigarI].count;
    }else if( cigarOperationStr == "I" ){
      if( (int) cigarString[cigarI].count >= O.minDelIns ) vai.nI += cigarString[cigarI].count;
    }else if( (cigarOperationStr == "S" or cigarOperationStr == "H") and (int) cigarString[cigarI].count > O.maxSoftClipped ){  
      // ACHTUNG, could be softclipped from either end, not clear if this is the best choice
      if( not O.genotypeRightBreakpoint ){
	if( cigarI == (int) length(cigarString) - 1 ){
	  vai.softClipped = true;
	}
      }else{
	if( cigarI == 0 ){
	  vai.softClipped = true;
	}
      }
    }
    if( O.verbose ) cerr << bar.qName << " " << cigarOperationStr << " " << cigarString[cigarI].count << " " << alignPos << " " << cigarI << " " << vai.nD << " " << vai.nI << endl;
    cigarI++;
   }  
  if( alignPos > var.beginPos + O.varWindow ) vai.alignsRight = true;
  if( O.verbose ) cerr << "examinSeq " << bar.qName << " " << vai.nD << " " << vai.nI << " " << vai.softClipped << endl;
  return 0;
}

// writes fasta file for BWA or minimap from bam alignment records
template<typename TOptions>
int writeFastaFromBars(   std::map< CharString, BamAlignmentRecord>& bars, int& maxBC, VcfRecord& var, TOptions& O )
{
  ofstream minimapF;
  stringstream ss;
  ss << "reads_to_minimap_" << O.randID << ".fa";
  string str = ss.str();
  minimapF.open( str );
  int barCount = 0;

  for( auto i = bars.begin(); i != bars.end() && barCount < maxBC; i++, barCount++ ){
    minimapF << ">" << (*i).second.qName << endl;
    if( O.cropRead ){
      CharString croppedSeq;
      cropSeq( (*i).second, var, O, croppedSeq);
      minimapF << croppedSeq << endl;
    }else{
      minimapF << (*i).second.seq << endl;
    }
  }
  minimapF.close();
  return 0;
}


int initializeBam(char* fileName, BamIndex<Bai> & bamIndex, BamFileIn & bamStream)
{
    if (!open(bamStream, fileName))
    {
        std::cerr << "ERROR: Could not open " << fileName << " for reading.\n";
        return 1;
    }

    char baifileName[strlen(fileName) + 10];
    strcpy(baifileName, fileName);
    strcat(baifileName, ".bai");
    if (!open(bamIndex, baifileName))
    {
        std::cerr << "ERROR: Could not read BAI index file " << fileName << "\n";
        return 1;
    }

    BamHeader header;
    readHeader(header, bamStream);
    return 0;
}


//Gets reads in the region overlapping the variant
template<typename TOptions>
int readBamRegion(BamIndex< Bai>& baiI, BamFileIn& bamS, CharString& chrom, VcfRecord& var,
		  TOptions& O, std::map< CharString, BamAlignmentRecord>& bars, std::map< CharString, varAlignInfo>& vaiMap )
{
  if (atEnd(bamS))
    return 0;

  int beg = var.beginPos-O.wSize;
  int end = var.beginPos+O.wSize; 

  if( O.genotypeRightBreakpoint ){
    beg = beg+length( var.ref );
    end = end+length( var.ref );
  }

  if( O.verbose ) cerr << "reading Bam region " << chrom << " " << var.beginPos << " "  << beg << " " << end << std::endl;
  
  int rID = 0;
  if (!getIdByName( rID, contigNamesCache(context( bamS )), chrom )){
    if( O.verbose ) cerr << "ERROR: Reference sequence named " << chrom << " not found in bam file." << std::endl;
    return 1;
  }

  BamAlignmentRecord record;
  varAlignInfo vai;

  
  if( O.verbose ) cerr << "reading " << beg << " " << record.seq << endl;

  //This should be a better implementation, but it is 2x slower
  /*  std::vector<seqan::BamAlignmentRecord> reads;
  viewRecords(reads, bamS, baiI, rID, beg, end);
  cout << "reads size " << reads.size() << endl;*/


  bool hasAlignments = false;
  unsigned regionBeg = std::max(1, beg - O.maxReadLenForCropping); // make sure we have the best implementation of this!!
  if (!jumpToRegion( bamS, hasAlignments, rID, regionBeg, end, baiI )){
    std::cerr << "ERROR: Could not jump to " << rID << " " << chrom << ":" << beg << "-" << end << "\n";
    return 1;
  }
  int stopReading = beg;
  if( O.genotypeRightBreakpoint ) stopReading = end;
 
  while (!atEnd(bamS)){
    readRecord(record, bamS);  
    if( (int) bars.size() >= O.maxBARcount or record.beginPos > stopReading or record.rID != rID )
      return 0;
    // Ignore the read if it does not stretch to the region we are interested in
    if( (int) (record.beginPos + length(record.seq)) < beg or (int) (record.beginPos + getAlignmentLengthInRef(record))  < beg or record.mapQ < O.minMapQ )
      continue;
    examineBamAlignment( record, var, vai, O );

    if( O.verbose ) cerr << "Read record " << record.qName << endl;

    // If we are on the next reference or at the end already then we stop.
    if (record.rID == -1 || record.rID > rID || record.beginPos >= end )
      break;
    // If we are left of the selected position then we skip this record.

    bool softClipRemove = false;
    if( not O.genotypeRightBreakpoint ){
      CharString cigarOperation = record.cigar[0].operation;
      string cigarOperationStr = toCString(cigarOperation);
      if( cigarOperationStr == "S" and  ( (int) record.cigar[0].count > O.maxSoftClipped) ){
	softClipRemove = true;
	if( O.verbose ) cerr << "SoftClip removed LeftBreakpoint " << record.qName << " " << O.genotypeRightBreakpoint << " " << record.cigar[0].count << " " << O.maxSoftClipped << endl;
      }
    }else if( O.genotypeRightBreakpoint ){
      CharString cigarOperation = record.cigar[length(record.cigar)-1].operation;
      string cigarOperationStr = toCString(cigarOperation);
      if( cigarOperationStr == "S" and ( (int) record.cigar[length(record.cigar)-1].count > O.maxSoftClipped) ){
	softClipRemove = true;
	if( O.verbose ) cerr << "SoftClip removed RightBreakpoint " << record.qName << " " << O.genotypeRightBreakpoint << " " << record.cigar[length(record.cigar)-1].count << " " << O.maxSoftClipped << endl;
      }
    }
    bool hardClipped = false;
    CharString cigarOperationL = record.cigar[0].operation;
    string cigarOperationStrL = toCString(cigarOperationL);
    CharString cigarOperationR = record.cigar[length(record.cigar)-1].operation;
    string cigarOperationStrR = toCString(cigarOperationR);
    if( (cigarOperationStrL == "H") or (cigarOperationStrL == "H") ){
      hardClipped = true;
      if( O.verbose ) cout << "Read " << record.qName << " is hardclipped at " << record.beginPos << endl;
    }
    if( (not softClipRemove) and (not hasFlagDuplicate( record )) and (not hasFlagQCNoPass( record )) and (not hardClipped)){
      bars[record.qName] = record;
      vaiMap[record.qName] = vai;
    }
    if( O.verbose ) cerr << "Finished soft clipping " << endl;
  }
  if( O.verbose ) cerr << "Exiting readBamRegion " << endl;
  return 0;
}

// Input: variant and VarAlignInfo records for each read overlapping variant
// Output: Relative genotype likelihoods in log_2 scale and read info counts
template< typename TOptions>
int multiUpdateVC( VcfRecord& var, map< CharString, varAlignInfo>& vais1, map< CharString, varAlignInfo>& vais2, 
		   vector< double>& vC, vector< int>& rI, TOptions& O, genotyping_model gtm  ){

  // Loop over all reads in bam file(s) 1 (deletion biased calls)
  for( auto i = vais1.begin(); i != vais1.end(); i++ ){
    if( (*i).second.alignmentAltPreference( O ) < 0 ) rI[0]++;
    if( (*i).second.alignmentAltPreference( O ) > 0 ) rI[1]++;
    rI[2]++;
    int update = 0;
    if( vais2.count( (*i).first ) == 0 ){
      //read does not occur in bam file(s) 2 (insertion biased calls)
      if( gtm == genotype_ad or gtm == genotype_joint ) update += (*i).second.alignmentAltPreference( O );
      if( gtm == genotype_va or gtm == genotype_joint ) update += O.overlapBits*((*i).second.supports( var, O ) - (*i).second.rejects( var, O ));
      if( O.verbose ) cerr << "Update VC bam1 " << (*i).first << " " << (*i).second.alignmentAltPreference( O ) << " " << (*i).second.supports( var, O ) << " " << (*i).second.rejects( var, O ) << " " << update << endl;
    }else{
      // read is also in the second file, current logic handles these independently
      varAlignInfo vai2 = vais2[(*i).first];
      if( gtm == genotype_ad or gtm == genotype_joint ) update += (*i).second.alignmentAltPreference( O ) + vai2.alignmentAltPreference( O); 
      if( gtm == genotype_va or gtm == genotype_joint ) update += O.overlapBits*((*i).second.supports( var, O ) - (*i).second.rejects( var, O ) + vai2.supports( var, O ) - vai2.rejects( var, O ));
      if( O.verbose ) cerr << "Update VC two bams " << (*i).first << " " << (*i).second.alignmentAltPreference( O ) << " " << vai2.alignmentAltPreference( O) << " " << (*i).second.supports( var, O ) << " " << (*i).second.rejects( var, O ) << " " << vai2.supports( var, O ) << " " << vai2.rejects( var, O ) << " " << update << endl;
    }      
    if( update > 1 ){
      vC[0] += update;
      vC[1] += 1;
    }else if( update < -1 ){
      vC[2] += (-update);
      vC[1] += 1;
    }
  }
  //Find reads in the second bam(s) not in the first
  for( auto i = vais2.begin(); i != vais2.end(); i++ ){
    if( (*i).second.alignmentAltPreference( O ) < 0 ) rI[0]++;
    if( (*i).second.alignmentAltPreference( O ) > 0 ) rI[1]++;
    rI[2]++;
    if( vais1.count( (*i).first ) == 0 ){
      int update = 0;
      if( gtm == genotype_ad or gtm == genotype_joint ) update += (*i).second.alignmentAltPreference( O );
      if( gtm == genotype_va or gtm == genotype_joint ) update += O.overlapBits*((*i).second.supports( var, O ) - (*i).second.rejects( var, O ));
      if( O.verbose ) cerr << "Update VC bam2 " << (*i).first << " " << (*i).second.alignmentAltPreference( O ) << " " << (*i).second.supports( var, O ) << " " << (*i).second.rejects( var, O ) << " " << update << endl;
      if( update > 1 ){
	vC[0] += update;
	vC[1] += 1;
      }else if( update < -1 ){
	vC[2] += -update;
	vC[1] += 1;
      }
    }
  }
  
  if( O.verbose ) cerr << "multiUpdateVC " << vC[0] << " " << vC[1] << " " << vC[2] << endl;
  return 0;
}


template<typename TOptions>
int computeVAs( VcfRecord& variant, map< CharString, varAlignInfo>& vais, vector< int>& VA, TOptions& O ){
  for( auto i = vais.begin(); i != vais.end(); i++ ){
    if( (*i).second.rejects( variant, O ))
      VA[0]++;
    if( (*i).second.supports( variant, O ))
      VA[1]++;
    if( (*i).second.aligns())
      VA[2]++;
  }
  return 0;
}


/** Input: bamStream, a VCF entry, reference fasta and alt fasta file if required by VCF entry
    Output: variant alignment info for each read near the VCF entry
 */
template<typename TOptions>
int LRprocessReads(VcfRecord & variant, VcfFileIn & vcfS,FaiIndex & faiI,
		   map< CharString, BamAlignmentRecord>& bars,
		   TOptions & O, map< CharString, varAlignInfo>& vais )
{
  TSequence refSeq;
  TSequence altSeq;
  CharString chrom;
  int32_t beginPos;

  // move outside of this function
  getLocRefAlt( variant, vcfS, faiI, chrom, beginPos, refSeq, altSeq, O );

  if( O.outputRefAlt ){
    cout << chrom << " " << beginPos + 1 << " " << variant.info << " " << refSeq << " " << altSeq << endl;
    return 0;
  }

  int nReads = (int) bars.size();
  if( O.verbose ) cerr << "LRprocessReads " << nReads << " " << vais.size() << endl;

  // communication is via files if we are using minimap or BWA for alignment
  if( O.aligner == aligner_minimap or O.aligner == aligner_bwa ){
    //minimap logic 
    writeFastaFromBars( bars, O.maxBARcount, variant, O );
    ofstream minimapF;
    stringstream ss2;
    ss2 << "variants_to_minimap_" << O.randID << ".fa";
    string str2 = ss2.str();
    minimapF.open( str2 );
    minimapF << ">Ref" << endl;
    minimapF << refSeq << endl;
    minimapF << ">Alt" << endl;
    minimapF << altSeq << endl;
    minimapF.close();
    
    std::stringstream cmd;
    cmd.str("");
    if( O.aligner == aligner_minimap ){
      //cmd << "./minimap2 -c reads_to_minimap.fa variants_to_minimap > minimap_out.paf";   
      // PAF might be more efficient more efficient, but cost is probably mainly related to creating cigar
      cmd << "minimap2 -a -A " << O.match << " -B " << -1*O.mismatch << " -O " << -1*O.gapOpen << " -E " << -1*O.gapExtend << " variants_to_minimap_" << O.randID << ".fa reads_to_minimap_" << O.randID << ".fa > mapping_out_" << O.randID << ".sam";   
    }else if( O.aligner == aligner_bwa){
      cmd << "bwa index variants_to_minimap_" << O.randID << ".fa";
      if (system(cmd.str().c_str()) != 0){
	std::cerr << "ERROR in bwa index " << std::endl;
	return 7;
      } 
      cmd.str("");
      cmd << "bwa mem -A " << O.match << " -B " << -1*O.mismatch << " -O " << -1*O.gapOpen << " -E " << -1*O.gapExtend << " variants_to_minimap_" << O.randID << ".fa reads_to_minimap_" << O.randID << ".fa > mapping_out_" << O.randID << ".sam";   
    }else{
      cout << "Unknown aligner " << endl;
      exit(-1);
    }
    if (system(cmd.str().c_str()) != 0){
      std::cerr << "ERROR in bwa " << std::endl;
      return 7;
    } 

    stringstream ss;
    ss << "mapping_out_" << O.randID << ".sam";
    string str = ss.str();

    BamFileIn bamFileIn(str.c_str());

    // Get rid of header
    BamHeader header;
    readHeader(header, bamFileIn);

    int refID = 0;
    if (!getIdByName( refID, contigNamesCache(context( bamFileIn )), "Ref" ) and O.verbose)
      cout << "getIdByName failed, could not find Ref" << endl;
    int altID = 0;
    if (!getIdByName( altID, contigNamesCache(context( bamFileIn )), "Alt" ) and O.verbose)
      cout << "getIdByName failed, could not find Alt" << endl;
    if( O.verbose ) cerr << "Got ids by name Ref " << refID << " Alt " << altID << endl;
    //    for( int i = 0; i < nReads; i++ )
    //refS[i] = altS[i] = NO_ALIGNMENT;
    BamAlignmentRecord record;
    while (!atEnd(bamFileIn))
    {
        readRecord(record, bamFileIn);
	BamTagsDict tagsDict(record.tags);
	unsigned tagIdx = 0;
	if (!findTagKey(tagIdx, tagsDict, "AS")){
	  if( O.verbose ) cerr << "Unknown key AS in bam tag\n";
	  //Then, we can read the value from the BamTagsDict using the function extractTagValue.
	}else{
	  int tagValInt = 0;
	  if (!extractTagValue(tagValInt, tagsDict, tagIdx))
	    std::cerr << "ERROR: There was an error extracting AS from tags!\n";
	  if( O.verbose ) cerr << record.qName << " " << record.rID << " " << tagValInt << endl;
	  if( (record.rID == refID) && (vais[record.qName].refS == NO_ALIGNMENT || vais[record.qName].refS < tagValInt) )
	    vais[record.qName].refS = (double) tagValInt;
	  else if( (record.rID == altID) && (vais[record.qName].altS == NO_ALIGNMENT || vais[record.qName].altS < tagValInt) )
	    vais[record.qName].altS = (double) tagValInt;
	}
    }
  }else if( O.aligner == aligner_seqan ){
    TAlign align;
    Score<int, Simple> scoringScheme(O.match, O.mismatch, O.gapExtend, O.gapOpen);

    for( auto i = bars.begin(); i != bars.end(); i++ ){
      CharString seqToAlign;
      if( O.cropRead ){
	cropSeq((*i).second, variant, O, seqToAlign);
      }else{
	seqToAlign = (*i).second.seq;
      }
      resize( rows( align ), 2 );
      assignSource(row(align,0), refSeq);
      assignSource(row(align,1), seqToAlign );
      vais[(*i).second.qName].refS = (int) globalAlignment(align, scoringScheme, AlignConfig< true, true, true, true>());
      assignSource(row(align,0), altSeq);
      vais[(*i).second.qName].altS = (int) globalAlignment(align, scoringScheme, AlignConfig< true, true, true, true>());
      if( O.verbose ) cerr << "Aligned " << (*i).second.qName << " " << vais[(*i).second.qName].refS  << " " << vais[(*i).second.qName].altS << endl;
    }
  }
  return 0;
}


// Open a bam file or a set of bam files if the filename does not end with .bam
int parseBamFileName( CharString& bfN, vector< BamIndex<Bai>>& bamIndexV, vector<std::unique_ptr<BamFileIn> >& bamStreamV )
{
  std::string bfNs( toCString(bfN));
  int bamStrLoc = bfNs.find(".bam" );  // ACHTUNG: This is not the best method, checks whether the first occurence of the string ".bam" in the filename is in last for letters of the name, code crashes for file names such as file.bamfile.bam
  if( bamStrLoc == (int) (length( bfN ) - 4) ){
    BamIndex<Bai> bamIndex;
    std::unique_ptr<BamFileIn> bamStream = std::unique_ptr<BamFileIn>(new BamFileIn);
    if (initializeBam(toCString(bfN), bamIndex, *(bamStream.get())) != 0)
      return 7;
    bamIndexV.push_back( bamIndex );
    bamStreamV.push_back( std::move(bamStream) );
  }else{
    ifstream fS;
    fS.open( toCString(bfN) );
    string bf;
    while( fS >> bf ){
      BamIndex<Bai> bamIndex;
      std::unique_ptr<BamFileIn> bamStream = std::unique_ptr<BamFileIn>(new BamFileIn);
      char* bfc = (char*) bf.c_str(); 
      if (initializeBam(bfc, bamIndex, (*bamStream)) != 0)
	return 7;
      bamIndexV.push_back( bamIndex );
      bamStreamV.push_back( std::move(bamStream) );
    }
  }
  return 0;
} 

int main(int argc, char const ** argv)
{

  std::random_device rd;
  std::default_random_engine e( rd() );
  // Parse the command line to get option values.
  LRCOptions O;
  O.randID = e(); 
  if( parseLRCArguments(argc, argv, O ) != ArgumentParser::PARSE_OK ) return 1;

  int nRecords = 0;  
  VcfFileIn vcfIn(toCString(O.vcfInFile));
  VcfHeader header;
  readHeader(header, vcfIn);
  vector< VcfRecord> records;
  records.resize(100);
  while (!atEnd(vcfIn)){
    readRecord(records[nRecords], vcfIn);
    nRecords++;
    if( nRecords >= (int) records.size() ) records.resize(2*records.size());
  }
  records.resize(nRecords);
  if( O.verbose ) std::cerr << "nRecords " << nRecords << " " << records.size() << std::endl;
  // Open the input VCF file and prepare output VCF stream.
  //  VcfFileIn vcfIn(toCString(O.vcfInFile));
  appendValue(header, VcfHeaderRecord("FORMAT", "<ID=GT,Number=1,Type=String,Description=\"Genotype\">"));
  appendValue(header, VcfHeaderRecord("FORMAT", "<ID=AD,Number=3,Type=Integer,Description=\"Allelic depths from alignment supporting ref and alt allele and total number of reads\">"));
  appendValue(header, VcfHeaderRecord("FORMAT", "<ID=VA,Number=3,Type=Integer,Description=\"Allelic depths from bam file supporting ref and alt allele and total number of reads\">"));
  appendValue(header, VcfHeaderRecord("FORMAT", "<ID=PL,Number=G,Type=Integer,Description=\"PHRED-scaled genotype likelihoods\">"));
  std::ofstream vcfStream(toCString(O.vcfOutFile));
  VcfFileOut vcfOut(vcfIn);
  open(vcfOut, vcfStream, Vcf());
  writeHeader(vcfOut, header);

  // Build an index of the fasta file (reference genome).
  FaiIndex faIndex;

  // Open the bam file(s). (A bam file needs the bai index and the bam file stream.)
  vector< BamIndex<Bai>> bamIndexV;
  vector< BamFileIn*> bamStreamV;
  //  parseBamFileName( O.bam, bamIndexV, bamStreamV );

  bool useBam2 = false;
  vector< BamIndex<Bai>> bamIndexV2;
  vector< BamFileIn*> bamStreamV2;
  if( O.bam2 != "" ){
    useBam2 = true;
    if( O.verbose ) cerr << "Using bam2 " << O.bam2 << endl;
    //    parseBamFileName( O.bam2, bamIndexV2, bamStreamV2 );
  }

  // Build an index of the insertion sequences' fasta file, if such a file is used.

  std::cerr << "Genotyping VCF records in \'" << O.vcfInFile << std::endl;

  if( O.verbose ) std::cerr << "setting number of threads " << O.nThreads << std::endl;
 
  int nThreads = omp_get_num_threads();
  if( O.verbose ) std::cerr << "number of threads set " << nThreads << std::endl;

  int num_threads{O.nThreads};     
  omp_set_num_threads(num_threads);     
  vector< vector< BamIndex<Bai>>>  bamIndexHandles;
  vector< vector< BamIndex<Bai>>>  bamIndex2Handles;
  std::vector< std::vector< std::unique_ptr<BamFileIn> > > bamHandles;     
  std::vector< std::vector< std::unique_ptr<BamFileIn> > > bam2Handles;     
  std::vector< FaiIndex> faIndexHandles;     
  bamHandles.resize(num_threads);     
  bam2Handles.resize(num_threads);     
  bamIndexHandles.resize(num_threads);     
  bamIndex2Handles.resize(num_threads);     
  faIndexHandles.resize(num_threads);
  if(O.verbose ) std::cerr << "created vectors " << std::endl;     
  for (int t = 0; t < num_threads; t++){
    if(O.verbose ) std::cerr << "Creating thread " << t << std::endl;     
    parseBamFileName( O.bam, bamIndexHandles[t], bamHandles[t] );
    if( useBam2 ) parseBamFileName( O.bam, bamIndex2Handles[t], bam2Handles[t] );
    if (!open(faIndexHandles[t], toCString(O.faFile))){
      if (!build(faIndexHandles[t], toCString(O.faFile))){
	std::cerr << "ERROR: Could not find nor build the index of " << O.faFile << std::endl;
	return 7;
      }
    }
  }

  #pragma omp parallel for
  for( int nR = 0; nR < nRecords; nR++ ){ 
    if( O.verbose ) std::cerr << "Genotyping " << nR << " thread " << omp_get_thread_num() << std::endl;
    map< CharString, varAlignInfo> vaisR;
    map< CharString, varAlignInfo> vais2R;

    CharString varChrom;  
    varChrom = contigNames(context(vcfIn))[records[nR].rID];

    // Implemented as a vector as we might want to have more than one output per marker
    std::vector< std::vector<double> > vC; //variant calls
    std::vector< std::vector<int> > AD;  // Allele depth counts
    std::vector< std::vector<int> > VA;  // variant alignment counts
    vC.resize( 1);
    AD.resize(1);
    VA.resize(1);
    if( O.gtModel == genotype_multi ){
      vC.resize( 3 );
      AD.resize( 3 );
      VA.resize( 3 );
    }
    for( int mI = 0; mI < (int) vC.size(); mI++ ){
      vC[mI].resize(3);
      AD[mI].resize(3);
      VA[mI].resize(3);
    }
    for( int i = 0; i < (int) bamIndexHandles[omp_get_thread_num()].size(); i++ ){
      std::map< CharString, BamAlignmentRecord> bars;
      readBamRegion( bamIndexHandles[omp_get_thread_num()][i], (*bamHandles[omp_get_thread_num()][i]), varChrom, records[nR], O, bars, vaisR );
      if( O.verbose ) cerr << "Finished reading bam region bam1 ";
      LRprocessReads(records[nR], vcfIn, faIndexHandles[omp_get_thread_num()], bars, O, vaisR );
    }
    if( useBam2 ){
      for( int i = 0; i < (int) bamIndex2Handles[omp_get_thread_num()].size(); i++ ){
	std::map< CharString, BamAlignmentRecord> bars;
	readBamRegion( bamIndex2Handles[omp_get_thread_num()][i], (*bam2Handles[omp_get_thread_num()][i]), varChrom, records[nR], O, bars, vaisR );
	LRprocessReads(records[nR], vcfIn, faIndexHandles[omp_get_thread_num()], bars, O, vais2R );
      }
    }
    if( O.gtModel == genotype_multi ){
      multiUpdateVC( records[nR], vaisR, vais2R, vC[0], AD[0], O, genotype_ad );
      multiUpdateVC( records[nR], vaisR, vais2R, vC[1], AD[1], O, genotype_va );
      multiUpdateVC( records[nR], vaisR, vais2R, vC[2], AD[2], O, genotype_joint );
      computeVAs( records[nR], vaisR, VA[1], O);
      if( useBam2 ) computeVAs( records[nR], vais2R, VA[1], O);
      computeVAs( records[nR], vaisR, VA[2], O);
      if( useBam2 ) computeVAs( records[nR], vais2R, VA[2], O);
    }else{
      multiUpdateVC( records[nR], vaisR, vais2R, vC[0], AD[0], O, O.gtModel );
      computeVAs( records[nR], vaisR, VA[0], O);
      if( useBam2 ) computeVAs( records[nR], vais2R, VA[0], O);
    }
    std::string gtString;
      
    for( int mI = 0; mI < (int) vC.size(); mI++ ){
      //      VcfRecord thisR = record;
      getGtString(vC[mI], AD[mI], VA[mI], gtString );
      appendValue(records[nR].genotypeInfos, gtString);
      records[nR].format = "GT:AD:VA:PL";
      //writeRecord(vcfOut, thisR);
    }
  }
  for( int nR = 0; nR < nRecords; nR++ ){
    writeRecord( vcfOut, records[nR] );
  }
  return 0;
}


