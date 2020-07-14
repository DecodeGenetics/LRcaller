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
#define NO_BEST -1

enum aligner_type{ aligner_minimap, aligner_bwa, aligner_seqan };

enum genotyping_model{ genotype_multi, genotype_ad, genotype_va, genotype_presence, genotype_va_old, genotype_joint};


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
  vector< int> randID;  // To avoid file name collision when runing bwa/minimap
  bool outputRefAlt;  // For debug purposes
  int maxSoftClipped; // Read cannot be softclipped at invariant end
  double maxAlignBits;  //needs to be added as option
  int overlapBits;  // needs to be added as option
  int varWindow;  //How far 
  int maxReadLenForCropping;  // Should not be needed
  int minMapQ; // minimum mapQ for a read to be used
  int minDelIns; // minum number of consecutive bp to be deleted or inserted to count it
  int minPresent; // minum number of del/ins for a variant to be considered present
  double refThreshFraction;  // Fraction of deletion/insertion length that can be deleted/inserted but still called ref
  double altThreshFraction;  // Fraction of deletion/insertion length that needs to be deleted/inserted to call alt
  double altThreshFractionMax;  // Fraction of deletion/insertion length can be at most be observed to be deleted/inserted to call alt
  genotyping_model gtModel;
  int nThreads;
  LRCOptions():  bam2(""),faFile("genome.fa"),wSize(500), maxBARcount(200),verbose(false), genotypeRightBreakpoint(false), minAlignScore(550),logScaleFactor(2.0),genotypingModel("multi"),match(1),mismatch(-1),gapOpen(-1),gapExtend(-1),cropRead(true),aligner(aligner_seqan),outputRefAlt(false), maxSoftClipped(500), maxAlignBits( 10.0), overlapBits(10), varWindow(100), maxReadLenForCropping(200000), minMapQ(30), minDelIns(6), minPresent(30), refThreshFraction(1.0/10.0), altThreshFraction(1.0/2.0), altThreshFractionMax(100.0), gtModel(genotype_multi), nThreads(1){};
  //regionWindowSize - how far we want to look in the bam file
};

template<typename TOptions>
int parseLRCArguments( int argc, char const ** argv, TOptions & O ){

  ArgumentParser parser("LRcaller");
  addUsageLine(parser, "[\\fIOPTIONS\\fP]  \"\\fIBAMFILE\\fP\"  \"\\fIVCF_FILE_IN\\fP\" \"\\fIVCF_FILE_OUT\\fP\" ");
  addDescription(parser, "Genotypes variants using long reads\n in bam file/file of bam files ");
 
  addArgument(parser, ArgParseArgument(ArgParseArgument::STRING, "BAMFILE bam file/file of bam files"));
  addArgument(parser, ArgParseArgument(ArgParseArgument::STRING, "VCF_FILE_IN - input vcf file"));
  addArgument(parser, ArgParseArgument(ArgParseArgument::STRING, "VCF_FILE_OUT - output genotyped vcf file"));
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

  addOption(parser, ArgParseOption("ma", "min_align_score", "Minimum alignment score", ArgParseArgument::INTEGER, "INT"));
  setDefaultValue(parser, "ma", O.minAlignScore);    

  addOption(parser, ArgParseOption("mapQ", "read_mapq", "Minimum read mapQ", ArgParseArgument::INTEGER, "INT"));
  setDefaultValue(parser, "mapQ", O.minMapQ);    

  addOption(parser, ArgParseOption("di", "min_del_ins", "Minimum del/ins in cigar to consider", ArgParseArgument::INTEGER, "INT"));
  setDefaultValue(parser, "di", O.minDelIns);    

  addOption(parser, ArgParseOption("mp", "min_present", "Minimum del/ins var a variant to be considered present", ArgParseArgument::INTEGER, "INT"));
  setDefaultValue(parser, "mp", O.minPresent);    

  addOption(parser, ArgParseOption("rt", "ref_thresh_fraction", "Threshold for fraction of del/ins bp to be considered ref", ArgParseArgument::DOUBLE, "DOUBLE"));
  setDefaultValue(parser, "rt", O.refThreshFraction);    

  addOption(parser, ArgParseOption("at", "alt_thresh_fraction", "Threshold for fraction of del/ins bp to be considered alt", ArgParseArgument::DOUBLE, "DOUBLE"));
  setDefaultValue(parser, "at", O.altThreshFraction);    

  addOption(parser, ArgParseOption("am", "alt_thresh_fraction_max", "Max threshold for fraction of del/ins bp to be considered alt", ArgParseArgument::DOUBLE, "DOUBLE"));
  setDefaultValue(parser, "am", O.altThreshFractionMax);    

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

  if (isSet(parser, "min_align_score"))
    getOptionValue(O.minAlignScore, parser, "min_align_score");

  if (isSet(parser, "read_mapq"))
    getOptionValue(O.minMapQ, parser, "read_mapq");

  if (isSet(parser, "min_del_ins"))
    getOptionValue(O.minDelIns, parser, "min_del_ins");

  if (isSet(parser, "min_present"))
    getOptionValue(O.minPresent, parser, "min_present");

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
  else if( GTmodelName == "va2" )
    O.gtModel = genotype_va_old;
  else if( GTmodelName == "p" )
    O.gtModel = genotype_presence;
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
  int nAlleles;
  vector< double> alignS;
  bool softClipped;
  bool alignsLeft;
  bool alignsRight;
  //int refAlignLen;

  // Supports and rejects logic is not correct for very long variants
  // Alternate allele is supported as judged from bam alignment record
  template<typename TOptions>
  bool supports(float refLen, float altLen, TOptions& O ){
    if( altLen > refLen){  // insertion (these are simplistic for insertion/deletion type variants)
      // doesn't work properly if alt and ref are of similar size

      if( (alignsLeft and alignsRight and ((float) nI > (altLen)*O.altThreshFraction) and ((float) nI < altLen*O.altThreshFractionMax)) or softClipped )
	return true;
      else
	return false;
    }else{
      if( (alignsLeft and alignsRight and ((float) nD > refLen*O.altThreshFraction) and ((float) nD < refLen*O.altThreshFractionMax)) or softClipped)
	return true;
      else
	return false;
    }
  }

  // Alternate allele is rejected as judged from bam alignment record
  template<typename TOptions>
  bool rejects( float refLen, float altLen, TOptions& O  ){
    if( altLen > refLen ){  // insertion
      if( (alignsLeft and alignsRight and ((float) nI < altLen*O.refThreshFraction)) and (not softClipped) )
	return true;
      else
	return false;
    }else{
      if( (alignsLeft and alignsRight and ((float) nD < refLen*O.refThreshFraction)) and (not softClipped) )
	return true;
      else
	return false;
    }
  }

  template<typename TOptions>
  bool present( TOptions& O ){
    if( nI >= O.minPresent or nD >= O.minPresent )
      return true;
    else 
      return false;
  }

  bool aligns(){
    return alignsLeft and alignsRight;
  }

  void reset(){
    nD = 0;
    nI = 0;
    for( int i = 0; i < (int) alignS.size(); i++ )
      alignS[i] = NO_ALIGNMENT;
    softClipped = false;
    alignsLeft = false;
    alignsRight = false;
  }

  // Likelihood of variant relative to the most likely, 
  // Updates pref values. A value x, represents that the allele is 2^-x less likely than the most likely allele 
  // Returns an index to the most likely allele, if one exists
  template<typename TOptions>
    int alignmentPreference( TOptions& O, vector< double>& pref ){
    int maxScore = alignS[0];
    int maxI = 0;
    for( int i = 0; i < nAlleles; i++ ){
      // pref[i] = 0;
      if( alignS[i] > maxScore ){
	maxI = i;
	maxScore = alignS[i];
      }
    }
    if( maxScore == NO_ALIGNMENT or maxScore <= O.minAlignScore ){
      return NO_BEST;
    }else{
      for( int i = 0; i < nAlleles; i++ ){
	float d = (maxScore - alignS[i])/O.logScaleFactor;
	if( alignS[i] == NO_ALIGNMENT or alignS[i] <= O.minAlignScore ){
	  d = (maxScore - O.minAlignScore)/O.logScaleFactor;
	}
	if( d > O.maxAlignBits ) d = O.maxAlignBits;
	if( d < 0 and O.verbose  ) cerr << "WTF negative alignment score";
	pref[i] += d;
      }
      return maxI;
    }
  }


  // Likelihood of variant relative to the most likely, 
  // Increments pref values. A value x, represents that the allele is 2^-x less likely than the most likely allele 
  // Returns an index to the most likely allele, if one exists
  template<typename TOptions>
  int vaPreference( TOptions& O, int refLen, vector< int> & altLens, vector< double>& pref ){

    if( softClipped ){  //Softclipped, doesn't support the reference, all other alleles are equally likely
      pref[0] += O.overlapBits;
      return NO_BEST;
    }

    if( not alignsLeft or not alignsRight ){
      return NO_BEST;
    }


    
    // Count number of deletions and insertions, find the variant that is closest in size

    int insDel = nI - nD;
    int minD = abs(insDel);
    int minDi = 0;

    for( int i = 1; i < nAlleles; i++ ){
      int cD = altLens[i-1] - refLen;
      if( abs( cD - insDel ) < minD ){
	minDi = i;
	minD = abs(cD-insDel);
      }
    }
    for( int i = 0; i < nAlleles; i++ ){
      if( i != minDi ) pref[i] += O.overlapBits;
    }
    return minDi;
  }


  //  varAlignInfo(): nD(0),nI(0), refS(NO_ALIGNMENT),altS(NO_ALIGNMENT),softClipped(false),alignsLeft(false),alignsRight(false){};
  varAlignInfo( int nAllelesIn ){
    nAlleles = nAllelesIn;
    nD = 0; nI = 0; 
    alignS.resize( nAlleles );
    for( int i = 0; i < nAlleles; i++ ){
      alignS[i] = NO_ALIGNMENT; 
    }
    softClipped = false; 
    alignsLeft = false; 
    alignsRight = false;
  }
  varAlignInfo(){
    varAlignInfo(2);
  }
};


/* Turns genotyping into string */
void
getGtString(std::vector<double> & lls, std::vector< int> & ads, std::vector< int>& vas, std::string & gtString )
{
  int gtLen = (int) lls.size();
  for( int i = 0; i < gtLen; i++ )
    lls[i] = -lls[i];   // Really silly hack for historical reasons, would confuse the hell out of me to fix it

  //  int max = 0;
  double maxP = lls[0];
  int a1 = 0;
  int a2 = 0;
  int maxA1 = 0;
  int maxA2 = 0;
  for( int i = 0; i < gtLen; i++ ){
    if( lls[i] > maxP ){
      maxP = lls[i];
      //      max = i;
      maxA1=a1;
      maxA2=a2;
    }
    if( a2 < a1 )
      a2++;
    else{
      a1++;
      a2=0;
    }
  }
  std::ostringstream buff;
  buff << maxA2 << "/" << maxA1 << ":";
  for( int i = 0; i < (int) ads.size()-1; i++ )
    buff << ads[i] << ",";
  buff << ads[ads.size()-1] << ":";

  for( int i = 0; i < (int) vas.size()-1; i++ )
    buff << vas[i] << ",";
  buff << vas[vas.size()-1] << ":";
  
  for( int i = 0; i < gtLen; i++ ){
    double lp = (lls[i]-maxP)/LG10;
    if( lp < LL_THRESHOLD ) lp = LL_THRESHOLD;
    buff << int( -10*lp );
    if( i != gtLen -1 ) buff << ",";
  }
  gtString = buff.str();
}


//Input: VcfRecord, reference and alternate reference
//Output: Chrom and position of variant, ref sequence and alt sequence
//TODO: Change this for a library that does this
template<typename TOptions>
int getLocRefAlt( VcfRecord & variant, VcfFileIn & vcfS,FaiIndex & faiI, 
		  CharString& chrom, int32_t& beginPos, TSequence& refSeq, vector<TSequence>& altSeqs, TOptions & O )
{
    TSequence ref = variant.ref;
    StringSet< CharString> altSet;
    strSplit(altSet,variant.alt,EqualsChar<','>());

    int nAlts = (int) length( altSet );
    //    CharString alt = variant.alt;

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
    
    vector< int> altLens;
    altLens.resize( nAlts );
    for( int i = 0; i < nAlts; i++ ){
      altLens[i] = (int) length( altSet[i] );
    }
    int refLen = length( variant.ref );
    for( int i = 0; i < nAlts; i++ ){
      if( not O.genotypeRightBreakpoint ){
	readRegion(altSeqs[i], faiI, idx, beginPos-O.wSize, beginPos ); // beginPos is included in the alt sequence
	TSequence post;
	if( altLens[i] < O.wSize ){
	  post = altSet[i];
	  TSequence post2;
	  readRegion(post2, faiI, idx, beginPos+refLen, beginPos+refLen + O.wSize - altLens[i] );
	  append( post, post2);
	}else{
	  post = infixWithLength( altSet[i], 0, O.wSize );
	}
	append( altSeqs[i], post );
      }else{	
	if( altLens[i] < O.wSize ){
	  readRegion(altSeqs[i], faiI, idx, beginPos-O.wSize+altLens[i], beginPos );
	  append( altSeqs[i], altSet[i] );
	}else{
	  altSeqs[i] = infixWithLength( altSet[i], altLens[i]-O.wSize, O.wSize );
	}
	TSequence post;
	readRegion(post, faiI, idx, beginPos+refLen, beginPos + refLen+O.wSize );
	append( altSeqs[i], post );
      }
    }
  
    if( O.verbose ){
      cerr << "Printing altSeq " << endl;
      for( int i = 0; i < nAlts; i++ )
	cerr << "altSeq " << i << " " << altSeqs[i] << endl;
      cerr << "Done printing altSeq " << endl;
    }
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
  int searchPos = var.beginPos-O.wSize;  // Want to change the search intervals for TRs
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
    if( O.verbose ) cerr << bar.qName << " readpos " << alignPos << " " << var.beginPos << " " << searchPos << " " << cigarI << " " << readPos << " " << cigarString[cigarI].count << " " << cigarOperationStr << endl;  
    cigarI++;
  }
  if( alignPos < searchPos and O.verbose ){
    cout  << "Read clipped " << alignPos << " " << var.beginPos << " " << searchPos << " " << cigarI << " " << length(cigarString) << endl;  
  }  

  if( cigarOperationStr == "S" or cigarOperationStr == "H") readPos = lReadPos;
  
  int rBeg,rEnd;
  if( O.genotypeRightBreakpoint ){
    //    if( alignPos >= searchPos ){
    //   rBeg = readPos - 2*O.wSize;
    // rEnd = readPos;
    //}else 
    if( alignPos >= searchPos - 2*O.wSize ){
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
    if( rShift < 0 and O.verbose ) cerr << "Poorly formatted read, case not accounted for " << bar.qName << " " << alignPos << " " << var.beginPos << " " << searchPos << " " << cigarI << " " << length(cigarString) << endl;
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
  int regionBeg = var.beginPos - O.varWindow;
  int regionEnd = var.beginPos+(int) length(var.ref)+O.varWindow;
  StringSet< CharString> infos;
  strSplit(infos,var.info,EqualsChar<';'>());
  for( int i = 0; i < (int) length( infos ); i++ ){
    StringSet< CharString> info2;
    strSplit( info2, infos[i],EqualsChar<'='>());
    int cVal;
    if( info2[0] == "TRRBEGIN" ){
      if( info2[1] != "." ){
	cVal = atoi(toCString(info2[1]))-O.varWindow;
	if( cVal < regionBeg )
	  regionBeg = cVal;
      }
    }
    if( info2[0] == "TRREND" ){
      if( info2[1] != "." ){
	cVal = atoi(toCString(info2[1]))+O.varWindow;
	if( cVal > regionEnd )
	  regionEnd = cVal;
      }
    }
    if( info2[0] == "REGBEGIN" ){
      if( info2[1] != "." ){
	cVal = atoi(toCString(info2[1]))-O.varWindow;
	if( cVal < regionBeg )
	  regionBeg = cVal;
      }
    }
    if( info2[0] == "REGEND" ){
      if( info2[1] != "." ){
	cVal = atoi(toCString(info2[1]))+O.varWindow;
	if( cVal > regionEnd )
	  regionEnd = cVal;
      }
    }
  }
  /*  if( infos.count( "TRRBEGIN" ) == 1 and infos.count( "TRREND" ) == 1 ){
    regionBeg = atoi(var.infos["TRRBEGIN"] ) - O.varWindow;
    CharString endStr = getValueById( var.genotypeInfos, trrEndID);
    regionEnd = atoi( var.infos["TRREND"] ) + O.varWindow;
    if( O.verbose ) cerr << "Found TRR " << regionBeg << " " << regionEnd << endl;
}else{*/
    if( O.verbose ) cerr << "TRR " << regionBeg << " " << regionEnd << endl;
    //  }
  

  if( alignPos < regionBeg ) vai.alignsLeft = true;

  // Find the first position that overlaps the window we are interested in
  while( alignPos < regionBeg and cigarI < (int) length( cigarString) ){
    cigarOperation = cigarString[cigarI].operation;
    cigarOperationStr = toCString(cigarOperation);
    if( cigarOperationStr == "M" or cigarOperationStr == "D" ){
      alignPos += cigarString[cigarI].count;
    }
    cigarI++;
  }

  // only counts the number of deleted bp in the window, probably fine to count a longer distance
  if( alignPos > regionBeg and cigarOperationStr == "D" and (alignPos - (regionBeg) >= (int) O.minDelIns) ){
    vai.nD = alignPos - (regionBeg);
  }

  while( alignPos < regionEnd and cigarI < (int) length( cigarString) ){
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
  if( alignPos > regionEnd ) vai.alignsRight = true;
  if( O.verbose ) cerr << "examinSeq " << bar.qName << " " << vai.nD << " " << vai.nI << " " << vai.softClipped << endl;
  return 0;
}

// writes fasta file for BWA or minimap from bam alignment records
template<typename TOptions>
int writeFastaFromBars(   std::map< CharString, BamAlignmentRecord>& bars, int& maxBC, VcfRecord& var, TOptions& O )
{
  ofstream minimapF;
  stringstream ss;
  ss << "reads_to_minimap_" << O.randID[omp_get_thread_num()] << ".fa";
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

  StringSet< CharString> altSet;
  strSplit(altSet,var.alt,EqualsChar<','>());

  int nAlts = (int) length( altSet );


  BamAlignmentRecord record;
  varAlignInfo vai( nAlts + 1);

  
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
      if( O.verbose ) cerr << "Read " << record.qName << " is hardclipped at " << record.beginPos << endl;
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
  int multiUpdateVC( VcfRecord& var, map< CharString, varAlignInfo>& vais, 
		     vector< double>& vC, vector< int>& rI, vector< int>& VAs, TOptions& O, genotyping_model gtm  ){

  StringSet< CharString> altSet;
  strSplit(altSet,var.alt,EqualsChar<','>());

  int nAlts = (int) length( altSet );

  vector< int> altLens;
  altLens.resize( nAlts );
  for( int i = 0; i < nAlts; i++ ){
    altLens[i] = (int) length( altSet[i] );
  }


  // Loop over all reads in bam file(s) 1 (deletion biased calls)
  for( auto i = vais.begin(); i != vais.end(); i++ ){
    vector< double> prefs;
    prefs.resize( nAlts+1 );
    //read does not occur in bam file(s) 2 (insertion biased calls)
    if( gtm == genotype_ad or gtm == genotype_joint ){
      int bestI = (*i).second.alignmentPreference( O, prefs );
      if( bestI != NO_BEST )
	rI[bestI]++;
      rI[rI.size()-1]++;
    }

    if( gtm == genotype_va or gtm == genotype_joint ){
      int bestI = (*i).second.vaPreference( O, length( var.ref ), altLens, prefs );
      if( bestI != NO_BEST )
	VAs[bestI]++;
      VAs[VAs.size()-1]++;
      if( O.verbose ) cerr << "va " << (*i).first << " " << (*i).second.nD << " " << (*i).second.nI << " " << prefs[0] << " " << prefs[1] << " " << bestI << endl;
    }

    if( gtm == genotype_va_old ){
      int bestI = 0;
      double bestScore = 0; //std::numeric_limits< float>::max();
      for( int iP = 0; iP < nAlts; iP++ ){
	double cScore = O.overlapBits*(-(*i).second.supports( (int) length(var.ref), altLens[iP], O )*1.0 + (*i).second.rejects( (int) length(var.ref), altLens[iP] , O )*1.0); // ACHTUNG: needs fixing
	prefs[iP+1] += cScore;
	if( cScore < bestScore ){
	  bestScore = cScore;
	  bestI = iP+1;
	}
      }
      VAs[bestI]++;
      VAs[VAs.size()-1]++;
      if( O.verbose ) cerr << "va_old " << (*i).first << " " << (*i).second.nD << " " << (*i).second.nI << " " << prefs[0] << " " << prefs[1] << " " << bestI << endl;      
    }

    if( gtm == genotype_presence ){
      if( (*i).second.present( O ) ){
	prefs[0] += O.overlapBits;
      }else{
	prefs[1] += O.overlapBits;
      }
      for( int iP = 2; iP <= nAlts; iP++ ){
	prefs[iP] += O.overlapBits;
      }
    }

    //    if( O.verbose ) cerr << "Update VC bam1 " << (*i).first << " " << prefs[0] << " " << prefs[1] << " " << (*i).second.supports( (int) length( var.ref ), (int) length( altSet[0] ), O ) << " " << (*i).second.rejects( (int) length( var.ref), (int) length( altSet[0]), O ) << endl;

  
    float minPref = std::numeric_limits< float>::max(); 
    float maxPref = std::numeric_limits< float>::min(); 
    for( int iP = 0; iP < nAlts+1; iP++ ){
      if( prefs[iP] < minPref )
	minPref = prefs[iP];
      if( prefs[iP] > maxPref )
	maxPref = prefs[iP];
    }
    for( int iP = 0; iP < nAlts+1; iP++ )
      prefs[iP] -= minPref;

#define MINIMUM_PREF_DIFF 2.0

    if( maxPref-minPref > MINIMUM_PREF_DIFF ){
      int vCI = 0;
      for( int a1 = 0; a1 < nAlts+1; a1++ ){
	for( int a2 = 0; a2 <= a1; a2++ ){
	  if( a1 != a2 ){
	    if( prefs[a1] == prefs[a2] ){
	      vC[vCI] += prefs[a1];
	    }else if( prefs[a1] > prefs[a2]+2 ){
	      vC[vCI] += prefs[a2]+1;
	    }else if( prefs[a2] > prefs[a1] +2 ){
	      vC[vCI] += prefs[a1]+1;
	    }else if( prefs[a1] > prefs[a2] ){
	      vC[vCI] += (prefs[a1]+prefs[a2])/2.0;
	    }
	  }else{
	    vC[vCI] += prefs[a1];
	  }
	  vCI++;
	}
      }
    }
  }

  if( O.verbose ) cerr << "multiUpdateVC " << vC[0] << " " << vC[1] << " " << vC[2] << endl;
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
  vector< TSequence> altSeqs;
  CharString chrom;
  int32_t beginPos;

  StringSet< CharString> altSet;
  strSplit(altSet,variant.alt,EqualsChar<','>());
  
  int nAlts = (int) length( altSet );
  altSeqs.resize( nAlts );
  if( O.verbose ) cerr << "nAlts " << nAlts << endl;

  // move outside of this function
  getLocRefAlt( variant, vcfS, faiI, chrom, beginPos, refSeq, altSeqs, O );

  if( O.outputRefAlt ){
    cout << chrom << " " << beginPos + 1 << " " << variant.info << " " << refSeq;
    for( int i = 0; i < nAlts; i++ )
      cout << " " << altSeqs[i];
    cout << endl;
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
    ss2 << "variants_to_minimap_" << O.randID[omp_get_thread_num()] << ".fa";
    string str2 = ss2.str();
    minimapF.open( str2 );
    minimapF << ">Ref" << endl;
    minimapF << refSeq << endl;
    for( int i = 0; i < nAlts; i++ ){ 
      minimapF << ">Alt" << i << endl;
      minimapF << altSeqs[i] << endl;
    }
    minimapF.close();
    
    std::stringstream cmd;
    cmd.str("");
    if( O.aligner == aligner_minimap ){
      //cmd << "./minimap2 -c reads_to_minimap.fa variants_to_minimap > minimap_out.paf";   
      // PAF might be more efficient more efficient, but cost is probably mainly related to creating cigar
      cmd << "minimap2 -a -A " << O.match << " -B " << -1*O.mismatch << " -O " << -1*O.gapOpen << " -E " << -1*O.gapExtend << " variants_to_minimap_" << O.randID[omp_get_thread_num()] << ".fa reads_to_minimap_" << O.randID[omp_get_thread_num()] << ".fa > mapping_out_" << O.randID[omp_get_thread_num()] << ".sam";   
    }else if( O.aligner == aligner_bwa){
      cmd << "bwa index variants_to_minimap_" << O.randID[omp_get_thread_num()] << ".fa"; 
      if (system(cmd.str().c_str()) != 0){
	std::cerr << "ERROR in bwa index " << std::endl;
	return 7;
      } 
      cmd.str("");
      cmd << "bwa mem -A " << O.match << " -B " << -1*O.mismatch << " -O " << -1*O.gapOpen << " -E " << -1*O.gapExtend << " variants_to_minimap_" << O.randID[omp_get_thread_num()]   << ".fa reads_to_minimap_" << O.randID[omp_get_thread_num()]  << ".fa > mapping_out_" << O.randID[omp_get_thread_num()]  << ".sam";   
    }else{
      cout << "Unknown aligner " << endl;
      exit(-1);
    }
    if (system(cmd.str().c_str()) != 0){
      std::cerr << "ERROR in bwa " << std::endl;
      return 7;
    } 

    stringstream ss;
    ss << "mapping_out_" << O.randID[omp_get_thread_num()]  << ".sam";
    string str = ss.str();

    BamFileIn bamFileIn(str.c_str());

    // Get rid of header
    BamHeader header;
    readHeader(header, bamFileIn);

    int refID = 0;
    if (!getIdByName( refID, contigNamesCache(context( bamFileIn )), "Ref" ) and O.verbose)
      cerr << "getIdByName failed, could not find Ref" << endl;
    vector<int> altIDs;
    altIDs.resize( nAlts );
    //    int altID = 0;
    for( int i = 0; i < nAlts; i++ ){
      std::ostringstream buff;
      buff << "Alt" << i;
      if (!getIdByName( altIDs[i], contigNamesCache(context( bamFileIn )), buff.str() ) and O.verbose)
	cerr << "getIdByName failed, could not find Alt" << endl;
      if( O.verbose ) cerr << "Got ids by name Ref " << refID << " Alt" << i << " " << altIDs[i] << endl;
    }
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
	  if( (record.rID == refID) && (vais[record.qName].alignS[0] == NO_ALIGNMENT || vais[record.qName].alignS[0] < tagValInt) )
	    vais[record.qName].alignS[0] = (double) tagValInt;
	  else{
	    for( int i = 0; i < nAlts; i++ ){
	      if( (record.rID == altIDs[i]) && (vais[record.qName].alignS[i+1] == NO_ALIGNMENT || vais[record.qName].alignS[i+1] < tagValInt) )
		vais[record.qName].alignS[i+1] = (double) tagValInt;
	    }
	  }
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
      vais[(*i).second.qName].alignS[0] = (int) globalAlignment(align, scoringScheme, AlignConfig< true, true, true, true>());
      for( int iAlt = 0; iAlt < nAlts; iAlt++ ){
	assignSource(row(align,0), altSeqs[iAlt]);
	vais[(*i).second.qName].alignS[iAlt+1] = (int) globalAlignment(align, scoringScheme, AlignConfig< true, true, true, true>());
	if( O.verbose ) cerr << "Aligned " << (*i).second.qName << " " << vais[(*i).second.qName].alignS[0]  << " " << vais[(*i).second.qName].alignS[iAlt+1] << " " << altSeqs[iAlt] << endl;
      }
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
  O.randID.resize( num_threads );
  for( int i = 0; i < num_threads; i++ ){
    O.randID[i] = e();
  }

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
 
    StringSet< CharString> altSet;
    strSplit(altSet,records[nR].alt,EqualsChar<','>());
    int nAlleles = (int) length( altSet ) + 1;
    
    CharString varChrom;  
    varChrom = contigNames(context(vcfIn))[records[nR].rID];

    // Implemented as a vector as we may have more than one output per marker
    std::vector< std::vector<double> > vC; //variant calls
    std::vector< std::vector<int> > AD;  // Allele depth counts
    std::vector< std::vector<int> > VA;  // Variant Alignment counts
    vC.resize(1);
    AD.resize(1);
    VA.resize(1);
    if( O.gtModel == genotype_multi ){
      vC.resize( 5 );
      AD.resize( 5 );
      VA.resize( 5 );
    }
    for( int mI = 0; mI < (int) vC.size(); mI++ ){
      vC[mI].resize(nAlleles*(nAlleles+1)/2);
      AD[mI].resize(nAlleles+1);
      VA[mI].resize(nAlleles+1);
    }
    for( int i = 0; i < (int) bamIndexHandles[omp_get_thread_num()].size(); i++ ){
      std::map< CharString, BamAlignmentRecord> bars;
      readBamRegion( bamIndexHandles[omp_get_thread_num()][i], (*bamHandles[omp_get_thread_num()][i]), varChrom, records[nR], O, bars, vaisR );
      if( O.verbose ) cerr << "Finished reading bam region bam1 " << endl;
      LRprocessReads(records[nR], vcfIn, faIndexHandles[omp_get_thread_num()], bars, O, vaisR );
    }
    if( useBam2 ){
      for( int i = 0; i < (int) bamIndex2Handles[omp_get_thread_num()].size(); i++ ){
	std::map< CharString, BamAlignmentRecord> bars;
	readBamRegion( bamIndex2Handles[omp_get_thread_num()][i], (*bam2Handles[omp_get_thread_num()][i]), varChrom, records[nR], O, bars, vaisR );
	LRprocessReads(records[nR], vcfIn, faIndexHandles[omp_get_thread_num()], bars, O, vaisR );
      }
    }
    if( O.gtModel == genotype_multi ){
      multiUpdateVC( records[nR], vaisR, vC[0], AD[0], VA[0], O, genotype_ad );
      multiUpdateVC( records[nR], vaisR, vC[1], AD[1], VA[1], O, genotype_va );
      multiUpdateVC( records[nR], vaisR, vC[2], AD[2], VA[2], O, genotype_joint );
      multiUpdateVC( records[nR], vaisR, vC[3], AD[3], VA[3], O, genotype_presence );
      multiUpdateVC( records[nR], vaisR, vC[4], AD[4], VA[4], O, genotype_va_old );
    }else{
      multiUpdateVC( records[nR], vaisR, vC[0], AD[0], VA[0], O, O.gtModel );
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


