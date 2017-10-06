/*
  John M. Gaspar (jsh58@wildcats.unh.edu)
  April 2015 (updated July 2016, Apr. 2017, June 2017)

  Analyzing paired-end reads for overlaps. Two modes:
  - 'stitch': producing a single, merged read for reads
     with sufficient overlaps
  - 'adapter-removal': removing adapters (3' overhangs
     of stitched alignment) from individual reads
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <getopt.h>
#include <zlib.h>
#include <omp.h>
#include "stitch_0_7.h"

/* void printVersion()
 * Print version and copyright.
 */
void printVersion(void) {
  fprintf(stderr, "NGmerge, version %s\n", VERSION);
  fprintf(stderr, "Copyright (C) 2017 John M. Gaspar (jsh58@wildcats.unh.edu)\n");
  exit(-1);
}

/* void usage()
 * Prints usage information.
 */
void usage(void) {
  fprintf(stderr, "Usage: ./NGmerge {-%c <file> -%c <file>", FIRST, SECOND);
  fprintf(stderr, " -%c <file>}  [optional arguments]\n", OUTFILE);
  fprintf(stderr, "Required arguments:\n");
  fprintf(stderr, "  -%c  <file>       Input FASTQ file with reads from forward direction\n", FIRST);
  fprintf(stderr, "  -%c  <file>       Input FASTQ file with reads from reverse direction\n", SECOND);
  fprintf(stderr, "  -%c  <file>       Output FASTQ file(s):\n", OUTFILE);
  fprintf(stderr, "                   - in 'stitch' mode (def.), the file of merged reads\n");
  fprintf(stderr, "                   - in 'adapter-removal' mode (-%c), the output files\n", ADAPTOPT);
  fprintf(stderr, "                     will be <file>%s and <file>%s\n", ONEEXT, TWOEXT);
  fprintf(stderr, "Alignment parameters:\n");
  fprintf(stderr, "  -%c  <int>        Minimum overlap of the paired-end reads (def. %d)\n", OVERLAP, DEFOVER);
  fprintf(stderr, "  -%c  <float>      Mismatches to allow in the overlapped region\n", MISMATCH);
  fprintf(stderr, "                     (a fraction of the overlap length; def. %.2f)\n", DEFMISM);
  fprintf(stderr, "  -%c               Use 'adapter-removal' mode (also sets -%c option)\n", ADAPTOPT, DOVEOPT);
  fprintf(stderr, "  -%c               Option to check for dovetailing (with 3' overhangs)\n", DOVEOPT);
  fprintf(stderr, "  -%c  <int>        Minimum overlap of dovetailed alignments (def. %d)\n", DOVEOVER, DEFDOVE);
  fprintf(stderr, "  -%c               Option to produce shortest stitched read\n", MAXOPT);
  fprintf(stderr, "I/O options:\n");
  fprintf(stderr, "  -%c  <file>       Log file for stitching results of each read pair\n", LOGFILE);
  fprintf(stderr, "  -%c  <file>       FASTQ files for reads that failed stitching\n", UNFILE);
  fprintf(stderr, "                     (output as <file>%s and <file>%s)\n", ONEEXT, TWOEXT);
  fprintf(stderr, "  -%c  <file>       Log file for dovetailed reads (adapter sequences)\n", DOVEFILE);
  fprintf(stderr, "  -%c  <file>       Log file for formatted alignments of merged reads\n", ALNFILE);
  fprintf(stderr, "  -%c/-%c            Option to gzip (-%c) or not (-%c) FASTQ output(s)\n", GZOPT, UNGZOPT, GZOPT, UNGZOPT);
  fprintf(stderr, "  -%c               Option to produce interleaved FASTQ output(s)\n", INTEROPT);
  fprintf(stderr, "  -%c               Option to use 'fastq-join' method for quality scores\n", FJOINOPT);
  fprintf(stderr, "  -%c  <int>        FASTQ quality offset (def. %d)\n", QUALITY, OFFSET);
  fprintf(stderr, "  -%c  <int>        Maximum input quality score (0-based; def. %d)\n", SETQUAL, MAXQUAL);
  fprintf(stderr, "  -%c  <int>        Number of threads to use (def. %d)\n", THREADS, DEFTHR);
  fprintf(stderr, "  -%c               Option to print status updates/counts to stderr\n", VERBOSE);
  exit(-1);
}

/* int error()
 * Prints an error message.
 */
int error(char* msg, enum errCode err) {
  fprintf(stderr, "Error! %s%s\n", msg, errMsg[err]);
  return -1;
}

/* void* memalloc()
 * Allocates a heap block.
 */
void* memalloc(int size) {
  void* ans = malloc(size);
  if (ans == NULL)
    exit(error("", ERRMEM));
  return ans;
}

/* float getFloat(char*)
 * Converts the given char* to a float.
 */
float getFloat(char* in) {
  char** endptr = NULL;
  float ans = strtof(in, endptr);
  if (endptr != '\0')
    exit(error(in, ERRFLOAT));
  return ans;
}

/* int getInt(char*)
 * Converts the given char* to an int.
 */
int getInt(char* in) {
  char** endptr = NULL;
  int ans = (int) strtol(in, endptr, 10);
  if (endptr != '\0')
    exit(error(in, ERRINT));
  return ans;
}

/* char rc(char)
 * Returns the complement of the given base.
 */
char rc(char in) {
  char out;
  if (in == 'A') out = 'T';
  else if (in == 'T') out = 'A';
  else if (in == 'C') out = 'G';
  else if (in == 'G') out = 'C';
  else if (in == 'N') out = 'N';
  else {
    char msg[4] = "' '";
    msg[1] = in;
    exit(error(msg, ERRUNK));
  }
  return out;
}

/* char* getLine()
 * Reads the next line from a file.
 */
char* getLine(char* line, int size, File in, bool gz) {
  if (gz)
    return gzgets(in.gzf, line, size);
  else
    return fgets(line, size, in.f);
}

/* void checkHeaders()
 * Ensure headers match (up to first space character);
 *   create consensus header.
 */
void checkHeaders(char* head1, char* head2, char* header) {
  bool ok = false;  // match boolean
  int j;
  for (j = 0; head1[j] != '\n' && head1[j] != '\0'; j++) {
    if (head1[j] != head2[j]) {
      if (ok)
        break;
      for ( ; head1[j] != '\n' && head1[j] != '\0'
        && head1[j] != ' '; j++) ;
      head1[j] = '\0';  // trim head1 for err msg
      exit(error(head1, ERRHEAD));
    } else if (head1[j] == ' ')
      ok = true;  // headers match
    header[j] = head1[j];
  }
  if (header[j - 1] == ' ')
    header[j - 1] = '\0'; // removing trailing space
  else
    header[j] = '\0';
}

/* void checkQual()
 * Check given quality string for offset errors.
 */
void checkQual(char* qual, int len, int offset,
    int maxQual) {
  for (int i = 0; i < len; i++)
    // error if qual < 0 or qual > maxQual
    if (qual[i] < offset || qual[i] > offset + maxQual)
      exit(error("", ERROFFSET));
}

/* void processSeq()
 * Process the given sequence; save length;
 *   for 2nd read, save reversed seq/qual.
 */
void processSeq(char** read, int* len, bool i,
    int j, int offset, int maxQual) {

  // remove new-line character and save length
  int k;
  for (k = 0; read[j][k] != '\n' && read[j][k] != '\0'; k++) ;
  read[j][k] = '\0';
  if (j == SEQ)
    *len = k;  // save read length
  else if (k != *len)
    exit(error("", ERRQUAL)); // seq/qual length mismatch

  // for 2nd read (i == true), save revComp(seq) or rev(qual)
  if (i) {
    int dest = j + EXTRA; // save to 'extra' field of read2
    int m = 0;
    if (j == SEQ) {
      dest++;  // increment b/c of fastq 'plus' line
      for (k--; k > -1; k--)
        read[dest][m++] = rc(read[j][k]);
    } else
      for (k--; k > -1; k--)
        read[dest][m++] = read[j][k];
    read[dest][m] = '\0';
  } else if (j == SEQ)
    // check 1st read's sequence for non-ACGTN chars
    for (int m = 0; m < k; m++)
      rc(read[j][m]);

  // check quality scores
  if (j == QUAL)
    checkQual(read[j], k, offset, maxQual);
}

/* bool loadReads()
 * Load a pair of reads. Check formatting, determine
 *   consensus header. Return false on EOF.
 */
bool loadReads(File in1, File in2, char** read1, char** read2,
    char* header, int* len1, int* len2, int offset,
    int maxQual, bool gz1, bool gz2) {

  // load both reads from input files (LOCK)
  bool flag = false;  // boolean for EOF
  #pragma omp critical
  for (int i = 0; i < 2; i++) {
    File in = in1;
    char** read = read1;
    bool gz = gz1;
    if (i) {
      in = in2;
      read = read2;
      gz = gz2;
    }

    // load read (4 lines)
    for (int j = 0; j < FASTQ; j++)
      if (getLine(read[j], MAX_SIZE, in, gz) == NULL) {
        if (j == 0) {
          if (i == 0) {
            flag = true;  // EOF
            break;
          } else {
            int k = 0;
            for ( ; read1[HEAD][k] != '\n' && read1[HEAD][k] != '\0'
              && read1[HEAD][k] != ' '; k++) ;
            read1[HEAD][k] = '\0';  // trim header for err msg
            exit(error(read1[HEAD], ERRHEAD));
          }
        } else
          exit(error("", ERRSEQ));
      }
    if (flag)
      break;

  }  // (UNLOCK)

  if (flag)
    return false;  // EOF

  // check fastq formatting
  if (read1[HEAD][0] != BEGIN || read1[PLUS][0] != PLUSCHAR
      || read2[HEAD][0] != BEGIN || read2[PLUS][0] != PLUSCHAR)
    exit(error("", ERRFASTQ));

  // process sequence/quality lines
  processSeq(read1, len1, false, SEQ, offset, maxQual);
  processSeq(read1, len1, false, QUAL, offset, maxQual);
  processSeq(read2, len2, true, SEQ, offset, maxQual);
  processSeq(read2, len2, true, QUAL, offset, maxQual);

  // check headers
  checkHeaders(read1[HEAD], read2[HEAD], header);

  return true;
}

/* float compare()
 * Compare two sequences. Return the fraction mismatch.
 */
float compare(char* seq1, char* seq2, int length,
    float mismatch, int overlap) {
  int mis = 0;       // number of mismatches
  int len = length;  // length of overlap, not counting Ns
  float allow = len * mismatch;
  for (int i = 0; i < length; i++) {
    // do not count Ns
    if (seq1[i] == 'N' || seq2[i] == 'N') {
      if (--len < overlap || mis > len * mismatch)
        return NOTMATCH;
      allow = len * mismatch;
    } else if (seq1[i] != seq2[i] && ++mis > allow)
      return NOTMATCH;
  }
  return (float) mis / len;
}

/* int findPos()
 * Find optimal overlapping position.
 *   Currently, quality scores are not considered
 *   (e.g. decreased penalty for a low-quality mismatch).
 */
int findPos (char* seq1, char* seq2, char* qual1,
    char* qual2, int len1, int len2, int overlap,
    bool dovetail, int doveOverlap, float mismatch,
    bool maxLen, float* best) {

  // check for regular (non-dovetailed) alignments
  int pos = len1 - overlap + 1;  // position of match
  int i = len1 - overlap;
  for ( ; i > -1 && len1 - i <= len2; i--) {
    // align sequences
    float res = compare(seq1 + i, seq2, len1 - i,
      mismatch, overlap);

    // compare result
    if (res < *best || (res == *best && !maxLen)) {
      *best = res;
      pos = i;
    }
    if (res == 0.0f && maxLen)
      return pos;  // shortcut for exact match
  }

  // check for dovetailing
  if (dovetail) {

    // if no regular alignment, reset i
    if (i == len1 - overlap)
      i = (len1 > len2 ? len1 - len2 - 1 : -1);

    // continue decrementing i
    for ( ; ; i--) {
      float res = NOTMATCH;
      if (i >= 0) {
        // read1 is longer, with 3' overhang
        if (len2 < doveOverlap)
          break;
        res = compare(seq1 + i, seq2, len2,
          mismatch, doveOverlap);

      } else if (len1 < len2 + i) {
        // read2 has 3' overhang, read1 determines overlap
        if (len1 < doveOverlap)
          break;
        res = compare(seq1, seq2 - i, len1,
          mismatch, doveOverlap);

      } else {
        // read2 has 3' overhang and determines overlap
        if (len2 + i < doveOverlap)
          break;
        res = compare(seq1, seq2 - i, len2 + i,
          mismatch, doveOverlap);
      }

      // compare result
      if (res < *best || (res == *best && !maxLen)) {
        *best = res;
        pos = i;
      }
      if (res == 0.0f && maxLen)
        return pos;  // shortcut for exact match
    }
  }

  return pos;
}

/* void printDove()
 * Log 3' overhangs of dovetailed reads.
 */
void printDove(File dove, char* header, char** read1,
    char** read2, int len1, int len2, int pos,
    omp_lock_t* lock) {
  if (len1 > len2 + pos || pos < 0) {
    omp_set_lock(lock);
    fprintf(dove.f, "%s\t%s\t%s\n", header + 1,
      len1 > len2 + pos ? read1[SEQ] + len2 + pos : "-",
      pos < 0 ? read2[SEQ] + len2 + pos : "-");
    omp_unset_lock(lock);
  }
}

/* void printGZNoAdapt()
 * Print the reads minus adapters (gzip output).
 */
void printGZNoAdapt(gzFile out1, gzFile out2,
    char** read1, char** read2, int end1, int end2) {

  // print fwd read
  gzprintf(out1, "%s", read1[HEAD]);
  for (int i = 0; i < end1; i++)
    gzputc(out1, read1[SEQ][i]);
  gzprintf(out1, "\n%s", read1[PLUS]);
  for (int i = 0; i < end1; i++)
    gzputc(out1, read1[QUAL][i]);
  gzputc(out1, '\n');

  // print rev read
  gzprintf(out2, "%s", read2[HEAD]);
  for (int i = 0; i < end2; i++)
    gzputc(out2, read2[SEQ][i]);
  gzprintf(out2, "\n%s", read2[PLUS]);
  for (int i = 0; i < end2; i++)
    gzputc(out2, read2[QUAL][i]);
  gzputc(out2, '\n');
}

/* void printNoAdapt()
 * Print the reads minus adapters.
 */
void printNoAdapt(FILE* out1, FILE* out2, char** read1,
    char** read2, int end1, int end2) {

  // print fwd read
  fprintf(out1, "%s", read1[HEAD]);
  for (int i = 0; i < end1; i++)
    fputc(read1[SEQ][i], out1);
  fprintf(out1, "\n%s", read1[PLUS]);
  for (int i = 0; i < end1; i++)
    fputc(read1[QUAL][i], out1);
  fputc('\n', out1);

  // print rev read
  fprintf(out2, "%s", read2[HEAD]);
  for (int i = 0; i < end2; i++)
    fputc(read2[SEQ][i], out2);
  fprintf(out2, "\n%s", read2[PLUS]);
  for (int i = 0; i < end2; i++)
    fputc(read2[QUAL][i], out2);
  fputc('\n', out2);
}

/* bool printResAdapt()
 * Control printing of reads minus adapters.
 *   Return 1 if adapter found, else 0.
 */
bool printResAdapt(File out1, File out2, File dove,
    bool doveOpt, char* header, char** read1, char** read2,
    int len1, int len2, int pos, float best, bool gz,
    omp_lock_t* lock) {

  bool adapter = false;
  int end1 = len1;
  int end2 = len2;

  // if found, identify locations of adapters
  if (len1 > len2 + pos || pos < 0) {
    adapter = true;
    if (len1 > len2 + pos)
      end1 = len2 + pos;
    if (pos < 0)
      end2 += pos;
    if (doveOpt)
      printDove(dove, header, read1, read2,
        len1, len2, pos, lock + DOVE);
  }

  // print output
  omp_set_lock(lock + OUT);
  if (gz)
    printGZNoAdapt(out1.gzf, out2.gzf, read1, read2,
      end1, end2);
  else
    printNoAdapt(out1.f, out2.f, read1, read2,
      end1, end2);
  omp_unset_lock(lock + OUT);

  return adapter;
}

/* void printAln2()
 * Printing details of stitch mismatches.
 */
void printAln2(File aln, char* header, char** read1,
    char** read2, int len1, int len2, int pos) {
  int i = pos;
  int j = 0;
  if (pos < 0) {
    j = -pos;
    i = 0;
  }
  while (i < len1 && j < len2) {
    if (read1[SEQ][i] == 'N' || read2[SEQ + EXTRA + 1][j] == 'N'
        || read1[SEQ][i] != read2[SEQ + EXTRA + 1][j])
      fprintf(aln.f, "%s\t%d\t%c\t%c\t%c\t%c\n",
        header + 1, i, read1[SEQ][i], read1[QUAL][i],
        read2[SEQ + EXTRA + 1][j], read2[QUAL + EXTRA][j]);
    i++;
    j++;
  }
}

/* void printAln()
 * Print nicely formatted alignment of stitched reads.
 */
void printAln(File aln, char* header, char** read1,
    char** read2, int len1, int len2, int pos) {
  fprintf(aln.f, "%s\n", header + 1);

  // print sequence alignment
  fprintf(aln.f, "seq_R1:  ");
  for (int i = 0; i > pos; i--)
    fputc(' ', aln.f);
  fprintf(aln.f, "%s\n", read1[SEQ]);

  // print '|' for matches, ':' for Ns
  fprintf(aln.f, "         ");
  int i;
  for (i = 0; i < abs(pos); i++)
    fputc(' ', aln.f);
  int j = 0;
  if (pos < 0) {
    j = -pos;
    i = 0;
  }
  while (i < len1 && j < len2) {
    fputc((read1[SEQ][i] == 'N' || read2[SEQ + EXTRA + 1][j] == 'N') ?
      ':' : (read1[SEQ][i] == read2[SEQ + EXTRA + 1][j] ?
      '|' : ' '), aln.f);
    i++;
    j++;
  }
  fputc('\n', aln.f);

  fprintf(aln.f, "seq_R2:  ");
  for (int i = 0; i < pos; i++)
    fputc(' ', aln.f);
  fprintf(aln.f, "%s\n\n", read2[SEQ + EXTRA + 1]);

  // print quality scores
  fprintf(aln.f, "qual_R1: ");
  for (int i = 0; i > pos; i--)
    fputc(' ', aln.f);
  fprintf(aln.f, "%s\n", read1[QUAL]);
  fprintf(aln.f, "qual_R2: ");
  for (int i = 0; i < pos; i++)
    fputc(' ', aln.f);
  fprintf(aln.f, "%s\n\n", read2[QUAL + EXTRA]);
}

/* void createSeq2()
 * Create stitched sequence (into seq1, qual1).
 *   Use empirical error profiles for quality scores.
 */
void createSeq2(char* seq1, char* seq2, char* qual1, char* qual2,
    int len1, int len2, int pos, int offset) {
  int len = len2 + pos;  // length of stitched sequence
  for (int i = 0; i < len; i++) {
    if (i - pos < 0)
      // 1st read only: continue
      continue;
    else if (i >= len1) {
      // 2nd read only: copy seq and qual
      seq1[i] = seq2[i-pos];
      qual1[i] = qual2[i-pos];
    } else if (seq2[i-pos] == 'N')
      // 2nd read 'N': continue
      continue;
    else if (seq1[i] == 'N') {
      // 1st read 'N': copy seq and qual
      seq1[i] = seq2[i-pos];
      qual1[i] = qual2[i-pos];
    } else if (seq1[i] != seq2[i-pos]) {
      // mismatch:
      //   - base matches higher quality score or equal
      //     quality score that is closer to 5' end
      //   - quality score copied from mism array
      if (qual1[i] < qual2[i-pos] ||
          (qual1[i] == qual2[i-pos] && i >= len2 - i + pos) )
        seq1[i] = seq2[i-pos];
      qual1[i] = mism[ (int) qual1[i] - offset ]
        [ (int) qual2[i-pos] - offset ] + offset;
    } else
      // match:
      //   - quality score copied from match array
      qual1[i] = match[ (int) qual1[i] - offset ]
        [ (int) qual2[i-pos] - offset ] + offset;

  }
  seq1[len] = '\0';
  qual1[len] = '\0';
}

/* void createSeq()
 * Create stitched sequence (into seq1, qual1).
 *   Use 'fastq-join' method for quality scores.
 */
void createSeq(char* seq1, char* seq2, char* qual1, char* qual2,
    int len1, int len2, int pos, int offset) {
  int len = len2 + pos;  // length of stitched sequence
  for (int i = 0; i < len; i++) {
    if (i - pos < 0)
      // 1st read only: continue
      continue;
    else if (i >= len1) {
      // 2nd read only: copy seq and qual
      seq1[i] = seq2[i-pos];
      qual1[i] = qual2[i-pos];
    } else if (seq1[i] != seq2[i-pos]) {
      // disagreements favor higher quality score or
      //   equal quality score that is not an 'N' or
      //   is closer to 5' end; reduce qual scores
      if (qual1[i] < qual2[i-pos] ||
          (qual1[i] == qual2[i-pos] &&
          ((i >= len2 - i + pos && seq2[i-pos] != 'N')
          || seq1[i] == 'N') ) ) {
        qual1[i] = (seq1[i] == 'N' ? qual2[i-pos]
          : qual2[i-pos] - qual1[i] + offset);
        seq1[i] = seq2[i-pos];
      } else if (seq2[i-pos] != 'N')
        qual1[i] -= qual2[i-pos] - offset;
    } else if (qual1[i] < qual2[i-pos])
      // seq agreement: use higher qual score
      qual1[i] = qual2[i-pos];

  }
  seq1[len] = '\0';
  qual1[len] = '\0';
}

/* void printRes()
 * Print stitched read.
 */
void printRes(File out, File log, bool logOpt, File dove,
    bool doveOpt, File aln, int alnOpt, char* header,
    char** read1, char** read2, int len1, int len2,
    int pos, float best, int offset, bool gz, bool fjoin,
    omp_lock_t* lock) {
  // log result
  if (logOpt) {
    omp_set_lock(lock + LOG);
    fprintf(log.f, "%s\t%d\t%d\t", header + 1,
      pos < 0 ? (len2+pos < len1 ? len2+pos : len1) :
      (len1-pos < len2 ? len1-pos : len2), len2 + pos);
    best ? fprintf(log.f, "%.3f", best) : fprintf(log.f, "0");
    fprintf(log.f, "\n");
    omp_unset_lock(lock + LOG);
  }
  if (doveOpt)
    printDove(dove, header, read1, read2, len1, len2,
      pos, lock + DOVE);
  if (alnOpt) {
    omp_set_lock(lock + ALN);
    if (alnOpt == 1)
      printAln(aln, header, read1, read2, len1, len2, pos);
    else if (alnOpt == 2)
      printAln2(aln, header, read1, read2, len1, len2, pos);
    omp_unset_lock(lock + ALN);
  }

  // print stitched sequence
  if (fjoin)
    createSeq(read1[SEQ], read2[SEQ + EXTRA + 1], read1[QUAL],
      read2[QUAL + EXTRA], len1, len2, pos, offset);
  else
    createSeq2(read1[SEQ], read2[SEQ + EXTRA + 1], read1[QUAL],
      read2[QUAL + EXTRA], len1, len2, pos, offset);
  omp_set_lock(lock + OUT);
  if (gz)
    gzprintf(out.gzf, "%s\n%s\n+\n%s\n", header,
      read1[SEQ], read1[QUAL]);
  else
    fprintf(out.f, "%s\n%s\n+\n%s\n", header,
      read1[SEQ], read1[QUAL]);
  omp_unset_lock(lock + OUT);

  // print to alignment output too
  if (alnOpt == 1) {
    omp_set_lock(lock + ALN);
    fprintf(aln.f, "merged\nseq:     ");
    for (int i = 0; i > pos; i--)
      fputc(' ', aln.f);
    fprintf(aln.f, "%s\n", read1[SEQ]);
    fprintf(aln.f, "qual:    ");
    for (int i = 0; i > pos; i--)
      fputc(' ', aln.f);
    fprintf(aln.f, "%s\n\n\n", read1[QUAL]);
    omp_unset_lock(lock + ALN);
  }
}

/* void printFail()
 * Print stitch failure reads.
 */
void printFail(File un1, File un2, bool unOpt,
    File log, bool logOpt, char* header, char** read1,
    char** read2, bool gz, omp_lock_t* outLock,
    omp_lock_t* logLock) {
  if (logOpt) {
    omp_set_lock(logLock);
    fprintf(log.f, "%s\t%s\n", header + 1, NA);
    omp_unset_lock(logLock);
  }
  if (unOpt) {
    omp_set_lock(outLock);
    if (gz) {
      gzprintf(un1.gzf, "%s%s\n%s%s\n", read1[HEAD],
        read1[SEQ], read1[PLUS], read1[QUAL]);
      gzprintf(un2.gzf, "%s%s\n%s%s\n", read2[HEAD],
        read2[SEQ], read2[PLUS], read2[QUAL]);
    } else {
      fprintf(un1.f, "%s%s\n%s%s\n", read1[HEAD],
        read1[SEQ], read1[PLUS], read1[QUAL]);
      fprintf(un2.f, "%s%s\n%s%s\n", read2[HEAD],
        read2[SEQ], read2[PLUS], read2[QUAL]);
    }
    omp_unset_lock(outLock);
  }
}

/* int readFile()
 * Analyzes the reads in a set of input files.
 *   Controls writing to the output file(s).
 *   Multithreaded.
 */
int readFile(File in1, File in2, File out, File out2,
    File un1, File un2, bool unOpt, File log,
    bool logOpt, int overlap, bool dovetail, int doveOverlap,
    File dove, bool doveOpt, File aln, int alnOpt,
    bool adaptOpt, float mismatch, bool maxLen,
    int* stitch, int offset, int maxQual,
    bool gz1, bool gz2, bool gzOut, bool fjoin,
    int threads) {

  // initialize omp locks -- out, un, log, dove, aln
  omp_lock_t lock[OMP_LOCKS];
  for (int i = 0; i < OMP_LOCKS; i++)
    omp_init_lock(&lock[i]);

  // process files in parallel
  int count = 0, stitchRed = 0;
  #pragma omp parallel num_threads(threads) reduction(+: count, stitchRed)
  {

    // allocate memory for both reads
    char** read1 = (char**) memalloc(FASTQ * sizeof(char*));
    char** read2 = (char**) memalloc((FASTQ + EXTRA) * sizeof(char*));
    for (int i = 0; i < FASTQ + EXTRA; i++) {
      if (i < FASTQ)
        read1[i] = (char*) memalloc(MAX_SIZE);
      // for 2nd read, save extra fields for revComp(seq) and rev(qual)
      read2[i] = (char*) memalloc(MAX_SIZE);
    }
    char* header = (char*) memalloc(MAX_SIZE); // consensus header

    // process reads
    int len1 = 0, len2 = 0; // lengths of reads
    while (loadReads(in1, in2, read1, read2, header,
        &len1, &len2, offset, maxQual, gz1, gz2)) {

      // find optimal overlap
      float best = 1.0f;
      int pos = findPos(read1[SEQ], read2[SEQ + EXTRA + 1],
        read1[QUAL], read2[QUAL + EXTRA], len1, len2, overlap,
        dovetail, doveOverlap, mismatch, maxLen, &best);

      // print result
      if (pos == len1 - overlap + 1) {
        // stitch failure
        if (adaptOpt)
          printFail(out, out2, 1, log, 0, header, read1,
            read2, gzOut, lock + OUT, lock + LOG);
        else
          printFail(un1, un2, unOpt, log, logOpt, header,
            read1, read2, gzOut, lock + UN, lock + LOG);
      } else {
        // stitch success
        if (adaptOpt) {
          stitchRed += printResAdapt(out, out2, dove, doveOpt,
            header, read1, read2, len1, len2, pos, best,
            gzOut, lock);
        } else {
          printRes(out, log, logOpt, dove, doveOpt, aln, alnOpt,
            header, read1, read2, len1, len2, pos, best, offset,
            gzOut, fjoin, lock);
          stitchRed++;
        }
      }

      count++;
    }

    // free memory
    free(header);
    for (int i = 0; i < FASTQ + EXTRA; i++) {
      if (i < FASTQ)
        free(read1[i]);
      free(read2[i]);
    }
    free(read1);
    free(read2);

  }  // END parallel

  // destroy omp locks
  for (int i = 0; i < 5; i++)
    omp_destroy_lock(&lock[i]);

  *stitch = stitchRed;
  return count;
}

/* void openWrite()
 * Open a file for writing (stdout if file is '-').
 */
void openWrite(char* outFile, File* out, bool gz) {
  if (outFile[0] == '-' && strlen(outFile) > 1)
    exit(error(outFile, ERRNAME));
  if (gz) {
    if (!strcmp(outFile + strlen(outFile) - strlen(GZEXT), GZEXT)
        || !strcmp(outFile, "/dev/null"))
      out->gzf = gzopen(outFile, "w");
    else if (!strcmp(outFile, "-"))
      out->gzf = gzdopen(fileno(stdout), "wb");
    else {
      // add ".gz" to outFile
      char* outFile2 = memalloc(strlen(outFile)
        + strlen(GZEXT) + 1);
      strcpy(outFile2, outFile);
      strcat(outFile2, GZEXT);
      out->gzf = gzopen(outFile2, "w");
      free(outFile2);
    }
    if (out->gzf == NULL)
      exit(error(outFile, ERROPENW));
  } else {
    out->f = (strcmp(outFile, "-") ?
      fopen(outFile, "w") : stdout);
    if (out->f == NULL)
      exit(error(outFile, ERROPENW));
  }
}

/* void openFiles()
 * Opens output files for the program,
 *   adjusting file names/extensions as needed.
 */
void openFiles(char* outFile, File* out, File* out2,
    char* unFile, File* un1, File* un2,
    char* logFile, File* log,
    char* doveFile, File* dove, bool dovetail,
    char* alnFile, File* aln,
    bool adaptOpt, bool gz, bool interOpt) {

  if (adaptOpt) {
    if (interOpt)
      openWrite(outFile, out, gz);
    else if (! strcmp(outFile, "-"))
      exit(error("stdout + \"_1.fastq\"", ERROPENW));
    else if (! strcmp(outFile, "/dev/null")) {
      openWrite(outFile, out, gz);
      openWrite(outFile, out2, gz);
    } else {
      // add "_1.fastq" and "_2.fastq" extensions
      int add = strlen(ONEEXT) > strlen(TWOEXT) ?
        strlen(ONEEXT) + 1 : strlen(TWOEXT) + 1;
      char* outFile2 = memalloc(strlen(outFile) + add);
      strcpy(outFile2, outFile);
      strcat(outFile2, ONEEXT);
      openWrite(outFile2, out, gz);
      strcpy(outFile2, outFile);
      strcat(outFile2, TWOEXT);
      openWrite(outFile2, out2, gz);
      free(outFile2);
    }

  } else {
    openWrite(outFile, out, gz);

    // open optional files
    if (unFile != NULL) {
      if (interOpt)
        openWrite(unFile, un1, gz);
      else if (! strcmp(unFile, "-"))
        exit(error("stdout + \"_1.fastq\"", ERROPENW));
      else {
        // add "_1.fastq" and "_2.fastq" extensions
        int add = strlen(ONEEXT) > strlen(TWOEXT) ?
          strlen(ONEEXT) + 1 : strlen(TWOEXT) + 1;
        char* unFile2 = memalloc(strlen(unFile) + add);
        strcpy(unFile2, unFile);
        strcat(unFile2, ONEEXT);
        openWrite(unFile2, un1, gz);
        strcpy(unFile2, unFile);
        strcat(unFile2, TWOEXT);
        openWrite(unFile2, un2, gz);
        free(unFile2);
      }
    }
    if (logFile != NULL) {
      openWrite(logFile, log, false);
      fprintf(log->f, "Read\tOverlapLen\tStitchedLen\tMismatch\n");
    }
    if (alnFile != NULL)
      openWrite(alnFile, aln, false);
  }

  if (dovetail && doveFile != NULL) {
    openWrite(doveFile, dove, false);
    fprintf(dove->f, "Read\tAdapter_R1\tAdapter_R2\n");
  }
}

/* bool openRead()
 * Open a file for reading (stdin if file is '-').
 *   Return 1 if gzip compressed.
 */
bool openRead(char* inFile, File* in) {

  // open file or stdin
  bool stdinBool = (strcmp(inFile, "-") ? false : true);
  FILE* dummy = (stdinBool ? stdin : fopen(inFile, "r"));
  if (dummy == NULL)
    exit(error(inFile, ERROPEN));

  // check for gzip compression: magic number 0x1F, 0x8B
  bool gzip = true;
  int save = 0;  // first char to pushback (for stdin)
  int i, j;
  for (i = 0; i < 2; i++) {
    j = fgetc(dummy);
    if (j == EOF)
      exit(error(inFile, ERROPEN));
    if ( (i && (unsigned char) j != 0x8B)
        || (! i && (unsigned char) j != 0x1F) ) {
      gzip = false;
      break;
    }
    if (! i)
      save = j;
  }

  // for stdin, push back chars
  if (stdinBool) {
    if (gzip)
      exit(error("", ERRGZIP));
    if (ungetc(j, dummy) == EOF)
      exit(error("", ERRUNGET));
    if (i && ungetc(save, dummy) == EOF)
      exit(error("", ERRUNGET));
  }

  // open file
  if (gzip) {
    if (fclose(dummy))
      exit(error("", ERRCLOSE));
    in->gzf = gzopen(inFile, "r");
    if (in->gzf == NULL)
      exit(error(inFile, ERROPEN));
  } else {
    if (! stdinBool)
      rewind(dummy);
    in->f = dummy;
  }

  return gzip;
}

/* void runProgram()
 * Controls the opening/closing of files,
 *   and analysis by readFile().
 */
void runProgram(char* outFile, char* inFile1,
    char* inFile2, bool inter, char* unFile,
    char* logFile, int overlap, bool dovetail,
    char* doveFile, int doveOverlap, char* alnFile,
    int alnOpt, bool adaptOpt, int gzOut, bool fjoin,
    bool interOpt, float mismatch, bool maxLen,
    int offset, int maxQual, bool verbose,
    int threads) {

  // get first set of input file names
  char* end1, *end2;
  char* file1 = strtok_r(inFile1, COM, &end1);
  char* file2 = file1;
  if (! inter)
    file2 = strtok_r(inFile2, COM, &end2);

  // loop through input files
  File out, out2, un1, un2, log, dove, aln; // output files
  int i = 0;  // number of files processed
  int tCount = 0, tStitch = 0;  // counting variables
  while (file1 && file2) {

    // open files
    File in1, in2;
    bool gz1 = openRead(file1, &in1);
    bool gz2 = gz1;
    if (! inter)
      gz2 = openRead(file2, &in2);

    // on first iteration, open output files
    if (! i) {
      if (gzOut == -1)
        gzOut = 0;
      else if (gz1 || gz2)
        gzOut = 1;
      openFiles(outFile, &out, &out2,
        unFile, &un1, &un2, logFile, &log,
        doveFile, &dove, dovetail, alnFile, &aln,
        adaptOpt, gzOut, interOpt);
    }

    // process files
    if (verbose)
      fprintf(stderr, "Processing files: %s,%s\n", file1,
        inter ? "(interleaved)" : file2);
    int stitch = 0;  // counting variable
    int count = readFile(in1, inter ? in1 : in2,
      out, interOpt ? out : out2,
      un1, interOpt ? un1 : un2, unFile != NULL,
      log, logFile != NULL,
      overlap, dovetail, doveOverlap, dove,
      dovetail && doveFile != NULL, aln, alnOpt,
      adaptOpt, mismatch, maxLen, &stitch,
      offset, maxQual, gz1, gz2, gzOut, fjoin, threads);
    tCount += count;
    tStitch += stitch;

    // log counts
    if (verbose) {
      fprintf(stderr, "  Fragments (pairs of reads) analyzed: %d\n", count);
      if (adaptOpt)
        fprintf(stderr, "  Adapters removed: %d\n", stitch);
      else
        fprintf(stderr, "  Successfully stitched: %d\n", stitch);
    }

    // close input files
    if ( (gz1 && gzclose(in1.gzf) != Z_OK) || (! gz1 && fclose(in1.f))
        || (! inter && ( (gz2 && gzclose(in2.gzf) != Z_OK)
        || (! gz2 && fclose(in2.f)) ) ) )
      exit(error("", ERRCLOSE));

    file1 = strtok_r(NULL, COM, &end1);
    file2 = file1;
    if (! inter)
      file2 = strtok_r(NULL, COM, &end2);
    i++;
  }

  if (verbose && i > 1) {
    fprintf(stderr, "Total counts\n");
    fprintf(stderr, "  Fragments (pairs of reads) analyzed: %d\n", tCount);
    if (adaptOpt)
      fprintf(stderr, "  Adapters removed: %d\n", tStitch);
    else
      fprintf(stderr, "  Successfully stitched: %d\n", tStitch);
  }

  // close files
  if ( ( gzOut && ( gzclose(out.gzf) != Z_OK ||
      (adaptOpt && ! interOpt && gzclose(out2.gzf) != Z_OK) ||
      (unFile != NULL && (gzclose(un1.gzf) != Z_OK ||
      (! interOpt && gzclose(un2.gzf) != Z_OK)) ) ) ) ||
      ( ! gzOut && ( fclose(out.f) ||
      (adaptOpt && ! interOpt && fclose(out2.f)) ||
      (unFile != NULL && (fclose(un1.f) ||
      (! interOpt && fclose(un2.f)) ) ) ) ) ||
      (logFile != NULL && fclose(log.f)) ||
      (dovetail && doveFile != NULL && fclose(dove.f)) ||
      (alnFile != NULL && fclose(aln.f)) )
    exit(error("", ERRCLOSE));
}

/* void getArgs()
 * Parse the command-line. Check for errors.
 */
void getArgs(int argc, char** argv) {

  // default parameters/filenames
  char* outFile = NULL, *inFile1 = NULL, *inFile2 = NULL,
    *unFile = NULL, *logFile = NULL, *doveFile = NULL,
    *alnFile = NULL;
  int overlap = DEFOVER, doveOverlap = DEFDOVE, gzOut = 0,
    offset = OFFSET, maxQual = MAXQUAL, threads = DEFTHR;
  float mismatch = DEFMISM;
  bool dovetail = false, adaptOpt = false, maxLen = true,
    diffOpt = false, interOpt = false, fjoin = false,
    verbose = false;

  // parse argv
  int c;
  while ( (c = getopt_long(argc, argv, OPTIONS, long_options, NULL)) != -1 )
    switch (c) {
      case HELP: usage(); break;
      case VERSOPT: printVersion(); break;
      case MAXOPT: maxLen = false; break;
      case DOVEOPT: dovetail = true; break;
      case ADAPTOPT: adaptOpt = true; break;
      case GZOPT: gzOut = 1; break;
      case UNGZOPT: gzOut = -1; break;
      case DIFFOPT: diffOpt = true; break;
      case INTEROPT: interOpt = true; break;
      case FJOINOPT: fjoin = true; break;
      case VERBOSE: verbose = true; break;
      case OUTFILE: outFile = optarg; break;
      case FIRST: inFile1 = optarg; break;
      case SECOND: inFile2 = optarg; break;
      case UNFILE: unFile = optarg; break;
      case LOGFILE: logFile = optarg; break;
      case DOVEFILE: doveFile = optarg; break;
      case ALNFILE: alnFile = optarg; break;
      case OVERLAP: overlap = getInt(optarg); break;
      case DOVEOVER: doveOverlap = getInt(optarg); break;
      case MISMATCH: mismatch = getFloat(optarg); break;
      case QUALITY: offset = getInt(optarg); break;
      case SETQUAL: maxQual = getInt(optarg); break;
      case THREADS: threads = getInt(optarg); break;
      default: exit(-1);
    }
  if (optind < argc)
    exit(error(argv[optind], ERRPARAM));

  // check for argument errors
  if (outFile == NULL || inFile1 == NULL) {
    error("", ERRFILE);
    usage();
  }
  bool inter = false;  // interleaved input
  if (inFile2 == NULL) {
    if (verbose)
      fprintf(stderr, "Warning: only one input file specified -- assuming interleaved\n");
    inter = true;
  }
  if (overlap <= 0 || doveOverlap <= 0)
    exit(error("", ERROVER));
  if (mismatch < 0.0f || mismatch >= 1.0f)
    exit(error("", ERRMISM));
  if (threads < 1)
    exit(error("", ERRTHREAD));

  // adjust parameters for adapter-removal mode
  if (adaptOpt) {
    dovetail = true;
    unFile = logFile = alnFile = NULL;
  }
  int alnOpt = (alnFile != NULL ? (diffOpt ? 2 : 1) : 0);

  // send arguments to runProgram()
  runProgram(outFile, inFile1, inFile2, inter, unFile,
    logFile, overlap, dovetail, doveFile, doveOverlap,
    alnFile, alnOpt, adaptOpt, gzOut, fjoin, interOpt,
    mismatch, maxLen, offset, maxQual, verbose, threads);
}

/* int main()
 * Main.
 */
int main(int argc, char* argv[]) {
  getArgs(argc, argv);
  return 0;
}
