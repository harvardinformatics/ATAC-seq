#!/usr/bin/env python

# John M. Gaspar (jsh58@wildcats.unh.edu)
# Oct. 2017

# Remove reference sequence(s) (e.g. chrM) from a SAM.
#   Adjust FLAG and other fields so that perhaps even
#   Picard won't object (unlikely, but...)

import sys
import gzip
import re

def openRead(filename):
  '''
  Open filename for reading. '-' indicates stdin.
    '.gz' suffix indicates gzip compression.
  '''
  if filename == '-':
    return sys.stdin
  try:
    if filename[-3:] == '.gz':
      f = gzip.open(filename, 'rb')
    else:
      f = open(filename, 'rU')
  except IOError:
    sys.stderr.write('Error! Cannot open %s for reading\n' % filename)
    sys.exit(-1)
  return f

def openWrite(filename):
  '''
  Open filename for writing. '-' indicates stdout.
    '.gz' suffix indicates gzip compression.
  '''
  if filename == '-':
    return sys.stdout
  try:
    if filename[-3:] == '.gz':
      f = gzip.open(filename, 'wb')
    else:
      f = open(filename, 'w')
  except IOError:
    sys.stderr.write('Error! Cannot open %s for writing\n' % filename)
    sys.exit(-1)
  return f

def procHeader(line, chrom):
  '''Check for removed chroms in header line.'''
  mat = re.search(r'@SQ\s+SN:(\S+)', line)
  if mat and mat.group(1) in chrom:
    return False
  return True

def unmapRead(fOut, spl, res):
  '''Adjust alignments of a read.'''
  flag = int(spl[1])
  ret = 0  # return value: 1 if alignment removed

  if res[0]:
    # unmap read
    if flag & 0x900:
      return ret   # skip if secondary/supplementary
    flag |= 0x4    # unmapped
    flag &= 0xFFD  # not properly paired
    spl[2] = '*'
    spl[3] = '0'
    spl[4] = '0'  # should be '255' for unavailable...
    spl[5] = '*'
    spl[8] = '0'
    del spl[11:]  # remove optional fields, if any
    spl[10] = spl[10].rstrip() + '\n'
    ret = 1

  if res[1]:
    # unmap mate
    flag |= 0x8    # mate unmapped
    flag &= 0xFFD  # not properly paired
    spl[6] = '*'
    spl[7] = '0'
    spl[8] = '0'

  # write result
  spl[1] = str(flag)
  fOut.write('\t'.join(spl))
  return ret

def procFiles(fIn, fOut, chrom):
  '''Process SAM files.'''
  total = unmap = 0
  for line in fIn:
    if line[0] == '@':
      if procHeader(line, chrom):
        fOut.write(line)
      continue

    # determine if refNames are on blacklist
    res = [False, False]
    spl = line.split('\t')
    if len(spl) < 11:
      sys.stderr.write('Error! Poorly formatted SAM record:\n' \
        + line)
      sys.exit(-1)
    if spl[2] in chrom:
      res[0] = True
    if spl[6] in chrom or (spl[6] == '=' and spl[2] in chrom):
      res[1] = True

    unmap += unmapRead(fOut, spl, res)
    total += 1

  sys.stderr.write('Alignments analyzed: %12d\n' % total)
  sys.stderr.write('   Removed from SAM: %12d\n' % unmap)

def main():
  '''Main.'''
  args = sys.argv[1:]
  if len(args) < 3:
    sys.stderr.write('Usage: removeChrom <inSAM> <outSAM> [<chrom>]+\n')
    sys.stderr.write('  [<chrom>]+  One or more chromosomes (reference\n')
    sys.stderr.write('              sequences) to be removed from the SAM\n')
    sys.exit(-1)

  fIn = openRead(args[0])
  fOut = openWrite(args[1])
  chrom = args[2:]

  procFiles(fIn, fOut, chrom)

  if fIn != sys.stdin:
    fIn.close()
  if fOut != sys.stdout:
    fOut.close()

if __name__ == '__main__':
  main()
