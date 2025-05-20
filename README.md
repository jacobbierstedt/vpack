# vpack
A single header library for compressing and decompressing SNV data using a 64-bit bit-packing scheme


## Overview
vpack uses a bit-packing scheme to compress SNV data into 64-bit integers, achieving ~2.5x compression
for SNVs in single-sample VCFs.

```
  sample ID       chrom      position               ref alt  gt
00000000000000000|00000|0000000000000000000000000000|00|00|000000000
  9 bits for Genotype information: {0, 1, /, |, .}
  2 bit encoding for Ref and Alt alleles {A, T, G, C}
  28 bit encoding for position (0...268Mb)
  5 bit chromosome encoding {0...31}
  17 bits remaining for numeric sample index or ID {0...131k}
```

## Quick Start
1. Clone this repo
2. Include `vpack.h`

