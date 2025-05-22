# vpack
A single header library for compressing and decompressing SNV data using a 64-bit bit-packing scheme

## Quick Start
1. Clone this repo
2. Include `vpack.h`

## Overview
vpack uses a bit-packing scheme to compress SNV data into 64-bit integers, achieving ~2.5x compression
for SNVs in single-sample VCFs.

There are structures and methods for compressing/decompressing single variants and genotypes into integers, as many genotypes.


### One SNP, one genotype
In some cases, you may want to compress a single diploid genotype with the variant information. The function `snvpack64` compresses the core information regarding a single genotype and a single variant into a 64-bit integer. To unpack, use `snvunpack`.

Most of the time, you have many genotypes and many variants. See [the API below.](#Many-variants-and-many-genotypes)

The bit-packing scheme is detailed below:
```
  sample ID       chrom      position               ref alt  gt
00000000000000000|00000|0000000000000000000000000000|00|00|000000000
  9 bits for Genotype information: {0, 1, /, |, .}
  2 bit encoding for Ref and Alt alleles {A, T, G, C}
  28 bit encoding for position (0...268Mb)
  5 bit chromosome encoding {0...31}
  17 bits remaining for numeric sample index or ID {0...131k}
```
Example usage:
```C
  int pos = 12345644;
  int sample_idx = 105045;
  int chrom = 22;
  char ref = 'A';
  char alt = 'G';
  char gts[3] = "G/A";
  vpack64_t v = snvpack64(sample_idx, chrom, pos, ref, alt, gts);
```

### Many variants and many genotypes
`vpack` provides an API for the compression and decompression of large arrays of variants and genotypes. The variant data is packed into a 64-bit integer with a similar scheme as above, and the diploid genotypes are sequentially packed into 64-bit integers (32 diploid genotypes per integer).

The bitpacking scheme is detailed below:
```
SNV bit packing without GT
  sample ID                chrom      position               ref alt
00000000000000000000000000|00000|0000000000000000000000000000|00|00
  2 bit encoding for Ref and Alt alleles {A, T, G, C}
  28 bit encoding for position (0...268Mb)
  5 bit chromosome encoding {0...31}
  26 bits remaining for numeric sample index or ID {}
```

Example packing and unpacking genotypes:
```C
  vpack_rec1_t rec = vpack_rec1_t_init(); // Zeroed initialization

  // Pack 2 samples
  vpack_rec(&rec, 'C', 'G'); // Sample_idx = 0
  vpack_rec(&rec, 'C', 'C'); // Sample_idx = 1

  // Unpack
  vunpack_rec(&rec);
  // Recall the genotype for sample_idx = 0
  vpack_gt_t genotype;
  vp_get_genotype(&rec, 0, &genotype);
  printf("%c/%c\n", genotype.a, genotype.b);
```
Example unpacking genotypes:
```C

vpack64_t v = 101;
int nsamples = 2;
vpack_rec1_t rec = vpack_rec1_load(v, nsamples);

vunpack_rec(&rec);

vpack_gt_t genotype;
vp_get_genotype(&rec, 0, &genotype);
printf("%c/%c\n", genotype.a, genotype.b)

```




