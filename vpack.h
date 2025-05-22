/*
  vpack - A single header library for compression and decompression
  of genomic variant data

  Copyright (C) 2025 Jacob Bierstedt

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program. If not, see <https://www.gnu.org/licenses/>.
*/
#ifndef VPACK_H
#define VPACK_H

#include <stdint.h>
#include <stdio.h>



/*
Single SNV bit packing strategy
  sample ID       chrom      position               ref alt  gt
00000000000000000|00000|0000000000000000000000000000|00|00|000000000

  9 bits for Genotype information: {0, 1, /, |, .}
  2 bit encoding for Ref and Alt alleles {A, T, G, C}
  28 bit encoding for position (0...268Mb)
  5 bit chromosome encoding {0...31}
  17 bits remaining for numeric sample index or ID {0...131k}

SNV bit packing without GT
  sample ID                chrom      position               ref alt
00000000000000000000000000|00000|0000000000000000000000000000|00|00
  2 bit encoding for Ref and Alt alleles {A, T, G, C}
  28 bit encoding for position (0...268Mb)
  5 bit chromosome encoding {0...31}
  26 bits remaining for numeric sample index or ID {}
*/

typedef uint64_t vpack64_t;
typedef uint32_t vpidx_t;

#define _DNA_8_MASK 0x07
#define _DECODE_8_MASK 0x03
#define VMASK_28 0x0FFFFFFF
#define VMASK_5  0x01F

#ifndef DNA_8
#define DNA_8
                                    /*  A     C  T        G*/
static const uint8_t dna_8[8]     = {4, 0, 4, 1, 3, 4, 4, 2};
static const uint8_t dec_dna_8[8] = {'A', 'C', 'G', 'T', 'N'};

#define ENCODE(c)    dna_8[c & _DNA_8_MASK]
#endif

#define VPACK_BASE(v, base)    (*v) <<= 2; (*v) |= ENCODE(base);
#define VUNPACK_BASE(v, base)  (*base) = dec_dna_8[(*v) & _DECODE_8_MASK]; (*v) >>= 2;


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*             SINGLE VARIANT PACKING
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*
  Bit pack formatted GT into 9 bits, with 3 bits representing 
  each original value
*/
static inline vpack64_t vpack_gt9(vpack64_t* v, uint8_t* gt) {
  for (int i=0; i<3; i++) {
    (*v) <<= 3;
    switch (gt[i]) {
      case '0': (*v) |= 0; break;
      case '1': (*v) |= 1; break;
      case '.': (*v) |= 2; break;
      case '/': (*v) |= 3; break;
      case '|': (*v) |= 4; break;
      default:             break;
    }
  }
}

/*
  @brief
  Unpack 9-bit GT to original values
*/
static inline void vunpack_gt9(vpack64_t* v, uint8_t* u) {
  for (int i=2; i>=0; i--) {
    switch ((*v) & 7)
    {
    case 0: u[i] = '0'; break;
    case 1: u[i] = '1'; break;
    case 2: u[i] = '.'; break;
    case 3: u[i] = '/'; break;
    case 4: u[i] = '|'; break;
    default:            break;
    }
    (*v) >>= 3;
  }

}

/*
  @brief
  Pack SNV variant information into a 64-bit integer
  @param sample_idx (int) sample index or ID
  @param chrom      (int) chromosome (0-32)
  @param pos        (int) position
  @param ref        (char) ref allele (len = 1)
  @param ref        (char) alt allele (len = 1)
  @param gt         (char*) GT in array, len>=3 (i.e. '0/1')
*/
static inline vpack64_t snvpack64(uint32_t sample_idx, uint32_t chrom, uint32_t pos, uint8_t ref, uint8_t alt, uint8_t* gt) {
  vpack64_t v = 0x0;
  v |= sample_idx;
  v <<= 5;
  v |= chrom;
  v <<= 28;
  v |= pos;
  VPACK_BASE(&v, ref);
  VPACK_BASE(&v, alt);
  vpack_gt9(&v, gt);
  return v;
}
/*
  @brief
  Unpack 64-bit integer to original SNV variant information
  @param v          (vpack64_t) 64-bit packed integer
  @param sample_idx (int) sample index or ID
  @param chrom      (int) chromosome (0-32)
  @param pos        (int) position
  @param ref        (char) ref allele (len = 1)
  @param ref        (char) alt allele (len = 1)
  @param gt         (char*) GT in array, len>=3 (i.e. '0/1')
*/
static inline void snvunpack64(vpack64_t* v, uint32_t* sample_idx, uint32_t* chrom, uint32_t* pos, uint8_t* ref, uint8_t* alt, uint8_t* gt) {
  vunpack_gt9(v, gt);
  VUNPACK_BASE(v, alt);
  VUNPACK_BASE(v, ref);
  (*pos) = (*v) & VMASK_28;
  (*v) >>= 28;
  (*chrom) = (*v) & VMASK_5;
  (*v) >>= 5;
  (*sample_idx) = (*v);
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*             MULTI-SAMPLE PACKING
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/*
  @brief
  Pack variant information into a 64-bit integer

  @params
  @param chrom chromosome
  @param pos   position
  @param ref   reference allele (A,T,C,G)
  @param alt   alt allele (A,T,C,G)
*/
static inline vpack64_t vpack64_loc(uint32_t chrom, uint32_t pos, uint8_t ref, uint8_t alt)
{
  vpack64_t v = 0x0;
  v |= chrom;
  v <<= 28;
  v |= pos;
  VPACK_BASE(&v, ref);
  VPACK_BASE(&v, alt);
  return v;
}
/*
  @brief
  Data structure for packing up to 32 diploid genotypes into a 
  64-bit integer.
  
  Must initialize to zero with `vpack_rec1_t_init()`
*/
typedef struct
{
  vpack64_t v;
  uint32_t offset;
  uint8_t gt[64];

  uint32_t u; // 0: initialized, 1: packed, 2: unpacked
} vpack_rec1_t;
/*
  @brief
  Initializer function for vpack_reck1_t, to ensure zeroed initialization
*/
static inline vpack_rec1_t vpack_rec1_t_init() {
    return (vpack_rec1_t){0};
}
/*
  @brief
  Data structure for a diploid genotype. `a`->`b` corresponds to the
  order in which the alleles were packed. See `vpack_rec` for details.
*/
typedef struct
{
  uint8_t a;
  uint8_t b;
} vpack_gt_t;
/*
  @brief
  Pack a single genotype into a record. Alleles are packed `a` then 
  `b`. The order of packed records should be kept track of. As these
  will be needed to reference the correct genotype after unpacking.
  
  Returns a status code 0 on success, -1 when out of bounds.

  @returns status  0: success, -1: error, out of bounds
*/
static inline int vpack_rec(vpack_rec1_t* rec, uint8_t a, uint8_t b)
{
  if (rec->offset >= 63) return -1;
  VPACK_BASE(&rec->v, a);
  rec->offset++; 
  VPACK_BASE(&rec->v, b);
  rec->offset++; 
  rec->u = 1;
  return 0;
}
/*
  @brief
  Unpack a packed record. See `vpack_rec` for info
  on packed records.

  @param rec vpack_rec1_t
*/
static inline void vunpack_rec(vpack_rec1_t* rec)
{
  for (int i = rec->offset - 1;i>= 0; --i) {
    rec->gt[i] = dec_dna_8[rec->v & 3];
    rec->v >>= 2;
  }
  rec->u = 2;
}
/*
  Get genotype from a packed record. Record must be packed and
  unpacked before use. See `vpack_rec` and `vunpack_rec`.

  @param
  @param rec           record
  @param sample_idx    sample index in the record
  @param gt            vpack_gt_t to hold genotype information for sample

  @returns
  @return status  0: success, -1: error, sample_idx out of bounds, record hasn't been unpacked
*/
static inline int vp_get_genotype(vpack_rec1_t* rec, uint32_t sample_idx, vpack_gt_t* gt)
{
  if(rec->u<2 || sample_idx >= (int)rec->offset/2) return -1;
  gt->a = rec->gt[sample_idx*2];
  gt->b = rec->gt[sample_idx*2 + 1];
  return 0;
}

#endif /* VPACK_H */