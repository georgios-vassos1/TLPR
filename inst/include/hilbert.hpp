#ifndef HILBERT_HPP
#define HILBERT_HPP

#include <cstdint>
#include <vector>

// Hilbert curve encoding / decoding for d-dimensional integer hypercubes.
//
// Algorithm: Skilling, J. (2004) "Programming the Hilbert curve",
//            AIP Conf. Proc. 707, pp. 381-387.
//
// Representation
// --------------
// Coordinates: each dimension i has a value in [0, 2^bitsPerDim - 1].
// Index:       a single uint64_t encoding all dimensions via MSB-first
//              bit interleaving — bit k of coordinate i maps to bit
//              (d*k + (d-1-i)) of the index.
//
// Constraints
// -----------
//   d * bitsPerDim <= 64   (index must fit in uint64_t)
//   bitsPerDim >= 1
//
// Typical usage in TLPR
// ----------------------
//   int b = hilbert::bitsFor(R);           // e.g. R=40 → b=6
//   uint64_t h = hilbert::encode(coords, b);
//   uint64_t bucket = hilbert::coarsen(h, shiftBits);

namespace hilbert {

// Return the number of bits needed to represent values in [0, maxVal].
// bitsFor(0) = bitsFor(1) = 1, bitsFor(2) = bitsFor(3) = 2, etc.
int bitsFor(int maxVal);

// Encode d-dimensional coordinates to a Hilbert index.
// coords[i] must be in [0, 2^bitsPerDim - 1].
// Requires d * bitsPerDim <= 64.
uint64_t encode(const int* coords, int d, int bitsPerDim);

// Convenience overload accepting std::vector<int>.
uint64_t encode(const std::vector<int>& coords, int bitsPerDim);

// Decode a Hilbert index to d-dimensional coordinates.
std::vector<int> decode(uint64_t h, int d, int bitsPerDim);

// Coarsen a Hilbert index by dropping the shiftBits least-significant bits,
// assigning the encoded point to an aggregation bucket.
inline uint64_t coarsen(uint64_t index, int shiftBits) {
  return index >> shiftBits;
}

} // namespace hilbert

#endif // HILBERT_HPP
