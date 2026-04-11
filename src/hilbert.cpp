#include <cassert>
#include <cstdint>
#include <vector>

#ifdef STANDALONE_BUILD
  #include "hilbert.hpp"
#else
  #include "../inst/include/hilbert.hpp"
#endif

namespace hilbert {

namespace detail {

// Transform coordinate array X[0..d-1] in-place to the transposed
// Hilbert representation (axes → index).
static void axesToTranspose(uint32_t* X, int bitsPerDim, int d) {
  const uint32_t M = 1U << (bitsPerDim - 1);
  for (uint32_t Q = M; Q > 1U; Q >>= 1U) {
    const uint32_t P = Q - 1U;
    for (int i = 0; i < d; i++) {
      if ((X[i] & Q) != 0U) {
        X[0] ^= P; // invert
      } else {
        const uint32_t t = (X[0] ^ X[i]) & P;
        X[0] ^= t;
        X[i] ^= t; // swap
      }
    }
  }
  // Gray encode
  for (int i = 1; i < d; i++) {
    X[i] ^= X[i - 1];
  }
  uint32_t t = 0U;
  for (uint32_t Q = M; Q > 1U; Q >>= 1U) {
    if ((X[d - 1] & Q) != 0U) {
      t ^= Q - 1U;
    }
  }
  for (int i = 0; i < d; i++) {
    X[i] ^= t;
  }
}

// Inverse: transposed Hilbert representation → coordinate array.
static void transposeToAxes(uint32_t* X, int bitsPerDim, int d) {
  const uint32_t N = 2U << (bitsPerDim - 1); // 2^bitsPerDim
  // Gray decode
  const uint32_t t0 = X[d - 1] >> 1U;
  for (int i = d - 1; i > 0; i--) {
    X[i] ^= X[i - 1];
  }
  X[0] ^= t0;
  // Undo excess work
  for (uint32_t Q = 2U; Q != N; Q <<= 1U) {
    const uint32_t P = Q - 1U;
    for (int i = d - 1; i >= 0; i--) {
      if ((X[i] & Q) != 0U) {
        X[0] ^= P;
      } else {
        const uint32_t t = (X[0] ^ X[i]) & P;
        X[0] ^= t;
        X[i] ^= t;
      }
    }
  }
}

// Pack transposed X[0..d-1] into a uint64_t Hilbert index.
// Bit k of X[i] maps to bit (d*k + (d-1-i)) of H (MSB-first interleaving).
static uint64_t transposeToIndex(const uint32_t* X, int bitsPerDim, int d) {
  uint64_t h = 0;
  for (int i = 0; i < d; i++) {
    for (int k = 0; k < bitsPerDim; k++) {
      if ((X[i] & (1U << k)) != 0U) {
        h |= 1ULL << (d * k + (d - 1 - i));
      }
    }
  }
  return h;
}

// Unpack a uint64_t Hilbert index into transposed X[0..d-1].
static void indexToTranspose(uint64_t h, uint32_t* X, int bitsPerDim, int d) {
  for (int i = 0; i < d; i++) {
    X[i] = 0U;
    for (int k = 0; k < bitsPerDim; k++) {
      if ((h & (1ULL << (d * k + (d - 1 - i)))) != 0U) {
        X[i] |= 1U << k;
      }
    }
  }
}

} // namespace detail

int bitsFor(int maxVal) {
  assert(maxVal >= 0);
  int b = 1;
  while ((1 << b) <= maxVal) {
    ++b;
  }
  return b;
}

uint64_t encode(const int* coords, int d, int bitsPerDim) {
  assert(d >= 1);
  assert(bitsPerDim >= 1);
  assert(d * bitsPerDim <= 64);
  std::vector<uint32_t> X(d);
  for (int i = 0; i < d; i++) {
    assert(coords[i] >= 0 && coords[i] < (1 << bitsPerDim));
    X[i] = static_cast<uint32_t>(coords[i]);
  }
  detail::axesToTranspose(X.data(), bitsPerDim, d);
  return detail::transposeToIndex(X.data(), bitsPerDim, d);
}

uint64_t encode(const std::vector<int>& coords, int bitsPerDim) {
  return encode(coords.data(), static_cast<int>(coords.size()), bitsPerDim);
}

std::vector<int> decode(uint64_t h, int d, int bitsPerDim) {
  assert(d >= 1);
  assert(bitsPerDim >= 1);
  assert(d * bitsPerDim <= 64);
  std::vector<uint32_t> X(d);
  detail::indexToTranspose(h, X.data(), bitsPerDim, d);
  detail::transposeToAxes(X.data(), bitsPerDim, d);
  std::vector<int> coords(d);
  for (int i = 0; i < d; i++) {
    coords[i] = static_cast<int>(X[i]);
  }
  return coords;
}

} // namespace hilbert
