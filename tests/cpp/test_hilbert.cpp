#include "/Users/georgios.vassos1/drayage/TLPR/inst/include/hilbert.hpp"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <vector>

static int failures = 0;
#define CHECK(cond, msg) \
  do { if (!(cond)) { printf("  FAIL: %s\n", msg); ++failures; } } while(0)

// ── reference implementation (Wikipedia d2xy / xy2d, 2D only) ─────────────────
// Used solely as ground truth; independent of Skilling's algorithm.
static void ref_rot(int n, int* x, int* y, int rx, int ry) {
  if (ry == 0) {
    if (rx == 1) { *x = n-1 - *x; *y = n-1 - *y; }
    int t = *x; *x = *y; *y = t;
  }
}
static int ref_xy2d(int n, int x, int y) {  // n must be power of 2
  int d = 0;
  for (int s = n/2; s > 0; s /= 2) {
    int rx = (x & s) > 0;
    int ry = (y & s) > 0;
    d += s * s * ((3 * rx) ^ ry);
    ref_rot(s, &x, &y, rx, ry);
  }
  return d;
}
static void ref_d2xy(int n, int d, int* x, int* y) {
  int t = d; *x = *y = 0;
  for (int s = 1; s < n; s *= 2) {
    int rx = 1 & (t / 2);
    int ry = 1 & (t ^ rx);
    ref_rot(s, x, y, rx, ry);
    *x += s * rx; *y += s * ry;
    t /= 4;
  }
}

// ── helpers ───────────────────────────────────────────────────────────────────
static int l1(const std::vector<int>& a, const std::vector<int>& b) {
  int s = 0;
  for (int i = 0; i < (int)a.size(); i++) s += std::abs(a[i]-b[i]);
  return s;
}
static int linf(const std::vector<int>& a, const std::vector<int>& b) {
  int m = 0;
  for (int i = 0; i < (int)a.size(); i++) m = std::max(m, std::abs(a[i]-b[i]));
  return m;
}

// ── 1. bitsFor ────────────────────────────────────────────────────────────────
void test_bitsFor() {
  CHECK(hilbert::bitsFor(0)  == 1, "bitsFor(0)");
  CHECK(hilbert::bitsFor(1)  == 1, "bitsFor(1)");
  CHECK(hilbert::bitsFor(2)  == 2, "bitsFor(2)");
  CHECK(hilbert::bitsFor(3)  == 2, "bitsFor(3)");
  CHECK(hilbert::bitsFor(4)  == 3, "bitsFor(4)");
  CHECK(hilbert::bitsFor(10) == 4, "bitsFor(10)");
  CHECK(hilbert::bitsFor(15) == 4, "bitsFor(15)");
  CHECK(hilbert::bitsFor(16) == 5, "bitsFor(16)");
  CHECK(hilbert::bitsFor(40) == 6, "bitsFor(40)");
  CHECK(hilbert::bitsFor(63) == 6, "bitsFor(63)");
  CHECK(hilbert::bitsFor(64) == 7, "bitsFor(64)");
  printf("bitsFor: OK\n");
}

// ── 2. bijectivity ────────────────────────────────────────────────────────────
// Verifies encode is a bijection onto [0, 2^(d*b)-1].
// Implies round-trip correctness (no need for a separate encode→decode test).
void test_bijectivity_2d() {
  int b = 4, total = (1<<b)*(1<<b);
  std::vector<uint64_t> hs;
  hs.reserve(total);
  for (int x = 0; x < (1<<b); x++)
    for (int y = 0; y < (1<<b); y++)
      hs.push_back(hilbert::encode({x, y}, b));
  std::sort(hs.begin(), hs.end());
  bool unique = std::adjacent_find(hs.begin(), hs.end()) == hs.end();
  CHECK(unique, "bijectivity 2D: no duplicate indices");
  CHECK(hs.front()==0 && hs.back()==(uint64_t)(total-1),
        "bijectivity 2D: index range is exactly [0, N-1]");
  printf("bijectivity 2D (b=%d, %d points)\n", b, total);
}

void test_bijectivity_3d() {
  int b = 3, d = 3, total = 1<<(d*b);
  std::vector<uint64_t> hs;
  hs.reserve(total);
  for (int x = 0; x < (1<<b); x++)
    for (int y = 0; y < (1<<b); y++)
      for (int z = 0; z < (1<<b); z++)
        hs.push_back(hilbert::encode({x,y,z}, b));
  std::sort(hs.begin(), hs.end());
  bool unique = std::adjacent_find(hs.begin(), hs.end()) == hs.end();
  CHECK(unique, "bijectivity 3D: no duplicate indices");
  CHECK(hs.front()==0 && hs.back()==(uint64_t)(total-1),
        "bijectivity 3D: index range is exactly [0, N-1]");
  printf("bijectivity 3D (b=%d, %d points)\n", b, total);
}

// ── 3. unit-step / face-adjacency (L1 = 1) ───────────────────────────────────
// The defining Hilbert property: consecutive curve indices are face-adjacent —
// exactly one coordinate changes, by exactly 1.  L1=1 is strictly stronger
// than L∞=1 (which allows diagonal steps).  This is the test from
// adishavit/hilbert ("miscount == 1") that caught our bit-packing inversion.
void test_unit_step_2d() {
  int b = 5, side = 1<<b;
  int max_l1 = 0;
  for (uint64_t h = 0; h+1 < (uint64_t)(side*side); h++) {
    int d = l1(hilbert::decode(h, 2, b), hilbert::decode(h+1, 2, b));
    max_l1 = std::max(max_l1, d);
  }
  CHECK(max_l1 == 1, "unit-step 2D: consecutive indices are face-adjacent (L1=1)");
  printf("unit-step 2D (b=%d, %d steps): max L1 between consecutive = %d\n",
         b, side*side-1, max_l1);
}

void test_unit_step_3d() {
  int b = 3, side = 1<<b;
  int max_l1 = 0;
  for (uint64_t h = 0; h+1 < (uint64_t)(side*side*side); h++) {
    int d = l1(hilbert::decode(h, 3, b), hilbert::decode(h+1, 3, b));
    max_l1 = std::max(max_l1, d);
  }
  CHECK(max_l1 == 1, "unit-step 3D: consecutive indices are face-adjacent (L1=1)");
  printf("unit-step 3D (b=%d, %d steps): max L1 = %d\n",
         b, side*side*side-1, max_l1);
}

// ── 4. ground-truth cross-check against reference 2D implementation ───────────
// The only test that can detect a self-consistently-wrong implementation.
// Uses the Wikipedia d2xy / xy2d algorithm (independently verified, not Skilling)
// as an oracle.  Self-consistency (encode→decode→encode) cannot catch this.
void test_ground_truth_2d() {
  // Spot-checks at specific (b, x, y, expected_h) tuples derived from
  // the Wikipedia reference independently.
  struct Case { int b, x, y; uint64_t h; };
  // Generated by running ref_xy2d and ref_d2xy at b=1..4:
  // b=1 (2x2):  (0,0)→0  (0,1)→1  (1,1)→2  (1,0)→3
  // b=2 (4x4):  (0,0)→0  (1,0)→1  (1,1)→2  (0,1)→3  (0,2)→4 ...
  // We generate these at runtime from the reference, then compare our encode.
  int fail = 0;
  for (int b = 1; b <= 5; b++) {
    int side = 1 << b;
    int total = side * side;
    for (int ref_d = 0; ref_d < total; ref_d++) {
      int rx, ry;
      ref_d2xy(side, ref_d, &rx, &ry);
      uint64_t our_h = hilbert::encode({rx, ry}, b);
      if (our_h != (uint64_t)ref_d) {
        if (fail < 3)
          printf("  mismatch b=%d: ref(%d,%d)→%d  ours→%llu\n",
                 b, rx, ry, ref_d, (unsigned long long)our_h);
        fail++;
      }
      // Also verify our decode matches reference
      auto our_pt = hilbert::decode((uint64_t)ref_d, 2, b);
      if (our_pt[0] != rx || our_pt[1] != ry) {
        if (fail < 3)
          printf("  decode mismatch b=%d d=%d: ref(%d,%d) ours(%d,%d)\n",
                 b, ref_d, rx, ry, our_pt[0], our_pt[1]);
        fail++;
      }
    }
  }
  CHECK(fail == 0, "ground-truth 2D: matches Wikipedia reference for b=1..5");
  printf("ground-truth 2D (b=1..5, %d total points, %d failures)\n",
         (1<<1)*(1<<1) + (1<<2)*(1<<2) + (1<<3)*(1<<3) +
         (1<<4)*(1<<4) + (1<<5)*(1<<5), fail);
}

// ── 5. coarsening consistency ─────────────────────────────────────────────────
// The property the TLPR state aggregation relies on: all points sharing a
// Hilbert bucket are within a bounded grid distance of each other.  The bound
// is 2^shiftBits - 1 (each curve step is L∞=1, and a bucket spans 2^s steps).
// Not tested by any of the surveyed libraries — specific to our aggregation use.
void test_coarsening_consistency() {
  int b = 4, side = 1<<b, total = side*side;

  for (int shiftBits = 1; shiftBits <= 4; shiftBits++) {
    int bound = (1 << shiftBits) - 1;
    int actual_max = 0;

    std::vector<std::pair<uint64_t, std::vector<int>>> bp;
    bp.reserve(total);
    for (int x = 0; x < side; x++)
      for (int y = 0; y < side; y++)
        bp.push_back({hilbert::coarsen(hilbert::encode({x,y}, b), shiftBits), {x,y}});
    std::sort(bp.begin(), bp.end(), [](const auto& a, const auto& b){ return a.first < b.first; });

    // Walk through each contiguous group (same bucket) and measure diameter
    int i = 0;
    while (i < total) {
      int j = i;
      while (j < total && bp[j].first == bp[i].first) j++;
      for (int a = i; a < j; a++)
        for (int c = a+1; c < j; c++)
          actual_max = std::max(actual_max, linf(bp[a].second, bp[c].second));
      i = j;
    }

    bool ok = actual_max <= bound;
    CHECK(ok, "coarsening: intra-bucket L∞ diameter ≤ 2^shiftBits - 1");
    printf("coarsening 2D shift=%d: max intra-bucket L∞ = %d (bound = %d): %s\n",
           shiftBits, actual_max, bound, ok?"OK":"FAIL");
  }
}

// ── 6. extreme dimensions ─────────────────────────────────────────────────────
// Verifies the worst-case TLPR topology (5x5, R=40) fits in uint64_t and that
// encode/decode round-trips at the all-zero and all-max boundary values in
// both the 5D origin and 10D destination subspaces.
void test_extreme_dimensions() {
  int R = 40, b = hilbert::bitsFor(R);
  CHECK(b == 6, "bitsFor(40)==6");
  CHECK(5*b  <= 64, "5D origin  (5x5 topology) fits in uint64_t");
  CHECK(10*b <= 64, "10D dest (5x S+/S- split) fits in uint64_t");

  for (int d : {5, 10}) {
    for (int val : {0, R}) {
      std::vector<int> coords(d, val);
      auto rec = hilbert::decode(hilbert::encode(coords, b), d, b);
      CHECK(rec == coords, "extreme: boundary round-trip");
    }
  }
  printf("extreme dimensions (5x5, R=%d, b=%d): OK\n", R, b);
}

int main() {
  test_bitsFor();
  test_bijectivity_2d();
  test_bijectivity_3d();
  test_unit_step_2d();
  test_unit_step_3d();
  test_ground_truth_2d();
  test_coarsening_consistency();
  test_extreme_dimensions();

  printf("\n%s (%d failure(s))\n",
         failures ? "SOME TESTS FAILED" : "All tests passed", failures);
  return failures ? 1 : 0;
}
