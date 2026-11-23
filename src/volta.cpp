// ---------------------------------------------------------------------------
// volta.cpp -- orderly generation of k-polycirculant graphs
//
// Copyright (C) 2026 Steven Van Overberghe
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License along
// with this program.  If not, see <https://www.gnu.org/licenses/>.
// ---------------------------------------------------------------------------
//
// One compile unit.  Everything that fixes the parameter set is a compile-time
// macro, so the program is compiled once per parameter set; the `volta` script
// next to the src/ directory does that (and caches the binaries), so normally
// there is no need to compile this file by hand.  By hand:
//
//   g++ -O3 -march=native -std=c++20 -pthread -DORBITS=5,5 volta.cpp -o volta-5-5
//   ./volta-5-5 -t 4
//
// Terminology: phi is the rotation (i,j) -> (i,j+1 mod m_i) of the orbits, and an
// "edge orbit" is always an orbit of <phi> on the vertex pairs (a bit of the word),
// never an orbit of the full automorphism group of a graph, which can be larger.
//
// Parameter set (compile time):
//   -DORBITS=m0,m1,...   orbit lengths (any order, sorted internally; at most 6
//                        orbits, every length divides the largest, which is <= 64)
//
// Forbidden subgraphs (doc/DESIGN.md 7.2, 7.7).  Each of the three is available
// for the graph and for its complement; all of them may be combined.
//   -DFORBID_CLIQUE=r    only K_r-free graphs (r >= 2)
//   -DFORBID_KME=r       only (K_r - e)-free graphs (3 <= r <= n): no complete graph
//                        on r vertices with one edge missing.  r = 3 is P_3 (the
//                        graphs are the matchings), r = 4 the diamond, i.e. no two
//                        triangles share an edge
//   -DFORBID_CYCLE=l     only C_l-free graphs (3 <= l <= n): no cycle of length
//                        exactly l (MINGIRTH below forbids every length < g instead)
//   -DFORBID_COCLIQUE=s  the same three in the complement: no independent set of
//   -DFORBID_COKME=r     size s (FORBID_CLIQUE=r with FORBID_COCLIQUE=s gives the
//   -DFORBID_COCYCLE=l   (r,s)-Ramsey graphs), no K_r - e and no C_l in the complement
//   -DFORWARD=1          forward checking on the undecided edge orbits (propagation
//                        of forced edges and forced non-edges, doc/DESIGN.md 7.3)
//
// Other properties:
//   -DMINDEG=d           only graphs of minimum degree >= d        (doc/DESIGN.md 7.4)
//   -DMAXDEG=d           only graphs of maximum degree <= d
//   -DREGULAR=d          only d-regular graphs (= MINDEG=d MAXDEG=d)
//   -DMAXCHI=c           only graphs of chromatic number <= c      (doc/DESIGN.md 7.5)
//   -DMINGIRTH=g         only graphs of girth >= g, i.e. without cycles shorter
//                        than g (g = 4: triangle-free, g = n+1: forests)  (7.6)
//
// Output and farming:
//   -DGRAPH6=1           write every graph to stdout, one graph6 line each (vertex
//                        (i,j) is vertex m_0+...+m_{i-1}+j); the summary then goes
//                        to stderr, as with nauty's geng.  Default: counts only.
//   -DPARTS=1            the binary accepts -p r/m and runs only the tasks with
//                        index r (mod m); the m parts together report every graph
//                        exactly once (doc/DESIGN.md 6)
//
// Run time options of the binary: -v verbose summary, -t threads, -s split depth,
// -p r/m (see usage()).
// ---------------------------------------------------------------------------

#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <array>
#include <vector>
#include <thread>
#include <atomic>
#include <chrono>
#include <algorithm>
#include <bit>
#include <cstring>

#include "bitset.h"

#ifndef VOLTA_VERSION
#define VOLTA_VERSION "1.0"
#endif
#ifndef ORBITS
#define ORBITS 5,5
#endif
#ifndef FORBID_CLIQUE
#define FORBID_CLIQUE 0     // r: no K_r in G (doc/DESIGN.md 7.2); 0 = off
#endif
#ifndef FORBID_COCLIQUE
#define FORBID_COCLIQUE 0   // s: no K_s in the complement (no independent set of size s)
#endif
#ifndef FORBID_KME
#define FORBID_KME 0        // r: no K_r minus an edge in G (doc/DESIGN.md 7.7); 0 = off
#endif
#ifndef FORBID_COKME
#define FORBID_COKME 0      // r: no K_r minus an edge in the complement
#endif
#ifndef FORBID_CYCLE
#define FORBID_CYCLE 0      // l: no cycle of length exactly l in G (doc/DESIGN.md 7.7); 0 = off
#endif
#ifndef FORBID_COCYCLE
#define FORBID_COCYCLE 0    // l: no cycle of length exactly l in the complement
#endif
#ifndef CLIQUE_SIMPLE_MAX
#define CLIQUE_SIMPLE_MAX 7 // largest forbidden clique handled by the plain clique search (7.2)
#endif
#ifndef FORWARD
#define FORWARD 0           // 1: forward checking (doc/DESIGN.md 7.3)
#endif
#ifndef MAXDEG
#define MAXDEG 0
#endif
#ifndef MINDEG
#define MINDEG 0
#endif
#ifndef REGULAR
#define REGULAR 0
#endif
#if REGULAR > 0
#undef MINDEG
#undef MAXDEG
#define MINDEG REGULAR
#define MAXDEG REGULAR
#endif
#ifndef GRAPH6
#define GRAPH6 0            // 1: write the graphs to stdout in graph6 format (doc/DESIGN.md 5)
#endif
#ifndef MAXCHI
#define MAXCHI 0            // c: only graphs of chromatic number <= c (doc/DESIGN.md 7.5); 0 = off
#endif
#ifndef MINGIRTH
#define MINGIRTH 0          // g: only graphs of girth >= g (doc/DESIGN.md 7.6); 0 = off
#endif
#ifndef PARTS
#define PARTS 0             // 1: run only the tasks r (mod m) given by -p r/m (doc/DESIGN.md 6)
#endif

namespace poly {

// ===========================================================================
// 1. Compile-time parameters and layout
// ===========================================================================

constexpr int RAWORB[] = { ORBITS };
constexpr int K = (int)(sizeof(RAWORB) / sizeof(RAWORB[0]));

constexpr int cgcd(int a, int b) { while (b) { int t = a % b; a = b; b = t; } return a; }

constexpr std::array<int, K> sortOrbits() {
    std::array<int, K> a{};
    for (int i = 0; i < K; ++i) a[i] = RAWORB[i];
    for (int i = 1; i < K; ++i) {            // insertion sort, non-increasing
        int v = a[i], j = i - 1;
        while (j >= 0 && a[j] < v) { a[j + 1] = a[j]; --j; }
        a[j + 1] = v;
    }
    return a;
}
constexpr std::array<int, K> MS = sortOrbits();   // m_0 >= m_1 >= ... >= m_{k-1}
constexpr int M0 = MS[0];                         // m

constexpr bool checkSizes() {
    for (int i = 0; i < K; ++i) { if (MS[i] < 1) return false; if (M0 % MS[i] != 0) return false; }
    return true;
}
static_assert(K >= 1 && K <= 6,  "k must satisfy 1 <= k <= 6");
static_assert(M0 <= 64,          "the largest orbit length must be at most 64");
static_assert(checkSizes(),      "every orbit length must be >= 1 and divide the largest one");

constexpr int NB = K * (K + 1) / 2;               // number of blocks

constexpr int blockLen(int a, int b) { return (a == b) ? MS[a] / 2 : cgcd(MS[a], MS[b]); }
constexpr int blockMod(int a, int b) { return (a == b) ? MS[a]     : cgcd(MS[a], MS[b]); }

constexpr int computeNPos() {
    int s = 0;
    for (int b = 0; b < K; ++b) for (int a = 0; a <= b; ++a) s += blockLen(a, b);
    return s;
}
constexpr int NPOS  = computeNPos();              // M, the length of the word
constexpr int NPOSA = NPOS > 0 ? NPOS : 1;

// Blocks are enumerated column by column: (0,0) | (0,1)(1,1) | (0,2)(1,2)(2,2) | ...
// Inside a block the word positions are s = 0..len-1 (off diagonal, s in Z_g) resp.
// s = 1..len (diagonal, the edge orbits {s,-s}).  A block is stored in one uint64_t
// with word position s in bit (len-1-s): the first word position is the most
// significant bit, so lexicographic order = unsigned integer order.
struct LayoutT {
    int      bidx[K][K];                          // block index of the pair (a,b)
    int      bA[NB], bB[NB], bLen[NB], bMod[NB], bStart[NB];
    int      posBlk[NPOSA], posA[NPOSA], posB[NPOSA];
    uint64_t posBit[NPOSA];
    int      posEdges[NPOSA], posDegA[NPOSA], posDegB[NPOSA];
    int      totalEdges;
};

constexpr LayoutT makeLayout() {
    LayoutT L{};
    int idx = 0, pos = 0;
    for (int b = 0; b < K; ++b) for (int a = 0; a <= b; ++a) {
        L.bidx[a][b] = idx; L.bidx[b][a] = idx;
        L.bA[idx] = a; L.bB[idx] = b;
        const int len = blockLen(a, b);
        L.bLen[idx] = len; L.bMod[idx] = blockMod(a, b); L.bStart[idx] = pos;
        for (int s = 0; s < len; ++s) {
            const int p = pos + s;
            L.posBlk[p] = idx; L.posA[p] = a; L.posB[p] = b;
            L.posBit[p] = uint64_t(1) << (len - 1 - s);
            if (a == b) {
                const int ss = s + 1;                       // word position s = 1..len
                const bool matching = (2 * ss == MS[a]);    // s = m/2: perfect matching
                L.posEdges[p] = matching ? MS[a] / 2 : MS[a];
                L.posDegA[p]  = matching ? 1 : 2;
                L.posDegB[p]  = 0;                          // posB == posA on the diagonal
            } else {
                const int g = cgcd(MS[a], MS[b]);
                L.posEdges[p] = MS[a] / g * MS[b];          // lcm(m_a,m_b) edges
                L.posDegA[p]  = MS[b] / g;
                L.posDegB[p]  = MS[a] / g;
            }
        }
        pos += len; ++idx;
    }
    L.totalEdges = 0;
    for (int p = 0; p < NPOS; ++p) L.totalEdges += L.posEdges[p];
    return L;
}
constexpr LayoutT L = makeLayout();

// Orbits of equal length form a contiguous interval: pi may only permute inside it.
constexpr std::array<int, K> makeLoi() {
    std::array<int, K> r{};
    for (int a = 0; a < K; ++a) { int lo = a; while (lo > 0 && MS[lo - 1] == MS[a]) --lo; r[a] = lo; }
    return r;
}
constexpr std::array<int, K> makeHii() {
    std::array<int, K> r{};
    for (int a = 0; a < K; ++a) { int hi = a; while (hi < K - 1 && MS[hi + 1] == MS[a]) ++hi; r[a] = hi; }
    return r;
}
constexpr std::array<int, K> LOI = makeLoi();
constexpr std::array<int, K> HII = makeHii();
// CANDM[b] = the orbits rho(b) may take: exactly the ones of the same length as m_b.
constexpr std::array<uint8_t, K> makeCandM() {
    std::array<uint8_t, K> r{};
    for (int b = 0; b < K; ++b) for (int i = LOI[b]; i <= HII[b]; ++i) r[b] |= uint8_t(1u << i);
    return r;
}
constexpr std::array<uint8_t, K> CANDM = makeCandM();

// --- vertices (for the forbidden-subgraph filters) ---------------------------
// Vertex (i,j) has number VOFF[i] + j; the orbits are consecutive intervals.
constexpr int computeNV() { int s = 0; for (int i = 0; i < K; ++i) s += MS[i]; return s; }
constexpr int NV = computeNV();
constexpr std::array<int, K + 1> makeVOff() {
    std::array<int, K + 1> r{};
    for (int i = 0; i < K; ++i) r[i + 1] = r[i] + MS[i];
    return r;
}
constexpr std::array<int, K + 1> VOFF = makeVOff();
constexpr int VBITS = NV <= 64 ? 64 : ((NV + 63) / 64) * 64;
using VSet = Bitset<VBITS>;                       // a set of vertices
using Adj  = std::array<VSet, NV>;                // adjacency rows

// --- the forbidden subgraphs, per side (0 = off) ---------------------------
// Three shapes, each of them for the graph of the decided edges (E) and for the
// graph of the decided non-edges (N): a clique K_r, a clique minus an edge
// K_r - e, and a cycle C_l.  The girth bound joins them on the E side (every
// cycle shorter than g).  All of them live in the same battery (2g).
constexpr int FC    = FORBID_CLIQUE;              // K_FC in G
constexpr int FKME  = FORBID_KME;                 // K_FKME - e in G
constexpr int FCYC  = FORBID_CYCLE;               // C_FCYC in G
constexpr int FCC   = FORBID_COCLIQUE;            // K_FCC in the complement
constexpr int FCKME = FORBID_COKME;               // K_FCKME - e in the complement
constexpr int FCCYC = FORBID_COCYCLE;             // C_FCCYC in the complement
static_assert(FC    == 0 || FC    >= 2, "FORBID_CLIQUE must be 0 (off) or at least 2");
static_assert(FCC   == 0 || FCC   >= 2, "FORBID_COCLIQUE must be 0 (off) or at least 2");
static_assert(FKME  == 0 || (FKME  >= 3 && FKME  <= NV), "FORBID_KME must be 0 (off) or between 3 and n");
static_assert(FCKME == 0 || (FCKME >= 3 && FCKME <= NV), "FORBID_COKME must be 0 (off) or between 3 and n");
static_assert(FCYC  == 0 || (FCYC  >= 3 && FCYC  <= NV), "FORBID_CYCLE must be 0 (off) or between 3 and n");
static_assert(FCCYC == 0 || (FCCYC >= 3 && FCCYC <= NV), "FORBID_COCYCLE must be 0 (off) or between 3 and n");
constexpr bool PART_MODE = PARTS != 0;            // -DPARTS=1: the binary runs only the tasks idx = r (mod m) of -p r/m (doc/DESIGN.md 6)
constexpr int  GIRTH  = MINGIRTH;                 // lower bound on the girth (0 = none)
static_assert(GIRTH >= 0 && GIRTH <= NV + 1, "MINGIRTH must be between 0 (off) and n+1 (forests)");
constexpr bool HAS_GIRTH = GIRTH >= 4;            // prune on short cycles (girth >= 3 holds trivially)
constexpr bool HAS_E  = FC  > 0 || FKME  > 0 || FCYC  > 0;   // a subgraph filter on the decided edges
constexpr bool HAS_N  = FCC > 0 || FCKME > 0 || FCCYC > 0;   // ... on the decided non-edges
constexpr bool E_SUB  = HAS_E || HAS_GIRTH;       // the E-side battery has something to test
constexpr bool FWD    = FORWARD != 0 && (E_SUB || HAS_N);   // forward checking (needs an orbit filter)
constexpr int  CHI     = MAXCHI;                  // bound on the chromatic number (0 = none)
static_assert(CHI == 0 || (CHI >= 2 && CHI <= 31), "MAXCHI must be 0 (off) or between 2 and 31");
constexpr bool HAS_CHI = CHI > 0;                 // prune on the chromatic number (2e)
constexpr bool BIP     = HAS_CHI && CHI == 2;     // bipartite: the incremental subgroup test
constexpr bool CHI_GEN = HAS_CHI && CHI >= 3;     // exact colouring of the decided edges
constexpr bool ADJ_E   = E_SUB || CHI_GEN;        // the adjacency of the decided edges is maintained

// --- units of Z_m ----------------------------------------------------------
constexpr int countUnits() { int c = 0; for (int t = 0; t < M0; ++t) if (cgcd(t, M0) == 1) ++c; return c; }
constexpr int NU = countUnits();
constexpr std::array<int, 64> makeUnits() {
    std::array<int, 64> u{}; int c = 0;
    for (int t = 0; t < M0; ++t) if (cgcd(t, M0) == 1) u[c++] = t;
    return u;
}
constexpr std::array<int, 64> UNITS = makeUnits();
constexpr std::array<int, 64> makeNeg() {           // index of -t
    std::array<int, 64> r{};
    for (int i = 0; i < NU; ++i) { const int v = (M0 - UNITS[i]) % M0;
        for (int j = 0; j < NU; ++j) if (UNITS[j] == v) r[i] = j; }
    return r;
}
constexpr std::array<int, 64> NEGU = makeNeg();

// --- distinct moduli -------------------------------------------------------
struct ModListT { int n; int val[NB > 0 ? NB : 1]; int ofOrbit[K]; int ofBlock[NB > 0 ? NB : 1]; };

constexpr ModListT makeDiagMods() {
    ModListT m{}; m.n = 0;
    for (int a = 0; a < K; ++a) {
        const int v = MS[a]; int f = -1;
        for (int j = 0; j < m.n; ++j) if (m.val[j] == v) f = j;
        if (f < 0) { f = m.n; m.val[m.n++] = v; }
        m.ofOrbit[a] = f;
    }
    return m;
}
constexpr ModListT makeOffMods() {
    ModListT m{}; m.n = 0;
    for (int bi = 0; bi < NB; ++bi) {
        if (L.bA[bi] == L.bB[bi]) { m.ofBlock[bi] = -1; continue; }
        const int v = L.bMod[bi]; int f = -1;
        for (int j = 0; j < m.n; ++j) if (m.val[j] == v) f = j;
        if (f < 0) { f = m.n; m.val[m.n++] = v; }
        m.ofBlock[bi] = f;
    }
    return m;
}
constexpr ModListT DIAGM = makeDiagMods();
constexpr ModListT OFFM  = makeOffMods();

constexpr int maxOf(int a, int b) { return a > b ? a : b; }
constexpr int maxDiagHalf() { int h = 0; for (int i = 0; i < DIAGM.n; ++i) h = maxOf(h, DIAGM.val[i] / 2); return h; }
constexpr int maxOffMod()   { int h = 0; for (int i = 0; i < OFFM.n;  ++i) h = maxOf(h, OFFM.val[i]);      return h; }

constexpr int NDM    = DIAGM.n;
constexpr int NOM    = maxOf(OFFM.n, 1);
constexpr int DBYTES = maxOf((maxDiagHalf() + 7) / 8, 1);
constexpr int OBYTES = maxOf((maxOffMod()   + 7) / 8, 1);

// ===========================================================================
// 2. Tables:  t*S by byte lookups, d mod g, block masks
// ===========================================================================

static uint64_t diagTab[NDM][NU][DBYTES][256];    // folded half word, modulus m_a
static uint64_t offTab [NOM][NU][OBYTES][256];    // full word, modulus gcd(m_a,m_b)
static uint8_t  MODT[65][64];                     // MODT[g][d] = d % g
static uint64_t GMASK[65];                        // low g bits

static void buildTables() {
    for (int g = 0; g <= 64; ++g) {
        GMASK[g] = (g >= 64) ? ~uint64_t(0) : (g == 0 ? uint64_t(0) : ((uint64_t(1) << g) - 1));
        for (int d = 0; d < 64; ++d) MODT[g][d] = uint8_t(g ? d % g : 0);
    }
    for (int mi = 0; mi < NDM; ++mi) {
        const int M = DIAGM.val[mi], h = M / 2;
        for (int ti = 0; ti < NU; ++ti) {
            const int t = UNITS[ti] % M;
            for (int by = 0; by < DBYTES; ++by) for (int bv = 0; bv < 256; ++bv) {
                uint64_t mask = 0;
                for (int q = 0; q < 8; ++q) if ((bv >> q) & 1) {
                    const int p = by * 8 + q; if (p >= h) continue;
                    const int s = h - p;                        // word position 1..h
                    int s2 = (s * t) % M; if (s2 > M - s2) s2 = M - s2;   // fold {x,-x}
                    mask |= uint64_t(1) << (h - s2);
                }
                diagTab[mi][ti][by][bv] = mask;
            }
        }
    }
    for (int mi = 0; mi < OFFM.n; ++mi) {
        const int M = OFFM.val[mi];
        for (int ti = 0; ti < NU; ++ti) {
            const int t = UNITS[ti] % M;
            for (int by = 0; by < OBYTES; ++by) for (int bv = 0; bv < 256; ++bv) {
                uint64_t mask = 0;
                for (int q = 0; q < 8; ++q) if ((bv >> q) & 1) {
                    const int p = by * 8 + q; if (p >= M) continue;
                    const int s  = M - 1 - p;                   // word position 0..M-1
                    const int s2 = (s * t) % M;
                    mask |= uint64_t(1) << (M - 1 - s2);
                }
                offTab[mi][ti][by][bv] = mask;
            }
        }
    }
}

static inline uint64_t mulDiag(int mi, int ti, uint64_t v) {
    uint64_t r = 0; int b = 0;
    while (v) { r |= diagTab[mi][ti][b][v & 0xFF]; v >>= 8; ++b; }
    return r;
}
static inline uint64_t mulOff(int mi, int ti, uint64_t v) {
    uint64_t r = 0; int b = 0;
    while (v) { r |= offTab[mi][ti][b][v & 0xFF]; v >>= 8; ++b; }
    return r;
}
// S -> S + r  is a cyclic right rotation in the reversed encoding.
static inline uint64_t rotr(uint64_t v, int r, int g) {
    return ((v >> r) | (v << ((g - r) & 63))) & GMASK[g];
}

// ===========================================================================
// 2b. Forbidden cliques: edge lists per edge orbit and the clique searches
// ===========================================================================

// The edges of every edge orbit p, as vertex pairs, consecutively: edges
// edgeStart[p] .. edgeStart[p+1]-1.  The first pair of an orbit is its
// representative ((a,0),(b,s)) resp. ((a,0),(a,s)).  Every pair of vertices lies
// in exactly one orbit, so the total is NV(NV-1)/2 = L.totalEdges.
constexpr int NPAIRS  = NV * (NV - 1) / 2;
constexpr int NPAIRSA = NPAIRS > 0 ? NPAIRS : 1;
static int      edgeStart[NPOSA + 1];
static uint16_t edgeU[NPAIRSA], edgeV[NPAIRSA];
// The neighbours of the two representative vertices u0, v0 of orbit p inside the
// orbit p itself (rows of the orbit's own adjacency), so that the common
// neighbourhood of (u0,v0) in adj + p is available without adding the orbit:
//   (adj[u0] | repRowU[p]) & (adj[v0] | repRowV[p]).
static VSet     repRowU[NPOSA], repRowV[NPOSA];
// The K_4-e test (2g) needs one more row of adj + p than those two: the row of
// the single vertex that closes a triangle on the representative edge.  Read
// from the orbit's pair list that is a scan of the orbit; as a table -- the row
// of every vertex in every orbit -- it is one lookup.  The table is built when
// the filter is on and it stays below ORBROW_MAX bytes (it is n times the pair
// lists, which is nothing on the sets that the filter is used on, and megabytes
// on the largest ones the program accepts).
#ifndef ORBROW_MAX
#define ORBROW_MAX (1 << 20)
#endif
constexpr bool ORBROW = (FKME == 4 || FCKME == 4)
                        && size_t(NPOSA) * NV * sizeof(VSet) <= size_t(ORBROW_MAX);
[[maybe_unused]] static VSet orbRowTab[ORBROW ? NPOSA : 1][ORBROW ? NV : 1];

[[maybe_unused]] static void buildEdgeLists() {
    int e = 0;
    for (int p = 0; p < NPOS; ++p) {
        edgeStart[p] = e;
        const int a = L.posA[p], b = L.posB[p];
        const int s = p - L.bStart[L.posBlk[p]] + (a == b ? 1 : 0);   // word position
        if (a == b) {
            const int m = MS[a];
            for (int j = 0; j < m; ++j) {
                if (2 * s == m && j >= m / 2) break;                   // matching: each edge once
                edgeU[e] = uint16_t(VOFF[a] + j); edgeV[e] = uint16_t(VOFF[a] + (j + s) % m); ++e;
            }
        } else {
            const int g = cgcd(MS[a], MS[b]);
            for (int j = 0; j < MS[a]; ++j)
                for (int t = 0; t < MS[b] / g; ++t) {
                    edgeU[e] = uint16_t(VOFF[a] + j);
                    edgeV[e] = uint16_t(VOFF[b] + (j + s) % g + t * g); ++e;
                }
        }
        if (e - edgeStart[p] != L.posEdges[p]) { fprintf(stderr, "edge list size mismatch\n"); abort(); }
        const int u0 = edgeU[edgeStart[p]], v0 = edgeV[edgeStart[p]];
        repRowU[p] = VSet{}; repRowV[p] = VSet{};
        for (int f = edgeStart[p]; f < e; ++f) {
            if (edgeU[f] == u0) repRowU[p].addBit(edgeV[f]);
            if (edgeV[f] == u0) repRowU[p].addBit(edgeU[f]);
            if (edgeU[f] == v0) repRowV[p].addBit(edgeV[f]);
            if (edgeV[f] == v0) repRowV[p].addBit(edgeU[f]);
        }
        if constexpr (ORBROW)
            for (int f = edgeStart[p]; f < e; ++f) {
                orbRowTab[p][edgeU[f]].addBit(edgeV[f]); orbRowTab[p][edgeV[f]].addBit(edgeU[f]);
            }
    }
    edgeStart[NPOS] = e;
    if (e != NPAIRS) { fprintf(stderr, "edge count mismatch\n"); abort(); }
}

[[maybe_unused]] static inline void addOrbit(Adj& adj, int p) {
    for (int e = edgeStart[p]; e < edgeStart[p + 1]; ++e) {
        adj[edgeU[e]].addBit(edgeV[e]); adj[edgeV[e]].addBit(edgeU[e]);
    }
}
[[maybe_unused]] static inline void removeOrbit(Adj& adj, int p) {
    for (int e = edgeStart[p]; e < edgeStart[p + 1]; ++e) {
        adj[edgeU[e]].removeBit(edgeV[e]); adj[edgeV[e]].removeBit(edgeU[e]);
    }
}
[[maybe_unused]] static inline void clearAdj(Adj& adj) { for (auto& r : adj) r = VSet{}; }
// The neighbours of the vertex w inside the edge orbit p (the row of w in the
// orbit alone), from the table or, when there is none, from the pair list.
[[maybe_unused]] static inline VSet orbitRow(int p, int w) {
    if constexpr (ORBROW) return orbRowTab[p][w];
    else {
        VSet r{};
        for (int e = edgeStart[p]; e < edgeStart[p + 1]; ++e) {
            if (edgeU[e] == w) r.addBit(edgeV[e]);
            if (edgeV[e] == w) r.addBit(edgeU[e]);
        }
        return r;
    }
}
// An absent graph (filter switched off) is an Empty; these overloads let the
// `if constexpr` branches that touch it compile without being instantiated.
struct Empty {};
[[maybe_unused]] static inline void addOrbit(Empty&, int) {}
[[maybe_unused]] static inline void removeOrbit(Empty&, int) {}
[[maybe_unused]] static inline void clearAdj(Empty&) {}

// --- does the graph adj contain a clique of size R inside the vertex set cand? ---

// (a) Plain recursion: pick a vertex, look for K_{R-1} among its neighbours in
//     cand, otherwise drop it.  Best for small R (nothing but ands and popcounts).
template<int R>
static inline bool hasKSimple(const Adj& adj, VSet cand) {
    if constexpr (R <= 0) return true;
    else if constexpr (R == 1) return cand.isNotEmpty();
    else {
        while (cand.popcount() >= R) {
            const int v = cand.findFirstSetBit();
            cand.removeBit(v);
            if (hasKSimple<R - 1>(adj, cand & adj[v])) return true;
        }
        return false;
    }
}

// Runtime-size entry to the plain recursion (need <= MAXR).
template<int MAXR>
static inline bool hasKSimpleRT(const Adj& adj, VSet cand, int need) {
    if constexpr (MAXR <= 0) { (void)adj; (void)cand; (void)need; return true; }
    else {
        if (need >= MAXR) return hasKSimple<MAXR>(adj, cand);
        return hasKSimpleRT<MAXR - 1>(adj, cand, need);
    }
}

// (b) Branch and bound with a greedy-colouring bound (Tomita's MCQ, decision
//     version): a greedy colouring of cand with fewer than `need` colours proves
//     that no K_need exists; otherwise the vertices are tried in decreasing
//     colour order, and as soon as the colour of the next vertex is below `need`
//     the remaining vertices (all of lower colour) cannot host a K_need.  The
//     last levels (need <= COLOUR_LEAF) and small candidate sets (fewer than
//     COLOUR_MINCAND vertices) go to the plain recursion, which is cheaper than a
//     colouring there.  MAXR bounds `need` (it is FS-2 at the top level).
#ifndef COLOUR_LEAF
#define COLOUR_LEAF 4
#endif
#ifndef COLOUR_MINCAND
#define COLOUR_MINCAND 16
#endif
template<int MAXR>
static bool hasKColour(const Adj& adj, VSet cand, int need) {
    if (need <= 0) return true;
    const int cnt = cand.popcount();
    if (cnt < need) return false;
    if (need <= COLOUR_LEAF || cnt < COLOUR_MINCAND) return hasKSimpleRT<MAXR>(adj, cand, need);
    int order[NV], col[NV];                       // up to NV colours
    int n = 0, k = 0;
    VSet U = cand;
    while (U.isNotEmpty()) {
        ++k;
        VSet Q = U;                                   // the vertices that may still get colour k
        while (Q.isNotEmpty()) {
            const int v = Q.findFirstSetBit();
            Q = Q & ~adj[v]; Q.removeBit(v);
            U.removeBit(v);
            order[n] = v; col[n] = k; ++n;
        }
    }
    if (k < need) return false;
    for (int i = n - 1; i >= 0; --i) {
        if (col[i] < need) return false;
        const int v = order[i];
        if (hasKColour<MAXR>(adj, cand & adj[v], need - 1)) return true;
        cand.removeBit(v);
    }
    return false;
}

// The forbidden clique has size FS; the search looks for K_{FS-2} in the common
// neighbourhood of the representative edge, so FS itself decides the algorithm.
template<int FS>
static inline bool hasK(const Adj& adj, VSet cand) {
    if constexpr (FS <= CLIQUE_SIMPLE_MAX) return hasKSimple<FS - 2>(adj, cand);
    else                                   return hasKColour<FS - 2>(adj, cand, FS - 2);
}

// --- (c) a clique minus an edge: K_R - e ------------------------------------
// K_R - e is a clique of R-2 vertices (its core) plus two vertices adjacent to
// all of it (its apexes; the missing edge is theirs).  So a graph contains it
// iff some clique of R-2 vertices has two common neighbours, and the search is
// the clique recursion above with the common neighbourhood carried along:
//   P  the vertices that may still extend the core -- a vertex once tried is
//      dropped, as in hasKSimple;
//   Q  the common neighbourhood of the core so far, where the apexes come from
//      -- nothing may be dropped from it.
// P is a subset of Q throughout, and NEED more core vertices plus two apexes,
// all of them distinct and in Q, need |Q| >= NEED + 2.
template<int NEED>
static bool hasKmeFrom(const Adj& adj, VSet P, VSet Q) {
    if constexpr (NEED <= 0) { (void)adj; (void)P; return Q.hasAtLeastWeight<2>(); }
    else {
        while (P.hasAtLeastWeight<NEED>() && Q.hasAtLeastWeight<NEED + 2>()) {
            const int v = P.findFirstSetBit();
            P.removeBit(v);
            if (hasKmeFrom<NEED - 1>(adj, P & adj[v], Q & adj[v])) return true;
        }
        return false;
    }
}
// Is there a K_R - e through the edge (u0,v0)?  W is the common neighbourhood of
// u0 and v0.  The edge is either inside the core or joins the core to an apex:
//   * both core: the other R-4 core vertices are adjacent to u0 and to v0, so
//     they lie in W, and so do the two apexes;
//   * u0 core, v0 apex: the other R-3 core vertices lie in W again (adjacent to
//     u0 as core members, to v0 as its neighbours), and the apexes lie in the
//     common neighbourhood of the core, which starts at N(u0) and holds v0 --
//     hence the "two apexes" test asks for one more vertex than v0 there.
template<int R>
static inline bool hasKmeThrough(const Adj& adj, int u0, int v0, VSet W) {
    if constexpr (R >= 4) { if (hasKmeFrom<R - 4>(adj, W, W)) return true; }
    return hasKmeFrom<R - 3>(adj, W, adj[u0]) || hasKmeFrom<R - 3>(adj, W, adj[v0]);
}

// --- (d) a cycle of a given length: C_L --------------------------------------
// A C_L through the edge (u0,v0) is that edge plus a path of L-1 edges from u0
// to v0 avoiding it: a depth-first search of fixed depth, one template level per
// edge.  The intermediate vertices avoid the visited set and the target, so the
// edge u0v0 is never used (u0 stays visited, v0 is entered only on the last
// step, from a vertex that is not u0 as soon as L >= 3).
template<int R>
static bool pathExact(const Adj& adj, int cur, int tgt, VSet vis) {
    if constexpr (R <= 1) { (void)vis; return adj[cur].hasBit(tgt); }
    else if constexpr (R == 2) return (adj[cur] & adj[tgt] & ~vis).isNotEmpty();
    else {
        VSet nxt = adj[cur] & ~vis;
        nxt.removeBit(tgt);
        for (int v : nxt.tracked()) {
            VSet vis2 = vis; vis2.addBit(v);
            if (pathExact<R - 1>(adj, v, tgt, vis2)) return true;
        }
        return false;
    }
}
template<int L>
static inline bool hasCycleThrough(const Adj& adj, int u0, int v0) {
    VSet vis{}; vis.addBit(u0);
    return pathExact<L - 1>(adj, u0, v0, vis);
}

// The adjacency of the decided edges (adjE) and of the decided non-edges (adjN),
// per thread, maintained incrementally by extend().  Both graphs are unions of
// edge orbits, hence phi-invariant; the filters that run on them are section 2g.
struct Forbid {
    [[no_unique_address]] std::conditional_t<ADJ_E, Adj, Empty> adjE;   // decided (and forced) edges
    [[no_unique_address]] std::conditional_t<HAS_N, Adj, Empty> adjN;

    // Forward checking (FWD): the two graphs then hold the decided AND the forced
    // orbits; state[q] = 0 undecided, 1 in adjE, 2 in adjN; every addition is
    // logged so that a node can undo its own additions on return.
    uint8_t state[NPOSA] = {};
    int16_t log[NPOSA]   = {};
    int     logN = 0;

    inline void pushE(int q) { addOrbit(adjE, q); state[q] = 1; log[logN++] = int16_t(q); }
    inline void pushN(int q) { addOrbit(adjN, q); state[q] = 2; log[logN++] = int16_t(q); }
    inline void noteE(int q) { state[q] = 1; log[logN++] = int16_t(q); }   // orbit already added
    inline void noteN(int q) { state[q] = 2; log[logN++] = int16_t(q); }
    inline void undo(int mark) {
        while (logN > mark) {
            const int q = log[--logN];
            if (state[q] == 1) removeOrbit(adjE, q); else removeOrbit(adjN, q);
            state[q] = 0;
        }
    }

    // Rebuild both graphs for a node taken from the task vector: positions
    // <= maxPos in X are edges, the others are non-edges.
    void rebuild(const uint64_t* blk, int maxPos) {
        clearAdj(adjE); clearAdj(adjN);
        for (int q = 0; q < NPOS; ++q) state[q] = 0;
        logN = 0;
        for (int q = 0; q <= maxPos; ++q) {
            if (blk[L.posBlk[q]] & L.posBit[q]) { addOrbit(adjE, q); state[q] = 1; }
            else                                { addOrbit(adjN, q); state[q] = 2; }
        }
    }
};

// ===========================================================================
// 2c. Degree bounds (-DMINDEG, -DMAXDEG, -DREGULAR)
// ===========================================================================
// Degrees are constant on orbits: deg[i] is the degree of every vertex of O_i,
// and every edge orbit adds a fixed amount to one or two of them (L.posDegA,
// L.posDegB).  Degrees only grow along a path of the search, so the upper bound
// is a test on the candidate, before its canonicity test (degMayAdd).  For the
// lower bound, lastP[i][r] is the largest position q such that the positions
// >= q can still add r to deg[i] (-1 if there is none): the candidate loop of a
// node stops at min_i lastP[i][MINDEG - deg[i]], since a later candidate could
// not lift orbit i to MINDEG any more, and neither could the node itself (which
// has all later positions as non-edges), so no graph is lost (doc/DESIGN.md 7.4).
constexpr int  DMIN     = MINDEG;
constexpr int  DMAX     = MAXDEG;                            // 0 = no upper bound
constexpr bool HAS_DEG  = DMIN > 0 || DMAX > 0;             // prune on degrees
static_assert(DMIN >= 0 && DMAX >= 0, "MINDEG and MAXDEG must be >= 0");
static_assert(DMAX == 0 || DMIN <= DMAX, "MINDEG must not exceed MAXDEG");
static_assert(DMIN == 0 || DMIN <= NV - 1, "MINDEG exceeds n-1");

struct DegT {
    int lastP[K][NV + 1];                // largest p such that the positions >= p add >= r to deg[i], or -1
};
constexpr DegT makeDegT() {
    DegT D{};
    for (int i = 0; i < K; ++i) {
        int gain[NPOSA + 1] = {};                                // gain[q] = what positions >= q add to deg[i]
        gain[NPOS] = 0;
        for (int q = NPOS - 1; q >= 0; --q)
            gain[q] = gain[q + 1] + (L.posA[q] == i ? L.posDegA[q] : 0) + (L.posB[q] == i ? L.posDegB[q] : 0);
        int p = NPOS - 1;
        for (int r = 0; r <= NV; ++r) { while (p >= 0 && gain[p] < r) --p; D.lastP[i][r] = p; }
    }
    return D;
}
constexpr DegT DT = makeDegT();

// The last candidate position of a node's loop, kept per depth in the thread's
// scratch (Canon::dg[nOrb]) so that Node itself is untouched: a child writes its
// own level, the parent's stays valid.
struct DegState {
    int16_t  degQ;        // last candidate position (later ones cannot reach MINDEG)
};
using DegArr = std::array<DegState, NPOSA + 1>;

// ===========================================================================
// 2d. graph6 output (-DGRAPH6=1)
// ===========================================================================
// Every reported graph is written to stdout as one graph6 line; the summary then
// goes to stderr (as with geng).  Vertex (i,j) is vertex VOFF[i] + j, so the
// orbits are the consecutive intervals 0..m_0-1, m_0..m_0+m_1-1, ...
// The body of a graph6 line is the upper triangle of the adjacency matrix, column
// by column, six bits per byte (the pair (u,v), u < v, is bit v(v-1)/2 + u of the
// stream, most significant bit first).  It is an OR: every edge orbit p has a fixed
// byte pattern g6orb[p] (the bits of its pairs), the body of G_X is the OR of the
// patterns of the orbits in X, plus 63 on every byte.  The patterns are kept as
// 64-bit words (their bytes in stream order), so a graph costs one word-OR per
// orbit and word; the +63 is one word addition, carry-free because every byte is
// at most 63.  Each thread fills its own buffer and hands it to stdout in one
// fwrite, which POSIX makes atomic against the other threads' stdio calls: lines
// never interleave, but with several threads their order depends on the schedule.
// Nothing of this exists in a build without GRAPH6: the writer is an empty member
// of Counters (G6Off, a type of its own so that it overlaps with the other empty
// member) and the calls are empty overloads.
constexpr bool G6     = GRAPH6 != 0;
constexpr int  G6HDR  = NV < 63 ? 1 : 4;                    // header bytes (NV <= 384 < 2^18)
constexpr int  G6BODY = (NPAIRS + 5) / 6;                   // body bytes
constexpr int  G6W    = (G6BODY + 7) / 8;                   // body words
constexpr int  G6LINE = G6HDR + G6BODY + 1;                 // with the newline
constexpr int  G6BUF  = G6LINE > (1 << 16) ? G6LINE : (1 << 16);   // per-thread buffer
static uint64_t g6orb[G6 ? NPOSA : 1][G6 ? G6W : 1];        // body bits of the pairs of orbit p

[[maybe_unused]] static void buildG6Orbits() {               // after buildEdgeLists()
    for (int p = 0; p < NPOS; ++p) {
        uint8_t* bytes = reinterpret_cast<uint8_t*>(g6orb[p]);   // byte i of the words = stream byte i
        for (int e = edgeStart[p]; e < edgeStart[p + 1]; ++e) {
            int u = edgeU[e], v = edgeV[e]; if (u > v) std::swap(u, v);
            const int k = v * (v - 1) / 2 + u;                  // bit k of the body
            bytes[k / 6] |= uint8_t(0x20 >> (k % 6));
        }
    }
}

struct G6Writer {
    char buf[G6BUF];
    int  len = 0;
    G6Writer() = default;
    G6Writer(const G6Writer&) = delete;                      // the buffer is not to be duplicated
    G6Writer& operator=(const G6Writer&) = delete;
    void flush() { if (len) { fwrite(buf, 1, size_t(len), stdout); len = 0; } }
    // The graph of X, given as the (full) word blk: its set bits are the edge orbits.
    void write(const uint64_t* blk) {
        if (len + G6LINE > G6BUF) flush();
        char* p = buf + len;
        if constexpr (G6HDR == 1) { *p++ = char(63 + NV); }
        else { *p++ = char(126); *p++ = char(63 + ((NV >> 12) & 63)); *p++ = char(63 + ((NV >> 6) & 63)); *p++ = char(63 + (NV & 63)); }
        uint64_t body[G6W];
        for (int i = 0; i < G6W; ++i) body[i] = 0;
        for (int bi = 0; bi < NB; ++bi) {
            for (uint64_t v = blk[bi]; v; v &= v - 1) {
                const int q = L.bStart[bi] + L.bLen[bi] - 1 - std::countr_zero(v);   // bit len-1-s is position s
                for (int i = 0; i < G6W; ++i) body[i] |= g6orb[q][i];
            }
        }
        for (int i = 0; i < G6W; ++i) body[i] += 0x3F3F3F3F3F3F3F3Full;   // +63 on every byte
        memcpy(p, body, G6BODY);
        p += G6BODY;
        *p++ = '\n';
        len = int(p - buf);
    }
};
struct G6Off {};                                             // the writer of a build without GRAPH6
[[maybe_unused]] static inline void g6Write(G6Writer& w, const uint64_t* blk) { w.write(blk); }
[[maybe_unused]] static inline void g6Write(G6Off&, const uint64_t*) {}
[[maybe_unused]] static inline void g6Flush(G6Writer& w) { w.flush(); }
[[maybe_unused]] static inline void g6Flush(G6Off&) {}

// ===========================================================================
// 2e. Chromatic number (-DMAXCHI=c)
// ===========================================================================
// chi(G) <= c is closed under edge deletion, so it is tested on G_X at every node
// and a failure cuts the subtree (doc/DESIGN.md 7.1, 7.5).
//
// c = 2, bipartite (BIP): no adjacency and no colouring.  Vertex (i,j) of G_X lifts
// to the vertices (i,x), x = j (mod m_i), of the graph G^ on k x Z_m in which
// (i,x) ~ (i',x') iff (x'-x) mod g_{ii'} in S_{ii'}.  The projection G^ -> G_X is a
// homomorphism along which every walk lifts, and a closed walk at (i,j) lifts to a
// walk from (i,x) to (i,x + t m_i).  In G^ the translation x -> x+1 is an
// automorphism, so G^ is the derived graph of the orbit graph Q (vertices = the
// orbits, an edge {i,i'} when S_{ii'} is non-empty, a loop at i when S_{ii} is)
// with the voltage *set* S_{ii'} + g_{ii'} Z_m on the dart i -> i'.  For a component
// of Q let W be the subgroup of Z_m x Z_2 of the pairs (net voltage, length mod 2)
// of its closed walks (abelian, hence independent of the base vertex).  G_X has an
// odd closed walk through orbit i iff (t m_i, 1) in W for some t; and since W
// contains (g, 0) for every edge of the component (one dart there with voltage v,
// back with -(v+g)) and g divides m_i, that is the case iff (0,1) in W.  So
//
//     G_X is bipartite  iff  (0,1) lies in the W of no component of Q.
//
// Relative to potentials (a_i, b_i) in Z_m x Z_2 that are consistent along a
// spanning tree, W is generated by (a_i + s - a_i', b_i + b_i' + 1) for s in S_{ii'}
// and (g_{ii'}, 0), and for a loop by (s, 1) and (m_i, 0).  A subgroup of Z_m x Z_2
// without (0,1) is the graph of a homomorphism on its projection <d>, i.e. a pair
// (d, chi) with d | m and chi = chi(d) = 1 only if m/d is even; adding a generator is
// one gcd and a parity case distinction (bipGen).  The state of a node (components,
// potentials, one (d, chi) per component: 5k bytes) lives per depth in the thread's
// scratch (Canon::bp) and is derived from the parent's by the generator, or the
// component merge, of the added orbit -- *before* the canonicity test, so a
// candidate that would close an odd cycle costs a 30-byte copy and a gcd.
//
// c >= 3 (CHI_GEN): the adjacency of the decided edges (Forbid::adjE, maintained as
// for FORBID_CLIQUE) and one c-colouring per depth.  A child inherits its parent's
// colouring when it is still proper; otherwise an exact backtracking search runs
// after the canonicity test (findColouring below).  Correct, not tuned.

constexpr std::array<std::array<uint8_t, 65>, 65> makeGcd() {
    std::array<std::array<uint8_t, 65>, 65> t{};
    for (int a = 0; a <= 64; ++a) for (int b = 0; b <= 64; ++b) t[a][b] = uint8_t(cgcd(a, b));
    return t;
}
constexpr auto GCDT = makeGcd();                             // gcd, with gcd(d, 0) = d

struct BipOff {};                                            // the state of a build without MAXCHI=2
struct BipState {
    uint8_t comp[K];   // component of the orbit, named by one of its orbits
    uint8_t pot[K];    // a_i in Z_m
    uint8_t par[K];    // b_i in Z_2
    uint8_t d[K];      // per component name: W = {(t d, t chi)}, d | m; d = m is the trivial group
    uint8_t chi[K];
};
using BipArr = std::array<BipState, NPOSA + 1>;

// W = {(t d, t chi)} without (0,1); add the generator (v, l), 0 <= v < m.  With
// d2 = gcd(d, v), d = alpha d2 and v = beta d2 (alpha, beta coprime, so not both
// even), the new group avoids (0,1) iff some c = chi'(d2) satisfies alpha c = chi
// and beta c = l (mod 2), and c = 1 only if m/d2 is even.  False: (0,1) entered.
static inline bool bipGen(uint8_t& d, uint8_t& chi, int v, int l) {
    const int d2 = GCDT[d][v];
    const int alpha = d / d2, beta = v / d2;
    int c;
    if (alpha & 1) { if ((beta & 1) ? (chi != l) : (l != 0)) return false; c = chi; }
    else           { if (chi) return false; c = l; }
    if (c && ((M0 / d2) & 1)) return false;
    d = uint8_t(d2); chi = uint8_t(c);
    return true;
}

// The state of X + {p} from the state of X, in place.  False: X + {p} is not bipartite.
static inline bool bipApply(BipState& st, int p) {
    const int bi = L.posBlk[p], a = L.posA[p], b = L.posB[p];
    const int s = p - L.bStart[bi] + (a == b ? 1 : 0);        // word position
    if (a == b) {                                              // loop: (s, 1) and (m_a, 0)
        const int ca = st.comp[a];
        return bipGen(st.d[ca], st.chi[ca], s, 1) && bipGen(st.d[ca], st.chi[ca], MS[a] % M0, 0);
    }
    const int ca = st.comp[a], cb = st.comp[b], g = L.bMod[bi];
    if (ca != cb) {
        // Two components: shift B's potentials so that the dart (a -> b, s) is a tree
        // edge (its generator vanishes), rename B, and W = W_A + W_B.
        const int delta = (st.pot[a] + s + M0 - st.pot[b]) % M0, eps = st.par[a] ^ st.par[b] ^ 1;
        for (int j = 0; j < K; ++j) if (st.comp[j] == cb) {
            st.pot[j] = uint8_t((st.pot[j] + delta) % M0); st.par[j] ^= uint8_t(eps); st.comp[j] = uint8_t(ca);
        }
        if (!bipGen(st.d[ca], st.chi[ca], st.d[cb] % M0, st.chi[cb])) return false;
    } else {
        const int v = (st.pot[a] + s + M0 - st.pot[b]) % M0, l = st.par[a] ^ st.par[b] ^ 1;
        if (!bipGen(st.d[ca], st.chi[ca], v, l)) return false;
    }
    return bipGen(st.d[ca], st.chi[ca], g % M0, 0);          // two voltages of one dart: (g, 0)
}

// The state of a node from its word (a cold start: the root, or a task taken from the
// task vector): the empty word's state, then the positions of X in order.
[[maybe_unused]] static inline void bipRebuild(BipArr& a, const uint64_t* blk, int maxPos, int nOrb) {
    BipState st{};
    for (int i = 0; i < K; ++i) { st.comp[i] = uint8_t(i); st.d[i] = uint8_t(M0); }
    for (int q = 0; q <= maxPos; ++q) if (blk[L.posBlk[q]] & L.posBit[q])
        if (!bipApply(st, q)) { fprintf(stderr, "bipartite state: a task is not bipartite\n"); abort(); }
    a[nOrb] = st;
}
[[maybe_unused]] static inline void bipRebuild(BipOff&, const uint64_t*, int, int) {}

// c >= 3.  Is the subgraph of adj induced on U C-colourable?  Exact backtracking;
// the uncoloured vertex with the fewest admissible colours is coloured next (a
// colour is admissible when the vertex has no neighbour in its class), and since
// the empty classes are interchangeable only the first of them is tried.
// cls[0..C-1] are the classes, `used` of them non-empty.
template<int C>
static bool colourRec(const Adj& adj, VSet U, std::array<VSet, C>& cls, int used) {
    if (U.isEmpty()) return true;
    int bestV = -1, bestN = C + 1;
    uint32_t bestM = 0;
    for (int v : U.tracked()) {
        uint32_t m = 0; int n = 0;
        for (int c = 0; c < used; ++c) if ((adj[v] & cls[c]).isEmpty()) { m |= 1u << c; ++n; }
        if (used < C) { m |= 1u << used; ++n; }
        if (n < bestN) { bestN = n; bestV = v; bestM = m; if (n == 0) return false; }
    }
    U.removeBit(bestV);
    for (uint32_t m = bestM; m; m &= m - 1) {
        const int c = std::countr_zero(m);
        cls[c].addBit(bestV);
        if (colourRec<C>(adj, U, cls, c == used ? used + 1 : used)) return true;
        cls[c].removeBit(bestV);
    }
    return false;
}
// Is the subgraph induced on U bipartite?  BFS two-colouring by levels, on bitsets.
[[maybe_unused]] static bool bipartiteSet(const Adj& adj, VSet U) {
    while (U.isNotEmpty()) {
        VSet A{}, B{}, cur{};
        cur.addBit(U.findFirstSetBit()); A = cur;
        bool inA = true;                                        // the colour of cur
        while (cur.isNotEmpty()) {
            VSet nb{};
            for (int v : cur.tracked()) nb |= adj[v];
            nb &= U;
            if ((nb & (inA ? A : B)).isNotEmpty()) return false;   // an edge inside a colour class
            cur = nb & ~(A | B);
            if (inA) B |= cur; else A |= cur;
            inA = !inA;
        }
        U &= ~(A | B);
    }
    return true;
}
template<int C>
[[maybe_unused]] static bool colourableSet(const Adj& adj, VSet U) {
    if constexpr (C == 2) return bipartiteSet(adj, U);
    else { std::array<VSet, C> cls{}; return colourRec<C>(adj, U, cls, 0); }
}
// A C-colouring of adj into cls, or false.  Necessary first: the neighbourhood of
// every vertex must be (C-1)-colourable, and by phi-invariance one vertex per orbit
// suffices; that is cheap (for C = 3 a BFS) and refutes most non-colourable graphs
// at once, where the exact search would prove non-colourability slowly.  Isolated
// vertices are left out of the search and put into class 0 afterwards.
constexpr int CHIA = CHI > 0 ? CHI : 1;
using Colouring = std::array<VSet, CHIA>;
using ChiArr    = std::array<Colouring, NPOSA + 1>;          // one colouring per depth
struct ColOff {};                                            // the colourings of a build without MAXCHI>=3
template<int C>
static bool findColouring(const Adj& adj, Colouring& out) {
    for (int i = 0; i < K; ++i) if (!colourableSet<C - 1>(adj, adj[VOFF[i]])) return false;
    VSet U{}, all{};
    for (int v = 0; v < NV; ++v) { all.addBit(v); if (adj[v].isNotEmpty()) U.addBit(v); }
    std::array<VSet, C> cls{};
    if (!colourRec<C>(adj, U, cls, 0)) return false;
    for (int c = 0; c < C; ++c) { out[c] = cls[c]; all &= ~cls[c]; }
    out[0] |= all;
    return true;
}
// Is the colouring proper for adj?  (No edge inside a class.)
template<int C>
static inline bool properColouring(const Adj& adj, const Colouring& cls) {
    for (int c = 0; c < C; ++c) for (int v : cls[c].tracked()) if ((adj[v] & cls[c]).isNotEmpty()) return false;
    return true;
}
// The child at depth nOrb: its parent's colouring usually survives the added orbit
// (then it is copied: NV word-ands), otherwise the exact search runs, and its
// result is the child's colouring.  False: the child is not C-colourable.
template<int C>
[[maybe_unused]] static inline bool chiCheck(ChiArr& a, int nOrb, const Adj& adj) {
    if (properColouring<C>(adj, a[nOrb - 1])) { a[nOrb] = a[nOrb - 1]; return true; }
    return findColouring<C>(adj, a[nOrb]);
}
template<int C, class A> [[maybe_unused]] static inline bool chiCheck(ColOff&, int, const A&) { return true; }
// A cold start (the root, or a task from the task vector): colour the graph outright.
[[maybe_unused]] static inline void chiRebuild(ChiArr& a, int nOrb, const Adj& adj) {
    if (!findColouring<CHIA>(adj, a[nOrb])) { fprintf(stderr, "colouring: a task is not colourable\n"); abort(); }
}
template<class A> [[maybe_unused]] static inline void chiRebuild(ColOff&, int, const A&) {}

// ===========================================================================
// 2f. Girth (-DMINGIRTH=g)
// ===========================================================================
// girth(G) >= g, i.e. no cycle of length < g, is closed under edge deletion, so it
// is tested on G_X in the prefilter slot (doc/DESIGN.md 7.1, 7.6), on the E-adjacency
// of the subgraph filters.  E is phi-invariant and has no short cycle (its parent
// passed), so a cycle of length < g in E + p contains an edge of p, and some
// phi-image of it contains the representative edge (u0,v0) of p.  A cycle of
// length L through u0v0 is that edge plus a path of length L-1 from u0 to v0
// avoiding it, hence
//
//     E + p has a cycle of length < g   iff   dist(u0,v0) <= g-2 in (E + p) - u0v0.
//
// A bidirectional BFS on bitsets decides this: the balls A around u0 and B around
// v0 grow alternately by one level (one row-OR per frontier vertex), and
// dist <= radius(A) + radius(B) iff the balls meet.  Level 1 leaves the edge out
// (N(u0) - v0, N(v0) - u0); the balls meet as soon as v0 enters A or u0 enters B,
// so neither endpoint is ever expanded and the removed edge is never used.  The
// two radius-1 balls meet iff u0 and v0 have a common neighbour, i.e. iff the
// orbit closes a triangle: that is stage 1 of the battery (2g), which runs from
// the two precomputed rows before the orbit is added -- so g = 4 never touches
// the adjacency, and for g >= 5 a candidate that closes a triangle costs no
// add/remove and the BFS below starts at radius 1 with two disjoint balls.
// With forward checking the forced edges pass the same test before they enter E,
// so E never holds a short cycle, and an undecided orbit that would close one is
// a forced non-edge.

// Does adj have a cycle of length < G through the edge (u0,v0)?  The edge is in
// adj, and the radius-1 balls fa = N(u0) - v0, fb = N(v0) - u0 are given and
// known to be disjoint (no triangle), so G >= 5 and the BFS starts at radius 1.
template<int G>
static inline bool shortCycleIn(const Adj& adj, int u0, int v0, VSet fa, VSet fb) {
    static_assert(G >= 5, "the triangle case is stage 1 of the battery");
    VSet A = fa, B = fb;                              // the balls
    A.addBit(u0); B.addBit(v0);
    for (int ra = 1, rb = 1; ra + rb < G - 2; ) {
        if (ra <= rb) {                               // grow A by one level
            VSet nb{};
            for (int v : fa.tracked()) nb |= adj[v];
            fa = nb & ~A;
            if (fa.isEmpty()) return false;           // A is the whole component of u0, and v0 is not in it
            A |= fa; ++ra;
            if ((fa & B).isNotEmpty()) return true;
        } else {                                      // grow B
            VSet nb{};
            for (int v : fb.tracked()) nb |= adj[v];
            fb = nb & ~B;
            if (fb.isEmpty()) return false;
            B |= fb; ++rb;
            if ((fb & A).isNotEmpty()) return true;
        }
    }
    return false;                                     // the balls of radii ra + rb = G-2 are disjoint
}

// ===========================================================================
// 2g. The forbidden-subgraph battery
// ===========================================================================
// One predicate per side: would adding the edge orbit p to adj create a graph
// the filters of that side forbid?  On the E side those are a clique K_CL, a
// clique minus an edge K_KME - e, a cycle C_CYC and -- from the girth bound --
// every cycle shorter than GTH; on the N side the same three (a girth bound has
// no complement form).  A zero switches a filter off; with all of them off the
// battery is an orbit addition, and a build without any filter never calls it.
//
// The battery is one test *through the representative edge*: adj is
// phi-invariant and contains none of the forbidden graphs (the parent was
// checked), so a new copy in adj + p contains an edge of p, and some phi-image
// of that copy contains the representative edge (u0,v0) of p (2b).  Every filter
// is therefore a question about the neighbourhood of that one edge:
//
//   K_CL       a K_{CL-2} in W = N(u0) & N(v0)                          (2b a,b)
//   K_KME - e  a core through the edge with two common neighbours       (2b c)
//   C_CYC      a path of CYC-1 edges from u0 to v0 avoiding the edge    (2b d)
//   girth      such a path of at most GTH-2 edges, by bidirectional BFS (2f)
//
// Stage 1 asks only the two rows of the representative edge, which repRowU /
// repRowV give *before* the orbit is added: N(u0) and N(v0) in adj + p, and
// their intersection W -- the vertices that close a triangle on the edge.  That
// settles K_3, C_3 and the triangle of the girth test (W non-empty), K_3 - e (a
// second edge at u0 or at v0), and K_4 - e but for the single case |W| = 1,
// which one more row -- the row of that one vertex, from orbitRow (2b) --
// finishes.  So a K_4-e filter never touches the adjacency, which is what makes
// it cheap enough to run on every undecided orbit under FORWARD=1.
//
// Stage 2 holds the searches that need the other rows of adj + p.  It adds the
// orbit once for all of them and removes it again when the answer is "yes", or
// when it is "no" and keep is false; with keep = true a passed orbit stays in
// adj, which is where the caller wants it (it descends into the child).
template<int CL, int KME, int CYC, int GTH>
static inline bool wouldForbid(Adj& adj, int p, bool keep) {
    if constexpr (CL == 2) { (void)adj; (void)p; (void)keep; return true; }   // no edge at all is allowed
    else {
        const int u0 = edgeU[edgeStart[p]], v0 = edgeV[edgeStart[p]];
        [[maybe_unused]] const VSet nu = adj[u0] | repRowU[p];    // N(u0) in adj + p
        [[maybe_unused]] const VSet nv = adj[v0] | repRowV[p];    // N(v0) in adj + p
        [[maybe_unused]] const VSet W  = nu & nv;                 // the triangles on (u0,v0)
        // --- stage 1: the two rows of the representative edge --------------
        if constexpr (KME == 3) {                                 // K_3 - e = P_3: a second edge at an end
            if (nu.hasAtLeastWeight<2>() || nv.hasAtLeastWeight<2>()) return true;
        }
        if constexpr (CL == 3 || CYC == 3 || GTH >= 4) {          // a triangle
            if (W.isNotEmpty()) return true;
        }
        if constexpr (KME == 4) {                                 // the diamond: two triangles on one edge
            if (W.isNotEmpty()) {
                if (W.hasAtLeastWeight<2>()) return true;         // both apexes sit on the edge itself
                const int  w  = W.findFirstSetBit();              // the edge lies in the one triangle u0 v0 w:
                const VSet nw = adj[w] | orbitRow(p, w);          // is one of its other two edges in a second one?
                if ((nu & nw).hasAtLeastWeight<2>() || (nv & nw).hasAtLeastWeight<2>()) return true;
            }
        }
        // --- stage 2: the searches that need the other rows of adj + p -----
        constexpr bool S2_CL  = CL  >= 4;                         // a K_{CL-2} in W
        constexpr bool S2_KME = KME >= 5;                         // a core of KME-3 core vertices in W, at least
        constexpr bool S2_CYC = CYC >= 4;                         // a path search
        constexpr bool S2_GTH = GTH >= 5;                         // the BFS, from the disjoint radius-1 balls
        if constexpr (S2_CL || S2_KME || S2_CYC || S2_GTH) {
            bool run = S2_CYC || S2_GTH;                          // the path searches always have work
            if constexpr (S2_CL)  run = run || W.hasAtLeastWeight<CL - 2>();
            if constexpr (S2_KME) run = run || W.hasAtLeastWeight<KME - 3>();
            if (run) {
                addOrbit(adj, p);
                bool r = false;
                if constexpr (S2_CL)  { if (!r && W.hasAtLeastWeight<CL  - 2>()) r = hasK<CL>(adj, W); }
                if constexpr (S2_KME) { if (!r && W.hasAtLeastWeight<KME - 3>()) r = hasKmeThrough<KME>(adj, u0, v0, W); }
                if constexpr (S2_CYC) { if (!r) r = hasCycleThrough<CYC>(adj, u0, v0); }
                if constexpr (S2_GTH) { if (!r) {
                    VSet fa = nu, fb = nv;                        // the radius-1 balls, without the edge
                    fa.removeBit(v0); fb.removeBit(u0);
                    r = shortCycleIn<GTH>(adj, u0, v0, fa, fb);
                } }
                if (r || !keep) removeOrbit(adj, p);
                return r;
            }
        }
        if (keep) addOrbit(adj, p);
        return false;
    }
}
template<int CL, int KME, int CYC, int GTH>
[[maybe_unused]] static inline bool wouldForbid(Empty&, int, bool) { return false; }

// The two sides as the search loops call them.
template<class A> [[maybe_unused]] static inline bool wouldForbidE(A& adj, int p, bool keep) {
    return wouldForbid<FC, FKME, FCYC, HAS_GIRTH ? GIRTH : 0>(adj, p, keep);
}
template<class A> [[maybe_unused]] static inline bool wouldForbidN(A& adj, int p, bool keep) {
    return wouldForbid<FCC, FCKME, FCCYC, 0>(adj, p, keep);
}

// ===========================================================================
// 3. Canonicity test
// ===========================================================================

// Per-thread scratch.  The two caches live across nodes and validate themselves
// by storing the block value they were computed from: between two siblings of the
// search tree exactly one block changes, so almost every entry stays valid.
// The caches carry no "initialised" flag: they are *pre-loaded with the entry for
// the empty word* (which is exactly the state of the root), so a zero `src` is a
// legitimate, valid key from the very first probe on.  A probe is then one load
// and one comparison, and both entries have a power-of-two size, so the index
// arithmetic is a shift instead of a multiply.
struct Canon {
    struct DInfo {                // t * S_{ii}          -- 16 bytes
        uint64_t src, img;
    };
    struct RInfo {                // t * Srel(i,j) and its rotations -- 32 bytes
        uint64_t src;             // block value it was built from
        uint64_t base;            // t * Srel(i,j)
        uint64_t maxv;            // lexicographically largest rotation of base
        uint64_t argmask;         // { r : rotr(base,r) == maxv }
    };
    struct L0Info {           // outcome of level 0
        uint64_t src[K];      // the diagonal blocks it was computed from
        uint8_t reject; int n;
        uint8_t t[NU * K], i0[NU * K];
    };
    const uint64_t* blk;      // target word (the current set X), block by block
    int      LB;              // index of the block containing max X
    int      curT;            // index of the current multiplier t
    int      rho[K];          // pi^{-1}
    int      dsh[K];          // d_a = c_{rho(a)} - c_{rho(0)}
    L0Info   l0;
    DInfo    dcache[NU][K];
    RInfo    rcache[NU][K][K];
    [[no_unique_address]] Forbid fb;   // forbidden-subgraph state (empty when no filter is on)
    [[no_unique_address]] std::conditional_t<HAS_DEG, DegArr, Empty> dg;   // degree state per depth (7.4)
    [[no_unique_address]] std::conditional_t<BIP, BipArr, BipOff> bp;      // bipartite state per depth (7.5)
    [[no_unique_address]] std::conditional_t<CHI_GEN, ChiArr, ColOff> ch;  // colouring per depth (7.5)

    // Entry of every cache for the empty word:  t*{} = {}, every rotation of {}
    // is {}, so maxv = 0 and argmask = all g rotations; level 0 accepts every
    // (t, i0) pair.  Must be called after buildTables() (it reads GMASK).
    Canon() : blk(nullptr), LB(0), curT(0), rho{}, dsh{} {
        for (int j = 0; j < K; ++j) l0.src[j] = 0;
        l0.reject = 0; l0.n = 0;
        for (int ti = 0; ti < NU; ++ti) for (int i0 = LOI[0]; i0 <= HII[0]; ++i0) {
            l0.t[l0.n] = uint8_t(ti); l0.i0[l0.n] = uint8_t(i0); ++l0.n;
        }
        for (int ti = 0; ti < NU; ++ti) {
            for (int i = 0; i < K; ++i) { dcache[ti][i].src = 0; dcache[ti][i].img = 0; }
            for (int i = 0; i < K; ++i) for (int j = 0; j < K; ++j) {
                RInfo& e = rcache[ti][i][j];
                e.src = 0; e.base = 0; e.maxv = 0;
                e.argmask = GMASK[i == j ? 0 : L.bMod[L.bidx[i][j]]];
            }
        }
    }
};

static inline uint64_t getDiag(Canon& c, int i) {
    Canon::DInfo& e = c.dcache[c.curT][i];
    const uint64_t src = c.blk[L.bidx[i][i]];
    if (e.src == src) [[likely]] return e.img;
    e.src = src; e.img = mulDiag(DIAGM.ofOrbit[i], c.curT, src);
    return e.img;
}

// t * Srel(i,j) together with its largest rotation and the argmax set.
// Srel(i,j) = S_{ij} for i<j and -S_{ji} for i>j; t*(-S) = (-t)*S, so one table
// serves both orientations.
static inline const Canon::RInfo& getRot(Canon& c, int i, int j) {
    Canon::RInfo& e = c.rcache[c.curT][i][j];
    const int bidx = L.bidx[i][j];
    const uint64_t src = c.blk[bidx];
    if (e.src == src) [[likely]] return e;
    const int g  = L.bMod[bidx];
    const uint64_t base = (i < j) ? mulOff(OFFM.ofBlock[bidx], c.curT,       src)
                                  : mulOff(OFFM.ofBlock[bidx], NEGU[c.curT], src);
    uint64_t best = base, am = 1;
    for (int r = 1; r < g; ++r) {
        const uint64_t v = rotr(base, r, g);
        if      (v >  best) { best = v; am = uint64_t(1) << r; }
        else if (v == best) { am |= uint64_t(1) << r; }
    }
    e.src = src; e.base = base; e.maxv = best; e.argmask = am;
    return e;
}

// Column B of the image word: choose rho(B) and d_B, compare blocks
// (0,B),(1,B),...,(B-1,B),(B,B).  Returns false as soon as some image word is
// lexicographically larger than the target (then X is not canonical).
template<int B>
[[gnu::always_inline]] static inline bool rec(Canon& c, Bitset<8> used) {
    if constexpr (B >= K) {
        return true;                                  // unreachable: the LB cut-off fires first
    } else {
        // The still free orbits of the right length, as a mask: no per-i branch.
        for (int i : (Bitset<8>{CANDM[B]} & ~used).tracked()) {
            Bitset<64> surv{};                        // candidate values of d_B; set by a = 0
            bool done = false;
            for (int a = 0; a < B; ++a) {
                const int bidx = L.bidx[a][B];
                const int g    = L.bMod[bidx];
                const uint64_t tgt = c.blk[bidx];
                const Canon::RInfo& ri = getRot(c, c.rho[a], i);
                if (a == 0) {
                    // d_0 = 0 and g = gcd(m_0,m_B) = m_B, so d_B runs over exactly the
                    // g rotations: the largest one decides, and the argmax set is the
                    // set of surviving d_B.  O(1) on a cache hit.
                    if (ri.maxv > tgt) [[unlikely]] return false;  // image > target
                    if (ri.maxv < tgt) [[likely]] { done = true; break; }
                    surv = Bitset<64>{ri.argmask};
                } else {
                    if (ri.maxv < tgt) [[likely]] { done = true; break; }  // no rotation reaches tgt
                    const int dag = MODT[g][c.dsh[a]];
                    Bitset<64> ns{};
                    for (int d : surv.tracked()) {
                        int r = MODT[g][d] - dag; if (r < 0) r += g;
                        const uint64_t v = rotr(ri.base, r, g);
                        if (v > tgt) [[unlikely]] return false;  // image > target
                        if (v == tgt) ns.addBit(d);   // otherwise this d is pruned
                    }
                    surv = ns;
                    if (surv.isEmpty()) { done = true; break; }
                }
                if (bidx == c.LB)   { done = true; break; }   // exact match, nothing left
            }
            if (done) continue;
            const int bdiag = L.bidx[B][B];
            const uint64_t v = getDiag(c, i), tgt = c.blk[bdiag];
            if (v > tgt) [[unlikely]] return false;
            if (v < tgt) [[likely]] continue;
            if (bdiag == c.LB) continue;
            c.rho[B] = i;
            const Bitset<8> used2 = used.with(i);
            for (int d : surv.tracked()) {
                c.dsh[B] = d;
                if (!rec<B + 1>(c, used2)) return false;
            }
        }
        return true;
    }
}

// Is the current set X the lexicographically largest word in its N-orbit?
// isCanonical and the levels rec<B> are inlined into the search loop (a few per
// cent; the compiler does that by itself only when there is a single caller).
[[gnu::always_inline]] static inline bool isCanonical(Canon& c) {
    const uint64_t tgt = c.blk[0];                    // block (0,0) is the first block
    Canon::L0Info& e = c.l0;
    // Level 0 compares t*S_{i0 i0} with block (0,0) for every unit t and every
    // candidate i0 = rho(0); it therefore depends on block (0,0) and on the
    // diagonal blocks of the orbits of maximal length (LOI[0] = 0).
    bool fresh = true;
    for (int j = 0; fresh && j <= HII[0]; ++j) fresh = (e.src[j] == c.blk[L.bidx[j][j]]);
    if (!fresh) {
        for (int j = 0; j <= HII[0]; ++j) e.src[j] = c.blk[L.bidx[j][j]];
        e.reject = 0; e.n = 0;
        for (int ti = 0; ti < NU && !e.reject; ++ti) {
            c.curT = ti;
            for (int i0 = LOI[0]; i0 <= HII[0]; ++i0) {
                const uint64_t v = getDiag(c, i0);
                if (v > tgt) { e.reject = 1; break; } // some t*S_{ii} beats block (0,0)
                if (v == tgt) { e.t[e.n] = uint8_t(ti); e.i0[e.n] = uint8_t(i0); ++e.n; }
            }
        }
    }
    if (e.reject) return false;
    if (c.LB == 0) return true;                       // nothing after block (0,0) to compare
    for (int x = 0; x < e.n; ++x) {
        c.curT = e.t[x];
        const int i0 = e.i0[x];
        c.rho[0] = i0; c.dsh[0] = 0;                  // c_{rho(0)} normalised to 0
        if (!rec<1>(c, Bitset<8>{}.with(i0))) return false;
    }
    return true;
}

// ===========================================================================
// 4. Search
// ===========================================================================

// The depth (number of edge orbits) at which the tree is cut into tasks for the
// worker threads (doc/DESIGN.md 6); set from the command line before the search.
static int splitDepth = 7;

struct Node {
    uint64_t blk[NB];
    int maxPos;      // -1 when X is empty
    int lastBlk;
    int nOrb;
    int nEdges;
    int deg[K];      // degree of a vertex of orbit i (constant on orbits)
};

// A counter that exists only in some builds (Empty otherwise), so that a build
// without the corresponding filter carries neither the counter nor its updates.
[[maybe_unused]] static inline void bump(uint64_t& a) { ++a; }
[[maybe_unused]] static inline void bump(Empty&) {}
[[maybe_unused]] static inline void addTo(uint64_t& a, uint64_t b) { a += b; }
[[maybe_unused]] static inline void addTo(Empty&, const Empty&) {}
[[maybe_unused]] static inline uint64_t cval(uint64_t a) { return a; }
[[maybe_unused]] static inline uint64_t cval(const Empty&) { return 0; }
// The same for the chromatic counter, whose absent form is a type of its own so that
// it overlaps with the other absent member of Counters (two Empty's could not).
struct ChiOff {};
struct MuteOff {};                                  // the mute flag of a build without PARTS
[[maybe_unused]] static inline bool muted(const MuteOff&) { return false; }
[[maybe_unused]] static inline bool muted(bool m) { return m; }
[[maybe_unused]] static inline void setMute(MuteOff&, bool) {}
[[maybe_unused]] static inline void setMute(bool& m, bool v) { m = v; }
[[maybe_unused]] static inline void bump(ChiOff&) {}
[[maybe_unused]] static inline void addTo(ChiOff&, const ChiOff&) {}
[[maybe_unused]] static inline uint64_t cval(const ChiOff&) { return 0; }

// alignas(64): the workers' Counters lie side by side in one std::vector, and
// without it a thread's total/nodes would share a cache line with the tail of
// its neighbour's (false sharing, which costs a factor 2 with many threads).
struct alignas(64) Counters {
    uint64_t total = 0;       // graphs reported (satisfying the forbidden-subgraph constraints)
    uint64_t nodes = 0;       // nodes of the (pruned) search tree
    uint64_t testsE = 0, cutE = 0;   // battery calls (2g) on the decided edges, and how many said "yes"
    uint64_t testsN = 0, cutN = 0;   // the same on the decided non-edges (a "yes" cuts the loop, or, in
                                     // the propagation of extendFC, forces the orbit the other way)
    uint64_t forcedE = 0, forcedN = 0, contra = 0;   // forward checking: forced orbits, dead nodes
    std::array<uint64_t, NPOSA + 1> byOrb{};
    std::vector<uint64_t> byEdges;
    [[no_unique_address]] std::conditional_t<HAS_CHI, uint64_t, ChiOff> cutC{};   // candidates / children cut by the chromatic bound (7.5)
    [[no_unique_address]] std::conditional_t<G6, G6Writer, G6Off> g6;             // GRAPH6: this thread's graph6 lines
    [[no_unique_address]] std::conditional_t<PART_MODE, bool, MuteOff> mute{};     // PARTS: report nothing (parts r > 0 above the split depth)
    Counters() : byEdges(size_t(L.totalEdges) + 1, 0) {}
    void add(const Counters& s) {
        total += s.total; nodes += s.nodes;
        testsE += s.testsE; cutE += s.cutE; testsN += s.testsN; cutN += s.cutN;
        forcedE += s.forcedE; forcedN += s.forcedN; contra += s.contra; addTo(cutC, s.cutC);
        for (int i = 0; i <= NPOS; ++i) byOrb[i] += s.byOrb[i];
        for (size_t i = 0; i < byEdges.size(); ++i) byEdges[i] += s.byEdges[i];
    }
};

// HAS_DEG: are the degrees of G_X within the bounds?  (The upper bound holds by
// construction -- a candidate that would pass it is never tried -- so only the
// lower bound can fail here.)
[[maybe_unused]] static inline bool degOK(const Node& nd) {
    for (int i = 0; i < K; ++i) if (nd.deg[i] < DMIN || (DMAX > 0 && nd.deg[i] > DMAX)) return false;
    return true;
}

static inline void report(const Node& nd, Counters& ct) {
    if (muted(ct.mute)) return;                          // PARTS: parts r > 0 report nothing above the split depth
    if constexpr (HAS_DEG) { if (!degOK(nd)) return; }
    ++ct.total; ++ct.byOrb[nd.nOrb]; ++ct.byEdges[nd.nEdges];
    g6Write(ct.g6, nd.blk);
}

// ============================ PREFILTER SLOT ===============================
// The child tests run right after the edge orbit p has been added to X and X has
// been found canonical: the degree limit below, the forbidden-subgraph battery of
// 2g, then the colouring of 2e; a failure cuts the whole subtree.  Tests that are
// cheaper than the canonicity test run before it, on the candidate (degMayAdd,
// chiMayAdd).  A new property goes into the same place (doc/DESIGN.md 7.1).
//
// All pairs in edge orbits <= p are decided at this point: those in X are edges
// of every graph reported below this node, the others are non-edges of all of
// them.  So every property closed under edge deletion (degree bounds,
// H-freeness) may be tested here on G_X, and every property closed under edge
// addition on the graph of decided non-edges.
//
// Available in nd (all maintained incrementally, no allocation):
//   nd.blk[]   the decided connection sets, block by block (reversed bit order)
//   nd.deg[i]  degree of a vertex of orbit i in G_X
//   nd.nEdges  number of edges of G_X          nd.nOrb  number of edge orbits
//   nd.maxPos  = p    nd.lastBlk = block of p  (L.bB[nd.lastBlk] = current column;
//              at the end of a column all pairs inside O_0..O_col are decided)
// ===========================================================================

// The degree state of a node (nd holds its degrees; the undecided positions are
// those > nd.maxPos): the last candidate position of its loop, written at its
// depth a[nd.nOrb].  A candidate is only worth trying while the positions from it
// on can still lift every degree to MINDEG.
[[maybe_unused]] static inline void degLimit(Empty&, const Node&) {}                // no degree bounds
template<class A>
static inline void degLimit(A& a, const Node& nd) {
    int Q = NPOS - 1;
    if constexpr (DMIN > 0) {
        for (int i = 0; i < K; ++i) if (nd.deg[i] < DMIN) Q = std::min(Q, DT.lastP[i][DMIN - nd.deg[i]]);
    }
    a[nd.nOrb].degQ = int16_t(Q);
}
// The loop limit of the node nd, from its DegState.
[[maybe_unused]] static inline int  degLast(const Empty&, const Node&) { return NPOS - 1; }
[[maybe_unused]] static inline int  degLast(const DegArr& a, const Node& nd) { return a[nd.nOrb].degQ; }
// May the position p (orbits a, b) be tried as an edge of a child of nd?  (The
// upper bound; degrees only grow.)
static inline bool degMayAdd(const Node& nd, int p) {
    if constexpr (DMAX > 0) {
        if (nd.deg[L.posA[p]] + L.posDegA[p] > DMAX || nd.deg[L.posB[p]] + L.posDegB[p] > DMAX) return false;
    }
    (void)nd; (void)p; return true;
}
// MAXCHI=2: does X + {p} stay bipartite?  Decided before the canonicity test; the
// child's state is written at its depth on the way (and simply overwritten by the
// next candidate when the answer is no or the child is not canonical).
[[maybe_unused]] static inline bool chiMayAdd(BipOff&, const Node&, int, Counters&) { return true; }
[[maybe_unused]] static inline bool chiMayAdd(BipArr& a, const Node& nd, int p, Counters& ct) {
    a[nd.nOrb + 1] = a[nd.nOrb];
    if (bipApply(a[nd.nOrb + 1], p)) return true;
    bump(ct.cutC); return false;
}

// ---------------------------------------------------------------------------
// Forward checking (-DFORWARD=1).  Before the candidate loop of a node, every
// undecided orbit q > max X is tried both ways against the current graphs (the
// battery of 2g on each side): if q as a non-edge creates a forbidden subgraph in
// N it is a *forced edge* of every valid graph below the node, if q as an edge
// creates one in E it is a *forced non-edge*; both at once is a contradiction and
// the node is dead.  Forced orbits are added to the graphs (they are in every
// valid completion, so every cut derived from them is sound) and the scan
// repeats until nothing changes.  The graphs are inherited by the children
// (forcedness is monotone), so a child only propagates what its own decision
// adds.  In the loop, a forced non-edge is never tried as a candidate (no
// canonicity test), and the first forced edge is the last candidate: no later
// orbit can be the first child on the path to a valid graph.
//
// nSame: the decided non-edges are exactly what they were when the parent's
// propagation reached its fixpoint, so "q as a non-edge is forbidden in N" is
// already known to be false for every still-undecided q; the first pass then skips
// that (expensive) side and only tests it on demand, when the edge side says
// "forced non-edge" (the contradiction check must stay exact, and so must the
// invariant that adjN holds no forbidden subgraph).  Later passes test both sides
// in full.
// ---------------------------------------------------------------------------
// aligned(64), here and on extend(): these are the hot loops, and where they start
// relative to a cache line changes the speed by a few percent (measured: 3-4 % on
// 5,5,5,5 and 8,8,8 from a change elsewhere that moved them by 16 bytes).  Aligned,
// their placement no longer depends on the size of the code before them.
template<bool COLLECT>
[[gnu::aligned(64)]] static void extendFC(Node& nd, Canon& c, Counters& ct, std::vector<Node>* tasks, bool nSame) {
    Forbid& f = c.fb;
    const int mark = f.logN;
    ++ct.nodes;
    const int base = nd.maxPos + 1;
    int last = NPOS - 1;                                   // last candidate of the loop
    // --- propagation ---------------------------------------------------------
    bool skipN = nSame;
    for (bool changed = true; changed; skipN = false) {
        changed = false;
        for (int q = base; q < NPOS; ++q) {
            if (f.state[q]) continue;
            bool fe = false, fn = false;                       // forced edge / forced non-edge
            if constexpr (E_SUB) { fn = wouldForbidE(f.adjE, q, false); ++ct.testsE; ct.cutE += fn; }
            if constexpr (HAS_N) { if (!skipN || fn) { fe = wouldForbidN(f.adjN, q, false); ++ct.testsN; ct.cutN += fe; } }
            if (fe && fn) { ++ct.contra; f.undo(mark); return; }      // dead node: no report, no children
            if (fe) { f.pushE(q); ++ct.forcedE; changed = true; }     // fn false: adjE + q is still clean
            if (fn) { f.pushN(q); ++ct.forcedN; if constexpr (HAS_N) changed = true; }   // fe false (and tested): adjN + q is
        }                                                             // still clean; a new non-edge matters to the N side only
    }
    if constexpr (DMIN > 0) {
        for (int i = 0; i < K; ++i) if (nd.deg[i] < DMIN) last = std::min(last, DT.lastP[i][DMIN - nd.deg[i]]);
    }
    if constexpr (!HAS_N) report(nd, ct);
    // --- candidate loop ------------------------------------------------------
    bool fin = true, nGrew = false;
    for (int p = base; p <= last; ++p) {
        const int st = f.state[p];
        if (st != 2 && degMayAdd(nd, p) && chiMayAdd(c.bp, nd, p, ct)) {   // not a forced non-edge: try it as the next edge
            const int bi = L.posBlk[p];
            nd.blk[bi] |= L.posBit[p];
            c.blk = nd.blk; c.LB = bi;
            if (isCanonical(c)) {
                const int sm = nd.maxPos, sl = nd.lastBlk;
                nd.maxPos = p; nd.lastBlk = bi; ++nd.nOrb; nd.nEdges += L.posEdges[p];
                nd.deg[L.posA[p]] += L.posDegA[p]; nd.deg[L.posB[p]] += L.posDegB[p];
                bool ok = true;
                if constexpr (ADJ_E) {
                    if (st != 1) f.pushE(p);          // no battery call here: at the fixpoint no undecided
                }                                     // orbit creates a forbidden subgraph, and p is one
                if constexpr (CHI_GEN) {              // p is the last log entry here
                    ok = chiCheck<CHI>(c.ch, nd.nOrb, f.adjE);
                    if (!ok) { bump(ct.cutC); if (st != 1) f.undo(f.logN - 1); }
                }
                if (ok) {
                    if constexpr (COLLECT) {
                        if (nd.nOrb >= splitDepth) tasks->push_back(nd);
                        else                       extendFC<true>(nd, c, ct, tasks, !nGrew);
                    } else {
                        extendFC<false>(nd, c, ct, tasks, !nGrew);
                    }
                }
                if constexpr (ADJ_E) { if (st != 1 && ok) f.undo(f.logN - 1); }
                nd.deg[L.posA[p]] -= L.posDegA[p]; nd.deg[L.posB[p]] -= L.posDegB[p];
                nd.nEdges -= L.posEdges[p]; --nd.nOrb; nd.maxPos = sm; nd.lastBlk = sl;
            }
            nd.blk[bi] ^= L.posBit[p];
        }
        if constexpr (HAS_N) {
            if (st == 1) { fin = false; break; }             // forced edge: last candidate, X not final
            if (st != 2) {                                   // p becomes a decided non-edge
                ++ct.testsN;
                if (wouldForbidN(f.adjN, p, true)) { ++ct.cutN; fin = false; break; }
                f.noteN(p); nGrew = true;
            }
        }
    }
    f.undo(mark);
    if constexpr (HAS_N) { if (fin) report(nd, ct); }
}

// The search loop of one node X.  Positions p > max X are tried in order; with
// the forbidden-subgraph filters every position that has been passed becomes a
// decided non-edge for the later candidates AND for X itself (all positions after
// max X are non-edges of G_X).  So:
//   * X + {p} is cut when the orbit p creates a forbidden subgraph on the decided
//     edges (the E-side battery of 2g);
//   * when the non-edge orbit p creates one on the decided non-edges, every
//     later candidate is cut as well (the loop stops) and X is not reported;
//   * X is reported iff the loop ran to the end -- the same walk that prunes the
//     candidates decides whether the complement of G_X is free of them.
// With degree bounds the loop stops at the node's degQ (later candidates cannot
// reach MINDEG; then X cannot either, so the shortened walk never loses a
// report), and a candidate that would pass MAXDEG is not tried.
// Without a filter on the non-edge side every node is reported.
template<bool COLLECT>
[[gnu::aligned(64)]] static void extend(Node& nd, Canon& c, Counters& ct, std::vector<Node>* tasks) {
    if constexpr (!HAS_N) report(nd, ct);
    ++ct.nodes;
    const int base = nd.maxPos + 1;
    const int last = degLast(c.dg, nd);      // HAS_DEG: later candidates cannot reach MINDEG
    int  pend = last + 1;                     // HAS_N: positions base..pend-1 were added to adjN
    bool fin  = true;                         // HAS_N: no non-edge orbit created a forbidden subgraph
    for (int p = base; p <= last; ++p) {
        const int bi = L.posBlk[p];
        if (degMayAdd(nd, p) && chiMayAdd(c.bp, nd, p, ct)) {
            nd.blk[bi] |= L.posBit[p];
            c.blk = nd.blk; c.LB = bi;
            if (isCanonical(c)) {
                const int sm = nd.maxPos, sl = nd.lastBlk;
                nd.maxPos = p; nd.lastBlk = bi; ++nd.nOrb; nd.nEdges += L.posEdges[p];
                nd.deg[L.posA[p]] += L.posDegA[p]; nd.deg[L.posB[p]] += L.posDegB[p];
                bool ok = true;
                if constexpr (HAS_DEG) degLimit(c.dg, nd);
                if constexpr (ADJ_E) {
                    if constexpr (E_SUB) {
                        ++ct.testsE;
                        ok = !wouldForbidE(c.fb.adjE, p, true);      // leaves p in adjE iff ok
                        if (!ok) ++ct.cutE;
                    } else addOrbit(c.fb.adjE, p);                   // no filter, only the colouring: p in adjE
                }
                if constexpr (CHI_GEN) {
                    if (ok) { ok = chiCheck<CHI>(c.ch, nd.nOrb, c.fb.adjE); if (!ok) { bump(ct.cutC); removeOrbit(c.fb.adjE, p); } }
                }
                if (ok) {
                    if constexpr (COLLECT) {
                        if (nd.nOrb >= splitDepth) tasks->push_back(nd);
                        else                       extend<true>(nd, c, ct, tasks);
                    } else {
                        extend<false>(nd, c, ct, tasks);
                    }
                }
                if constexpr (ADJ_E) { if (ok) removeOrbit(c.fb.adjE, p); }
                nd.deg[L.posA[p]] -= L.posDegA[p]; nd.deg[L.posB[p]] -= L.posDegB[p];
                nd.nEdges -= L.posEdges[p]; --nd.nOrb; nd.maxPos = sm; nd.lastBlk = sl;
            }
            nd.blk[bi] ^= L.posBit[p];
        }
        if constexpr (HAS_N) {
            ++ct.testsN;
            if (wouldForbidN(c.fb.adjN, p, true)) { ++ct.cutN; pend = p; fin = false; break; }   // p not added
        }
    }
    if constexpr (HAS_N) {
        for (int q = base; q < pend; ++q) removeOrbit(c.fb.adjN, q);
        if (fin) report(nd, ct);
    }
}

} // namespace poly

// ===========================================================================
// 5. Main
// ===========================================================================

// The summary goes to stdout, unless stdout carries the graphs (GRAPH6).
#if GRAPH6
#define PRINT(...) fprintf(stderr, __VA_ARGS__)
#else
#define PRINT(...) printf(__VA_ARGS__)
#endif

#define VOLTA_STR2(...) #__VA_ARGS__
#define VOLTA_STR(...) VOLTA_STR2(__VA_ARGS__)

static void usage(const char* prog, FILE* out) {
    fprintf(out,
        "usage: %s [-v] [-t THREADS] [-s SPLITDEPTH]%s\n"
        "\n"
        "Generates the k-polycirculant graphs of one parameter set, which is fixed\n"
        "at compile time (this binary: ORBITS=%s).  Normally run through the\n"
        "`volta` script, which compiles and caches such binaries.\n"
        "\n"
        "  -v             verbose summary: the parameters and statistics of the search,\n"
        "                 and the graphs counted by number of edge orbits and of edges\n"
        "  -t THREADS     worker threads (default 1; 0 = one per available CPU)\n"
        "  -s SPLITDEPTH  depth at which the search tree is cut into tasks for the\n"
        "                 threads (default 7); the depth of a node is its number of\n"
        "                 edge orbits, i.e. of orbits of phi on the edges\n"
        "%s"
        "  -h             this text\n"
        "  -V             version\n",
        prog, poly::PART_MODE ? " [-p R/M]" : "", VOLTA_STR(ORBITS),
        poly::PART_MODE ? "  -p R/M         run only part R of M (0 <= R < M); the M parts together\n"
                          "                 report every graph exactly once\n" : "");
}

// The summary.  By default: what was generated, under which conditions, how many
// graphs, and the time; -v adds the parameters of the search, its statistics and
// the counts by number of edge orbits and of edges.  Cold and out of line: it runs
// once.
[[gnu::cold, gnu::noinline]]
static void summary(const poly::Counters& master, bool verbose, int usedThreads, size_t nTasks,
                    [[maybe_unused]] long partRes, [[maybe_unused]] long partMod, double secs) {
    using namespace poly;
    PRINT("volta %s   orbits = ", VOLTA_VERSION);
    for (int i = 0; i < K; ++i) PRINT(i ? ",%d" : "%d", MS[i]);
    PRINT("   n = %d", NV);
    if constexpr (PART_MODE) PRINT("   part %ld/%ld", partRes, partMod);
    PRINT("\n");

    constexpr bool ANY_COND = FC > 0 || FKME > 0 || FCYC > 0 || FCC > 0 || FCKME > 0 || FCCYC > 0
                              || DMIN > 0 || DMAX > 0 || CHI > 0 || GIRTH > 3;
    if constexpr (ANY_COND) {
        const char* sep = "conditions: ";
        auto cond = [&](const char* fmt, int x) { PRINT("%s", sep); PRINT(fmt, x); sep = ", "; };
        if (FC    > 0) cond("no K%d", FC);
        if (FKME  > 0) cond("no K%d-e", FKME);
        if (FCYC  > 0) cond("no C%d", FCYC);
        if (FCC   > 0) cond("no K%d in the complement", FCC);
        if (FCKME > 0) cond("no K%d-e in the complement", FCKME);
        if (FCCYC > 0) cond("no C%d in the complement", FCCYC);
        if (GIRTH > 3) cond("girth >= %d", GIRTH);
        if (CHI   > 0) cond("chromatic number <= %d", CHI);
        if (DMIN > 0 && DMIN == DMAX) cond("%d-regular", DMIN);
        else {
            if (DMIN > 0) cond("minimum degree >= %d", DMIN);
            if (DMAX > 0) cond("maximum degree <= %d", DMAX);
        }
        PRINT("\n");
    }

    if (verbose) {
        PRINT("orbits of phi on the vertex pairs (word length) = %d   units of Z_%d = %d   blocks = %d\n",
              NPOS, M0, NU, NB);
        PRINT("threads = %d   split depth = %d   tasks = %zu\n", usedThreads, splitDepth, nTasks);
        if constexpr (PART_MODE)
            PRINT("part %ld/%ld: the tasks with index = %ld (mod %ld), %zu of %zu; the graphs above the split depth are %s\n",
                  partRes, partMod, partRes, partMod,
                  nTasks > size_t(partRes) ? (nTasks - 1 - size_t(partRes)) / size_t(partMod) + 1 : 0,
                  nTasks, partRes == 0 ? "counted in this part" : "counted in part 0");
#if GRAPH6
        PRINT("graph6: vertex j of orbit i (both from 0, the orbits longest first) is number m_0+...+m_{i-1}+j\n");
#endif
        PRINT("nodes visited = %llu", (unsigned long long)master.nodes);
        if (HAS_CHI) PRINT("   chromatic cuts = %llu", (unsigned long long)cval(master.cutC));
        if (E_SUB) PRINT("   edge checks = %llu (cut %llu)", (unsigned long long)master.testsE, (unsigned long long)master.cutE);
        if (HAS_N) PRINT("   non-edge checks = %llu (cut %llu)", (unsigned long long)master.testsN, (unsigned long long)master.cutN);
        PRINT("\n");
        if (FWD) PRINT("forward checking: forced edge orbits = %llu   forced non-edge orbits = %llu   dead nodes = %llu\n",
                        (unsigned long long)master.forcedE, (unsigned long long)master.forcedN, (unsigned long long)master.contra);
    }

    PRINT("graphs (N-classes) = %llu\n", (unsigned long long)master.total);

    if (verbose && master.total > 0) {
        PRINT("by number of edge orbits (orbits of phi on the edges):\n");
        for (int i = 0; i <= NPOS; ++i) if (master.byOrb[i])
            PRINT("  %3d : %llu\n", i, (unsigned long long)master.byOrb[i]);
        PRINT("by number of edges:\n");
        for (size_t i = 0; i < master.byEdges.size(); ++i) if (master.byEdges[i])
            PRINT("  %3zu : %llu\n", i, (unsigned long long)master.byEdges[i]);
    }

    PRINT("time = %.3f s\n", secs);
}

// A whole decimal number in [lo, hi], or false.
static bool parseNum(const char* s, long lo, long hi, long& v) {
    if (!s || !*s) return false;
    char* end = nullptr;
    errno = 0;
    const long x = strtol(s, &end, 10);
    if (errno || *end || x < lo || x > hi) return false;
    v = x;
    return true;
}

int main(int argc, char** argv) {
    using namespace poly;
    long nThreads = 1, split = 7;
    bool verbose = false;
    // -DPARTS=1: -p r/m (0 <= r < m) selects the tasks idx = r (mod m) of the task
    // vector; the graphs above the split depth are reported by part 0 only, so the m
    // parts together report every graph exactly once (doc/DESIGN.md 6).
    [[maybe_unused]] long partRes = 0, partMod = 1;
    for (int i = 1; i < argc; ++i) {
        const char* a = argv[i];
        const char* val = (i + 1 < argc) ? argv[i + 1] : nullptr;
        if (!strcmp(a, "-h") || !strcmp(a, "--help")) { usage(argv[0], stdout); return 0; }
        if (!strcmp(a, "-V") || !strcmp(a, "--version")) { printf("volta %s\n", VOLTA_VERSION); return 0; }
        if (!strcmp(a, "-v") || !strcmp(a, "--verbose")) { verbose = true; continue; }
        if (!strcmp(a, "-t")) {
            if (!parseNum(val, 0, 4096, nThreads)) { fprintf(stderr, "%s: -t expects a number of threads (0 = all CPUs)\n", argv[0]); return 2; }
            ++i;
        } else if (!strcmp(a, "-s")) {
            if (!parseNum(val, 1, 1000000, split)) { fprintf(stderr, "%s: -s expects a positive split depth\n", argv[0]); return 2; }
            ++i;
        } else if (!strcmp(a, "-p")) {
            if (!PART_MODE) { fprintf(stderr, "%s: this binary was not built for parts (volta --part R/M builds one)\n", argv[0]); return 2; }
            char extra;
            if (!val || sscanf(val, "%ld/%ld%c", &partRes, &partMod, &extra) != 2 || partMod < 1 || partRes < 0 || partRes >= partMod) {
                fprintf(stderr, "%s: -p expects R/M with 0 <= R < M, e.g. -p 3/64\n", argv[0]); return 2;
            }
            ++i;
        } else {
            fprintf(stderr, "%s: unknown argument '%s'\n\n", argv[0], a);
            usage(argv[0], stderr);
            return 2;
        }
    }
    if (nThreads == 0) nThreads = std::max(1u, std::thread::hardware_concurrency());
    splitDepth = int(split);

    buildTables();
    if constexpr (ADJ_E || HAS_N || G6) buildEdgeLists();
    if constexpr (G6) buildG6Orbits();

    Node root{};
    root.maxPos = -1; root.lastBlk = 0; root.nOrb = 0; root.nEdges = 0;

    const auto t0 = std::chrono::steady_clock::now();

    Counters master;
    std::vector<Node> tasks;
    int usedThreads = 1;

    if (nThreads <= 1 && !PART_MODE) {
        Canon c{};
        bipRebuild(c.bp, root.blk, root.maxPos, root.nOrb); chiRebuild(c.ch, 0, c.fb.adjE);
        if constexpr (FWD) extendFC<false>(root, c, master, nullptr, false);
        else { degLimit(c.dg, root); extend<false>(root, c, master, nullptr); }
    } else {
        { Canon c{};
          setMute(master.mute, partRes != 0);                   // PARTS: only part 0 reports above the split depth
          bipRebuild(c.bp, root.blk, root.maxPos, root.nOrb); chiRebuild(c.ch, 0, c.fb.adjE);
          if constexpr (FWD) extendFC<true>(root, c, master, &tasks, false);
          else { degLimit(c.dg, root); extend<true>(root, c, master, &tasks); }
          setMute(master.mute, false); }
        usedThreads = (int)std::min<size_t>(size_t(nThreads), std::max<size_t>(tasks.size(), 1));
        std::atomic<size_t> next{0};
        std::vector<Counters> sub(usedThreads);
        std::vector<std::thread> th;
        for (int w = 0; w < usedThreads; ++w) {
            th.emplace_back([&, w]() {
                Canon c{};
                for (;;) {
                    const size_t idx = next.fetch_add(1, std::memory_order_relaxed);
                    if (idx >= tasks.size()) break;
                    if constexpr (PART_MODE) { if (long(idx % size_t(partMod)) != partRes) continue; }
                    Node nd = tasks[idx];
                    if constexpr (ADJ_E || HAS_N) c.fb.rebuild(nd.blk, nd.maxPos);
                    bipRebuild(c.bp, nd.blk, nd.maxPos, nd.nOrb); chiRebuild(c.ch, nd.nOrb, c.fb.adjE);
                    if constexpr (FWD) extendFC<false>(nd, c, sub[w], nullptr, false);
                    else { degLimit(c.dg, nd); extend<false>(nd, c, sub[w], nullptr); }
                }
                g6Flush(sub[w].g6);
            });
        }
        for (auto& t : th) t.join();
        for (auto& s : sub) master.add(s);
    }
    g6Flush(master.g6);
    if (fflush(stdout) != 0 || ferror(stdout)) { fprintf(stderr, "volta: error writing the graphs: %s\n", strerror(errno)); return 1; }

    const auto t1 = std::chrono::steady_clock::now();
    const double secs = std::chrono::duration<double>(t1 - t0).count();

    summary(master, verbose, usedThreads, tasks.size(), partRes, partMod, secs);
    return 0;
}
