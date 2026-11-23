// ---------------------------------------------------------------------------
// bitset.h -- a fixed-size bitset for volta
//
// Copyright (C) 2023-2026 Steven Van Overberghe
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
// Scalable bitset: supports maxSize ∈ {8, 16, 32} for compact storage,
// and any positive multiple of 64 (64, 128, 192, 256, ...) for larger sets.
// Single-word (≤64) operations compile to bare scalar instructions.
// Multi-word operations are constexpr loops, fully unrolled at -O2.
//

#ifndef CODING_BITSET_H
#define CODING_BITSET_H

#include <cstdint>
#include <array>
#include <vector>
#include <iostream>

// ---------------------------------------------------------------------------
// Forward declarations
// ---------------------------------------------------------------------------
template<int maxSize>
struct BitsetIterator;

template<int maxSize>
struct TrackedBitsetIterator;

template<int maxSize>
struct TrackedBitsetRange;

// ---------------------------------------------------------------------------
// Bitset
// ---------------------------------------------------------------------------
template<int maxSize = 64>
struct Bitset {

    static_assert(maxSize == 8 || maxSize == 16 || maxSize == 32 ||
                  (maxSize >= 64 && maxSize % 64 == 0),
                  "maxSize must be 8, 16, 32, or a positive multiple of 64");

    // --- Type traits -------------------------------------------------------

    static constexpr bool MultiWord = (maxSize > 64);
    static constexpr int  NumWords  = MultiWord ? (maxSize / 64) : 1;

    /// Scalar type used for single-word / compact storage.
    using scalar_t = std::conditional_t<maxSize <= 8,  uint8_t,
                     std::conditional_t<maxSize <= 16, uint16_t,
                     std::conditional_t<maxSize <= 32, uint32_t, uint64_t>>>;

    /// Actual storage: a native integer for ≤64 bits, an array otherwise.
    using storage_t = std::conditional_t<MultiWord,
                                         std::array<uint64_t, NumWords>,
                                         scalar_t>;

    /// Backward-compatible alias (used by e.g. pipeFormat).
    using bitword = storage_t;

    // --- Data --------------------------------------------------------------

    storage_t w;

    // === Static helpers ====================================================

    /// Return a Bitset with only bit i set.
    constexpr static Bitset bitAt(int i) {
        if constexpr (MultiWord) {
            Bitset r{};
            r.w[i / 64] = uint64_t{1} << (i % 64);
            return r;
        } else {
            return Bitset{scalar_t(scalar_t{1} << i)};
        }
    }

    /// Return a Bitset with all maxSize bits set.
    constexpr static Bitset allOnes() {
        if constexpr (MultiWord) {
            Bitset r{};
            for (int i = 0; i < NumWords; ++i)
                r.w[i] = ~uint64_t{0};
            return r;
        } else {
            return Bitset{scalar_t(~scalar_t{0})};
        }
    }

    /// Construct a bitset of which the first 'length' bits are set.
    constexpr static Bitset bitRange(int length) {
        if (length == 0)       return Bitset{};
        if (length == maxSize) return allOnes();
        if constexpr (MultiWord) {
            Bitset r{};
            int full = length / 64;
            int rem  = length % 64;
            for (int i = 0; i < full; ++i)
                r.w[i] = ~uint64_t{0};
            if (rem > 0)
                r.w[full] = (uint64_t{1} << rem) - 1;
            return r;
        } else {
            return Bitset{scalar_t((scalar_t{1} << length) - 1)};
        }
    }

    /// Construct a bitset of which bits [startIncl, stopExcl) are set.
    constexpr static Bitset bitRange(int startIncl, int stopExcl) {
        return bitRange(stopExcl - startIncl) << startIncl;
    }

    /// Construct a bitset for which all bits ≥ startIncl are set.
    constexpr static Bitset fromRange(int startIncl) {
        return allOnes() << startIncl;
    }

    // === Accessors =========================================================

    constexpr bool hasBit(int i) const {
        if constexpr (MultiWord) {
            return (w[i / 64] & (uint64_t{1} << (i % 64))) != 0;
        } else {
            return (w & (scalar_t{1} << i)) != scalar_t{0};
        }
    }

    constexpr bool operator[](int i) const { return hasBit(i); }

    // === Mutators ==========================================================

    void addBit(int i) {
        if constexpr (MultiWord) {
            w[i / 64] |= uint64_t{1} << (i % 64);
        } else {
            w |= scalar_t{1} << i;
        }
    }

    void add(int i) { addBit(i); }

    void removeBit(int i) {
        if constexpr (MultiWord) {
            w[i / 64] &= ~(uint64_t{1} << (i % 64));
        } else {
            w &= ~(scalar_t{1} << i);
        }
    }

    void remove(int i) { removeBit(i); }

    void toggleBit(int i) {
        if constexpr (MultiWord) {
            w[i / 64] ^= uint64_t{1} << (i % 64);
        } else {
            w ^= scalar_t{1} << i;
        }
    }

    /// Return a copy with bit i cleared.
    constexpr Bitset without(int i) const {
        if constexpr (MultiWord) {
            Bitset r = *this;
            r.w[i / 64] &= ~(uint64_t{1} << (i % 64));
            return r;
        } else {
            return Bitset{scalar_t(w & ~(scalar_t{1} << i))};
        }
    }

    /// Return a copy with bit i set.
    constexpr Bitset with(int i) const {
        if constexpr (MultiWord) {
            Bitset r = *this;
            r.w[i / 64] |= uint64_t{1} << (i % 64);
            return r;
        } else {
            return Bitset{scalar_t(w | (scalar_t{1} << i))};
        }
    }

    // === Weight queries ====================================================

    constexpr bool isEmpty() const {
        if constexpr (MultiWord) {
            for (int i = 0; i < NumWords; ++i)
                if (w[i] != 0) return false;
            return true;
        } else {
            return w == scalar_t{0};
        }
    }

    constexpr bool isNotEmpty() const { return !isEmpty(); }

    constexpr int popcount() const {
        if constexpr (MultiWord) {
            int c = 0;
            for (int i = 0; i < NumWords; ++i)
                c += __builtin_popcountll(w[i]);
            return c;
        } else if constexpr (maxSize <= 32) {
            return __builtin_popcount(static_cast<unsigned>(w));
        } else {
            return __builtin_popcountll(w);
        }
    }

    constexpr int size() const { return popcount(); }

    /// Check if the number of set bits equals a runtime value.
    constexpr bool hasWeight(int weight) const {
        return popcount() == weight;
    }

    /// Check if the number of set bits equals a compile-time value.
    template<int weight>
    constexpr bool hasWeight() const {
        static_assert(weight >= 0);
        if constexpr (weight == 0) {
            return isEmpty();
        } else if constexpr (weight == 1) {
            return isNotEmpty() && withoutLowest().isEmpty();
        } else {
            return hasWeight(weight);
        }
    }

    /// Check if the number of set bits is at least a compile-time value.
    template<int weight>
    constexpr bool hasAtLeastWeight() const {
        static_assert(weight >= 0);
        if constexpr (weight <= 0) {
            return true;
        } else if constexpr (weight == 1) {
            return isNotEmpty();
        } else if constexpr (weight == 2) {
            return withoutLowest().isNotEmpty();
        } else {
            return popcount() >= weight;
        }
    }

    // === Bit scanning ======================================================

    /// Index of the first (lowest) set bit.  Undefined if empty.
    constexpr int findFirstSetBit() const {
        if constexpr (MultiWord) {
            for (int i = 0; i < NumWords; ++i)
                if (w[i] != 0)
                    return i * 64 + __builtin_ctzll(w[i]);
            __builtin_unreachable();
        } else if constexpr (maxSize <= 32) {
            return __builtin_ctz(static_cast<unsigned>(w));
        } else {
            return __builtin_ctzll(w);
        }
    }

    /// Return this bitset with the lowest set bit cleared.
    constexpr Bitset withoutLowest() const {
        if constexpr (MultiWord) {
            Bitset r = *this;
            for (int i = 0; i < NumWords; ++i) {
                if (r.w[i] != 0) {
                    r.w[i] &= r.w[i] - 1;
                    return r;
                }
            }
            return r;  // was empty
        } else {
            return Bitset{scalar_t(w & (w - 1))};
        }
    }

    /// Return only the lowest set bit (as a Bitset).
    constexpr Bitset lowestOrderBit() const {
        return (*this) ^ withoutLowest();
    }

    // === Complement ========================================================

    /// Return the bits in [0..n) not present in this bitset.
    constexpr Bitset complement(int n) const {
        return bitRange(n) & (~(*this));
    }

    /// If you know this set is a subset of [0..n).
    constexpr Bitset complementByXOR(int n) const {
        return bitRange(n) ^ (*this);
    }

    // === Iteration =========================================================

    /// Find first set bit, clear it, write its index to 'index'.
    /// Returns whether the bitset was non-empty.
    template<bool safe = true>
    bool consume_next_bit(int& index) {
        bool cont = isNotEmpty();
        if constexpr (safe) {
            if (!cont) return false;
        }
        index = findFirstSetBit();
        *this = withoutLowest();
        return cont;
    }

    constexpr BitsetIterator<maxSize> begin() const;
    constexpr static BitsetIterator<maxSize> end();

    /// Return a range adapter that uses a word-tracking iterator.
    /// For multi-word Bitsets, this avoids rescanning exhausted words.
    /// For single-word Bitsets, this is identical to the default iterator.
    constexpr TrackedBitsetRange<maxSize> tracked() const;

    // === Constructors ======================================================

    constexpr Bitset() : w{} {}

    /// Construct from raw storage (scalar or array).
    constexpr Bitset(const storage_t& word) : w(word) {}

    /// Multi-word: allow construction from a single uint64_t (sets lowest word).
    template<bool M = MultiWord, typename = std::enable_if_t<M>>
    constexpr Bitset(uint64_t val) : w{} { w[0] = val; }

    /// Construct from a vector of bit indices.
    explicit Bitset(const std::vector<int>& indices) : w{} {
        for (int i : indices) addBit(i);
    }

    /// Construct from a std::array of bit indices (constexpr).
    template<typename T, unsigned long size>
    constexpr explicit Bitset(const std::array<T, size>& indices) : w{} {
        for (unsigned long i = 0; i < size; ++i)
            setBitConstexpr(indices[i]);
    }

    // === Bitwise operators =================================================

    constexpr Bitset operator&(const Bitset& o) const {
        if constexpr (MultiWord) {
            Bitset r{};
            for (int i = 0; i < NumWords; ++i)
                r.w[i] = w[i] & o.w[i];
            return r;
        } else {
            return Bitset{scalar_t(w & o.w)};
        }
    }

    constexpr void operator&=(const Bitset& o) {
        if constexpr (MultiWord) {
            for (int i = 0; i < NumWords; ++i) w[i] &= o.w[i];
        } else {
            w &= o.w;
        }
    }

    constexpr Bitset operator|(const Bitset& o) const {
        if constexpr (MultiWord) {
            Bitset r{};
            for (int i = 0; i < NumWords; ++i)
                r.w[i] = w[i] | o.w[i];
            return r;
        } else {
            return Bitset{scalar_t(w | o.w)};
        }
    }

    constexpr void operator|=(const Bitset& o) {
        if constexpr (MultiWord) {
            for (int i = 0; i < NumWords; ++i) w[i] |= o.w[i];
        } else {
            w |= o.w;
        }
    }

    constexpr Bitset operator^(const Bitset& o) const {
        if constexpr (MultiWord) {
            Bitset r{};
            for (int i = 0; i < NumWords; ++i)
                r.w[i] = w[i] ^ o.w[i];
            return r;
        } else {
            return Bitset{scalar_t(w ^ o.w)};
        }
    }

    constexpr void operator^=(const Bitset& o) {
        if constexpr (MultiWord) {
            for (int i = 0; i < NumWords; ++i) w[i] ^= o.w[i];
        } else {
            w ^= o.w;
        }
    }

    constexpr Bitset operator~() const {
        if constexpr (MultiWord) {
            Bitset r{};
            for (int i = 0; i < NumWords; ++i)
                r.w[i] = ~w[i];
            return r;
        } else {
            return Bitset{scalar_t(~w)};
        }
    }

    constexpr Bitset operator<<(int x) const {
        if constexpr (MultiWord) {
            if (x <= 0) return *this;
            if (x >= maxSize) return Bitset{};
            Bitset r{};
            int wordShift = x / 64;
            int bitShift  = x % 64;
            if (bitShift == 0) {
                for (int i = NumWords - 1; i >= wordShift; --i)
                    r.w[i] = w[i - wordShift];
            } else {
                for (int i = NumWords - 1; i > wordShift; --i)
                    r.w[i] = (w[i - wordShift] << bitShift)
                           | (w[i - wordShift - 1] >> (64 - bitShift));
                r.w[wordShift] = w[0] << bitShift;
            }
            return r;
        } else {
            return Bitset{scalar_t(w << x)};
        }
    }

    constexpr Bitset operator>>(int x) const {
        if constexpr (MultiWord) {
            if (x <= 0) return *this;
            if (x >= maxSize) return Bitset{};
            Bitset r{};
            int wordShift = x / 64;
            int bitShift  = x % 64;
            int limit = NumWords - wordShift;
            if (bitShift == 0) {
                for (int i = 0; i < limit; ++i)
                    r.w[i] = w[i + wordShift];
            } else {
                for (int i = 0; i < limit - 1; ++i)
                    r.w[i] = (w[i + wordShift] >> bitShift)
                           | (w[i + wordShift + 1] << (64 - bitShift));
                r.w[limit - 1] = w[NumWords - 1] >> bitShift;
            }
            return r;
        } else {
            return Bitset{scalar_t(w >> x)};
        }
    }

    // === Comparison ========================================================

    constexpr bool operator==(const Bitset& o) const {
        if constexpr (MultiWord) {
            for (int i = 0; i < NumWords; ++i)
                if (w[i] != o.w[i]) return false;
            return true;
        } else {
            return w == o.w;
        }
    }

    constexpr bool operator!=(const Bitset& o) const {
        return !(*this == o);
    }

    // === Output ============================================================

    friend std::ostream& operator<<(std::ostream& os, const Bitset& set) {
        for (int i = 0; i < maxSize; ++i)
            os << (set[i] ? '1' : '0');
        return os;
    }

private:
    /// Helper for constexpr array-constructor (no if-constexpr in member init).
    constexpr void setBitConstexpr(int i) {
        if constexpr (MultiWord) {
            w[i / 64] |= uint64_t{1} << (i % 64);
        } else {
            w |= scalar_t{1} << i;
        }
    }
};

// ---------------------------------------------------------------------------
// BitsetIterator
// ---------------------------------------------------------------------------
template<int maxSize = 64>
struct BitsetIterator {
    Bitset<maxSize> current;

    constexpr bool operator!=(const BitsetIterator& other) const {
        return current != other.current;
    }

    void operator++() {
        current = current.withoutLowest();
    }

    int operator*() const {
        return current.findFirstSetBit();
    }
};

template<int maxSize>
constexpr BitsetIterator<maxSize> Bitset<maxSize>::begin() const {
    return BitsetIterator<maxSize>{*this};
}

template<int maxSize>
constexpr BitsetIterator<maxSize> Bitset<maxSize>::end() {
    return BitsetIterator<maxSize>{Bitset<maxSize>{}};
}

// ---------------------------------------------------------------------------
// TrackedBitsetIterator
// ---------------------------------------------------------------------------
// For multi-word Bitsets, tracks which word the current lowest bit lives in,
// so operator++ / operator* avoid scanning from word 0 each time.
// For single-word Bitsets, the word index is a zero-size empty struct,
// making this identical to the basic BitsetIterator.
// ---------------------------------------------------------------------------

namespace bitset_detail {
    struct EmptyTag {};
}

template<int maxSize>
struct TrackedBitsetIterator {
    using Bs = Bitset<maxSize>;
    static constexpr bool MW = Bs::MultiWord;
    static constexpr int  NW = Bs::NumWords;

    Bs current;

    [[no_unique_address]]
    std::conditional_t<MW, int, bitset_detail::EmptyTag> wordIdx_;

    /// Construct from a Bitset (begin position).
    constexpr explicit TrackedBitsetIterator(const Bs& bs) : current(bs), wordIdx_{} {
        if constexpr (MW) {
            wordIdx_ = 0;
            while (wordIdx_ < NW && current.w[wordIdx_] == 0)
                ++wordIdx_;
        }
    }

    /// Construct an end sentinel.
    constexpr TrackedBitsetIterator() : current{}, wordIdx_{} {
        if constexpr (MW) wordIdx_ = NW;
    }

    constexpr bool operator!=(const TrackedBitsetIterator& other) const {
        if constexpr (MW) {
            return wordIdx_ != other.wordIdx_;
        } else {
            return current != other.current;
        }
    }

    void operator++() {
        if constexpr (MW) {
            // Clear lowest bit in the active word.
            current.w[wordIdx_] &= current.w[wordIdx_] - 1;
            // Advance to the next non-zero word if this one is exhausted.
            while (wordIdx_ < NW && current.w[wordIdx_] == 0)
                ++wordIdx_;
        } else {
            current = current.withoutLowest();
        }
    }

    int operator*() const {
        if constexpr (MW) {
            return wordIdx_ * 64 + __builtin_ctzll(current.w[wordIdx_]);
        } else {
            return current.findFirstSetBit();
        }
    }
};

// ---------------------------------------------------------------------------
// TrackedBitsetRange — lightweight adapter for range-for
// ---------------------------------------------------------------------------
template<int maxSize>
struct TrackedBitsetRange {
    Bitset<maxSize> bs;

    constexpr TrackedBitsetIterator<maxSize> begin() const {
        return TrackedBitsetIterator<maxSize>{bs};
    }

    constexpr TrackedBitsetIterator<maxSize> end() const {
        return TrackedBitsetIterator<maxSize>{};
    }
};

template<int maxSize>
constexpr TrackedBitsetRange<maxSize> Bitset<maxSize>::tracked() const {
    return TrackedBitsetRange<maxSize>{*this};
}

#endif //CODING_BITSET_H
