#include <bits/stdc++.h>
#include <assert.h>

#include "rxdp.hpp"

using namespace std;

static inline uint64_t verifyRXDP_fixed_inputs(int n, int k, uint64_t alpha, uint64_t beta) {
    uint64_t bound_n = (1ULL << n);
    uint64_t mask_n = bound_n - 1;
    uint64_t n_tests = 0;
    
    vector<uint64_t> ddt(bound_n);

    // compute the distribution with fixed input differences
    for(uint64_t x = 0; x < bound_n; x++){
        for(uint64_t y = 0; y < bound_n; y++){
            uint64_t value1 = CSHL((x+y) & mask_n,k,n);
            uint64_t value2 = ((CSHL(x,k,n) ^ alpha) + (CSHL(y,k,n) ^ beta)) & mask_n;
            uint64_t delta = value1 ^ value2;

            ddt[delta]++;
        }
    }

    for(uint64_t delta = 0; delta < bound_n; delta++) {
        if (checkCriteria(alpha,beta,delta,n,k) != (ddt[delta] > 0)) {
            printf("verifyRXDP failed validity!\n");
            exit(-1);
        }
        n_tests++;

        // note: getRXDiffCount tests validity inside too
        // so we can call it even if it's impossible
        uint64_t theory_count = getRXDiffCount(alpha,beta,delta,n,k);
        uint64_t exper_count = ddt[delta];
        if (theory_count != exper_count){
            printf("verifyRXDP failed probability!\n");
            exit(-1);
        }
    }
    assert(n_tests == 1ull << n);
    return n_tests;
}

static void verifyRXDP_exhaustive(int n, int k) {
    uint64_t bound_n = (1ULL << n);
    
    #pragma omp parallel for
    for(uint64_t alpha = 0; alpha < bound_n; alpha++){
        for(uint64_t beta = 0; beta < bound_n; beta++){
            verifyRXDP_fixed_inputs(n, k, alpha, beta);
        }
    }
    printf("verifyRXDP_exhaustive n=%d k=%d finished successfully!\n", n, k);
}

static void verifyRXDP_random(int n, int k, uint64_t count) {
    uint64_t bound_n = (1ULL << n);
    uint64_t mask_n = bound_n - 1;
    
    for(uint64_t itr = 0; itr < count; itr++) {
        uint64_t alpha = random() & mask_n;
        uint64_t beta = random() & mask_n;
        verifyRXDP_fixed_inputs(n, k, alpha, beta);
    }
    printf("verifyRXDP_random n=%d k=%d count=%lu finished successfully!\n", n, k, count);
}

static inline bool LipmaaMoriai_validity(int n, uint64_t alpha, uint64_t beta, uint64_t delta, uint64_t bar) {
    assert(alpha < (1ull << n));
    assert(beta < (1ull << n));
    assert(delta < (1ull << n));
    assert(bar < 2);
    uint64_t u = (alpha^beta^delta) ^ shl(alpha^beta^delta,1,n);
    uint64_t v = shl( (alpha^delta) | (beta^delta) ,1,n);

    u ^= bar; // Lemma 1 or Lemma 2 if bar=1
    return ((u&v) == u);
}

static inline uint64_t LipmaaMoriai_count(int n, uint64_t alpha, uint64_t beta, uint64_t delta, uint64_t bar) {
    assert(alpha < (1ull << n));
    assert(beta < (1ull << n));
    assert(delta < (1ull << n));
    assert(bar < 2);
    if (!LipmaaMoriai_validity(n, alpha, beta, delta, bar)) {
        return 0;
    }
    uint64_t mask = (1ULL << n) - 1;
    uint64_t tmp = (((alpha ^ delta) | (beta ^ delta)) << 1)&mask;
    uint d = __builtin_popcountll(tmp);
    return (1ULL << (2*n-d));
}

static inline uint64_t Rn_lemma(int n, uint64_t alpha, uint64_t beta, uint64_t delta, uint64_t bar) {
    assert(alpha < (1ull << n));
    assert(beta < (1ull << n));
    assert(delta < (1ull << n));
    assert(bar < 2);
    if (!LipmaaMoriai_validity(n, alpha, beta, delta, bar)) {
        return 0;
    }
    assert( (1 & (alpha ^ beta ^ delta)) == bar );

    assert((LipmaaMoriai_count(n, alpha, beta, delta, bar) & 1) == 0);
    assert(LipmaaMoriai_count(n, alpha, beta, delta, bar) > 0);

    // all-zeros or all-ones depending on bar
    uint64_t target_xor = bar ? ((1ull << n) - 1) : 0;
    if ((alpha ^ beta ^ delta) != target_xor) {
        return LipmaaMoriai_count(n, alpha, beta, delta, bar) / 2;
    }
    else {
        return LipmaaMoriai_count(n, alpha, beta, delta, bar) / 2 + (1ull << (n - 1));
    }
}

static inline void verifyDiff_fixed_inputs(int n, uint64_t alpha, uint64_t beta) {
    uint64_t bound_n = (1ULL << n);
    uint64_t mask_n = bound_n - 1;

    // bar corresponds to adding 1 (or not) in the first input pair's addition
    for(uint64_t bar = 0; bar < 2; bar++) {
        vector<uint64_t> ddt(bound_n);
        vector<uint64_t> ddt_no_carry(bound_n); // no output carry in the first input pair's addition
        vector<uint64_t> ddt_carry(bound_n); // no output carry in the first input pair's addition

        for(uint64_t x = 0; x < bound_n; x++){
            for(uint64_t y = 0; y < bound_n; y++){
                uint64_t value1 = (x+y) & mask_n;
                uint64_t value2 = ((x^alpha) + (y^beta) + bar) & mask_n;
                uint64_t delta = value1 ^ value2;

                ddt[delta]++;
                if (x + y < bound_n) {
                    ddt_no_carry[delta]++;
                }
                else {
                    ddt_carry[delta]++;
                }
            }
        }

        for(uint64_t delta = 0; delta < bound_n; delta++) {
            bool valid = LipmaaMoriai_validity(n, alpha, beta, delta, bar);
            if (valid != (ddt[delta] > 0)) {
                printf("LipmaaMoriai (Lemma 1/2) failed validity!\n");
                exit(-1);
            }
            // LipmaaMoriai_count and Rn_lemma checks validity so we can check it again
            // if (!valid) {
            //     continue;
            // }
            
            uint64_t theory_XDS = LipmaaMoriai_count(n, alpha, beta, delta, bar);
            uint64_t exper_XDS = ddt[delta];
            if (theory_XDS != exper_XDS){
                printf("LipmaaMoriai (Lemma 1/2) failed probability!\n");
                exit(-1);
            }

            // Tn = Rn if output carry = 0
            uint64_t theory_Tn_carry0 = Rn_lemma(n, alpha, beta, delta, bar);
            uint64_t exper_Tn_carry0 = ddt_no_carry[delta];

            if (theory_Tn_carry0 != exper_Tn_carry0){
                printf("Rn lemma failed probability!\n");
                exit(-1);
            }

            // Tn = complement of Rn if output carry = 1
            uint64_t theory_Tn_carry1 = theory_XDS - Rn_lemma(n, alpha, beta, delta, bar);
            uint64_t exper_Tn_carry1 = ddt_carry[delta];
            if (theory_Tn_carry1 != exper_Tn_carry1){
                printf("Rn lemma failed probability!\n");
                exit(-1);
            }
        }
    }
}

static void verifyDiff_exhaustive(int n) {
    uint64_t bound_n = (1ULL << n);
    
    #pragma omp parallel for
    for(uint64_t alpha = 0; alpha < bound_n; alpha++){
        for(uint64_t beta = 0; beta < bound_n; beta++){
            verifyDiff_fixed_inputs(n, alpha, beta);
        }
    }
    printf("verifyDiff_exhaustive n=%d finished successfully!\n", n);
}

static void verifyDiff_random(int n, uint64_t count) {
    uint64_t bound_n = (1ULL << n);
    uint64_t mask_n = bound_n - 1;
    
    for(uint64_t itr = 0; itr < count; itr++) {
        uint64_t alpha = random() & mask_n;
        uint64_t beta = random() & mask_n;
        // printf("%016llx %016llx\n", alpha, beta);
        verifyDiff_fixed_inputs(n, alpha, beta);
    }
    printf("verifyDiff_random n=%d count=%lu finished successfully!\n", n, count);
}

static inline void verifyTheory_fixed_valid(int n, uint64_t alpha, uint64_t beta, uint64_t delta, uint64_t bar, uint64_t x, uint64_t y) {
    uint64_t bound_n = (1ULL << n);
    uint64_t mask_n = bound_n - 1;

    // Proposition 2 (XDS bar = complement all diffs)
    // note: we substitute bar and ~bar to ensure that we are in the valid case
    uint64_t theory_XDS_complement = LipmaaMoriai_count(n, mask_n & ~alpha, mask_n & ~beta, mask_n & ~delta, bar ^ 1);
    uint64_t theory_XDS_bar = LipmaaMoriai_count(n, alpha, beta, delta, bar);
    assert(theory_XDS_complement > 0);
    assert(theory_XDS_bar > 0);
    if (theory_XDS_complement != theory_XDS_bar){
        printf("Propositon 2 failed probability!\n");
        exit(-1);
    }

    // Lemma 3
    if (n > 1 && bar == 0) {
        uint64_t bound_n1 = mask_n >> 1;
        uint64_t alpha_n1 = alpha & bound_n1;
        uint64_t beta_n1 = beta & bound_n1;
        uint64_t delta_n1 = delta & bound_n1;

        uint64_t theory_XDS_n = LipmaaMoriai_count(n, alpha, beta, delta, 0); // bar = 0
        uint64_t theory_XDS_n1 = LipmaaMoriai_count(n-1, alpha_n1, beta_n1, delta_n1, 0); // bar = 0
        assert(theory_XDS_n > 0);
        assert(theory_XDS_n1 > 0);

        uint64_t a_msb = alpha_n1 >> (n - 2);
        uint64_t b_msb = beta_n1 >> (n - 2);
        uint64_t d_msb = delta_n1 >> (n - 2);
        if (a_msb == b_msb && b_msb == d_msb) {
            if (theory_XDS_n != 4 * theory_XDS_n1) {
                printf("Lemma 3 failed case 1\n");
                exit(-1);
            }
        }
        else {
            if (theory_XDS_n != 2 * theory_XDS_n1) {
                printf("Lemma 3 failed case 2\n");
                exit(-1);
            }
        }
    }

    // Lemma 4
    if (n > 1 && bar == 0) {
        uint64_t mask_n1 = mask_n >> 1;

        uint64_t a_msb = alpha >> (n - 1);
        uint64_t b_msb = beta >> (n - 1);
        uint64_t d_msb = delta >> (n - 1);

        uint64_t carry1 = ((x & mask_n1) + (y & mask_n1)) > mask_n1;
        uint64_t carry2 = (((x ^ alpha) & mask_n1) + ((y ^ beta) & mask_n1)) > mask_n1;
        if ((carry1 ^ carry2) != (a_msb ^ b_msb ^ d_msb)) {
            printf("Lemma 4 failed\n");
            exit(-1);
        }
    }

    // Lemma 5
    if (n > 1 && bar == 0) {
        uint64_t a_msb = alpha >> (n - 1);
        uint64_t b_msb = beta >> (n - 1);
        uint64_t d_msb = delta >> (n - 1);

        uint64_t a_premsb = 1 & (alpha >> (n - 2));
        uint64_t b_premsb = 1 & (beta >> (n - 2));
        uint64_t d_premsb = 1 & (delta >> (n - 2));

        if ((a_premsb == 0 && b_premsb == 0 && d_premsb == 0) && ((a_msb ^ b_msb ^ d_msb) == 1)) {
            printf("Lemma 5 failed case 1\n");
            exit(-1);
        }
        if ((a_premsb == 1 && b_premsb == 1 && d_premsb == 1) && ((a_msb ^ b_msb ^ d_msb) == 0)) {
            printf("Lemma 5 failed case 2\n");
            exit(-1);
        }
    }
}

static void verifyTheory_random(int n, uint64_t count) {
    uint64_t bound_n = (1ULL << n);
    uint64_t mask_n = bound_n - 1;
    
    for(uint64_t itr = 0; itr < count; itr++) {
        uint64_t alpha = random() & mask_n;
        uint64_t beta = random() & mask_n;
        for(uint bar = 0; bar < 2; bar++) {
            uint64_t x = random() & mask_n;
            uint64_t y = random() & mask_n;

            uint64_t delta = ( (x + y + bar) ^ ((x ^ alpha) + (y ^ beta)) ) & mask_n;
            verifyTheory_fixed_valid(n, alpha, beta, delta, bar, x, y);
        }
    }
    printf("verifyTheory_random n=%d count=%lu finished successfully!\n", n, count);
}


static inline void verifyTnC_fixed_valid(int n, uint64_t a, uint64_t b, uint64_t delta, uint64_t Delta) {
    // uint64_t bound_n = (1ULL << n);
    // uint64_t mask_n = bound_n - 1;

    // verify definition match
    for(uint64_t C = 0; C < 2; C++) {
        for(int i = 0; i <= n; i++) {
            uint64_t theory_TnC = RXDPconstTcounter_i(a, b, delta, Delta, C, n, i);
            uint64_t exper_TnC;
            if (i == 0) {
                exper_TnC = 1-C;
            }
            else if (1 <= i && i <= n-1) {
                uint64_t bound_i = (1ULL << i);
                uint64_t mask_i = bound_i - 1;

                exper_TnC = 0;
                for (uint64_t x = 0; x < bound_i; x++) {
                    uint64_t Delta_cand = (((x ^ delta) + a) ^ (x + b)) & mask_i;
                    if (Delta_cand != (Delta & mask_i))
                        continue;
                    uint64_t c_i = ((x & mask_i) + (b & mask_i)) >= bound_i;
                    if (c_i != C)
                        continue;
                    uint64_t c_prime_i = (((x ^ delta) & mask_i) + (a & mask_i)) >= bound_i;
                    uint64_t out_c_diff_i = ((a ^ b ^ delta ^ Delta) >> i) & 1;
                    if ((c_i ^ c_prime_i) != (out_c_diff_i))
                        continue;
                    exper_TnC++;
                }
                if (exper_TnC != theory_TnC) {
                    printf("Theorem 3 failed case 2\n");
                    exit(-1);
                }
            }
        }
    }
}

static void verifyTnC_random(int n, uint64_t count) {
    uint64_t bound_n = (1ULL << n);
    uint64_t mask_n = bound_n - 1;
    
    for(uint64_t itr = 0; itr < count; itr++) {
        uint64_t a = random() & mask_n;
        uint64_t b = random() & mask_n;
        uint64_t delta = random() & mask_n;
        for(uint bar = 0; bar < 2; bar++) {
            uint64_t x = random() & mask_n;

            uint64_t Delta = (((x ^ delta) + a) ^ (x + b)) & mask_n;
            for(uint64_t C = 0; C < 2; C++) {
                verifyTnC_fixed_valid(n, a, b, delta, Delta);
            }
        }
    }
    printf("verifyTnC_random n=%d count=%lu finished successfully!\n", n, count);
}

static inline uint64_t verifyRXDPconst_fixed_inputs(int n, int k, uint64_t a, uint64_t delta) {
    uint64_t bound_n = (1ULL << n);
    uint64_t mask_n = bound_n - 1;
    uint64_t n_tests = 0;
    assert(a < bound_n);
    assert(delta < bound_n);
    
    vector<uint64_t> ddt(bound_n);

    // compute the distribution with fixed input differences
    for(uint64_t x = 0; x < bound_n; x++){
        uint64_t value1 = ((CSHL(x, k, n) ^ delta) + a) & mask_n;
        uint64_t value2 = CSHL((x+a) & mask_n, k, n) & mask_n;
        uint64_t Delta = value1 ^ value2;

        ddt[Delta]++;
    }

    for(uint64_t Delta = 0; Delta < bound_n; Delta++) {
        n_tests++;

        // note: no validity tests, probability computation may return 0
        uint64_t theory_count = getRXDiffConstCount(a,delta,Delta,n,k);
        uint64_t exper_count = ddt[Delta];
        if (theory_count != exper_count){
            printf("verifyRXDPconst failed probability!\n");
            exit(-1);
        }
    }
    assert(n_tests == 1ull << n);
    return n_tests;
}

static void verifyRXDPconst_exhaustive(int n, int k) {
    uint64_t bound_n = (1ULL << n);
    
    #pragma omp parallel for
    for(uint64_t a = 0; a < bound_n; a++){
        for(uint64_t delta = 0; delta < bound_n; delta++){
            verifyRXDPconst_fixed_inputs(n, k, a, delta);
        }
    }
    printf("verifyRXDPconst_exhaustive n=%d k=%d finished successfully!\n", n, k);
}

static void verifyRXDPconst_random(int n, int k, uint64_t count) {
    uint64_t bound_n = (1ULL << n);
    uint64_t mask_n = bound_n - 1;
    
    for(uint64_t itr = 0; itr < count; itr++) {
        uint64_t a = random() & mask_n;
        uint64_t delta = random() & mask_n;
        verifyRXDPconst_fixed_inputs(n, k, a, delta);
    }
    printf("verifyRXDPconst_random n=%d k=%d count=%lu finished successfully!\n", n, k, count);
}

int main(int argc, char *argv[]) {
    srandom(time(NULL));
    srandom(random() ^ getpid());

    if (argc == 1 || atoi(argv[1]) <= 1) {
        printf("Verifying validity criteria and main theorem exhaustively for n <= 8\n");
        for(uint n = 2; n <= 8; n++){ // 8
            for(uint k = 1; k < n; k++){
                verifyRXDP_exhaustive(n, k);
            }
        }
        printf("\n");
    }

    if (argc == 1 || atoi(argv[1]) <= 2) {
        printf("Verifying validity criteria and main theorem for random 50 input differences for n <= 12\n");
        for(uint n = 2; n <= 12; n++){ // 12+
            #pragma omp parallel for
            for(uint k = 1; k < n; k++){
                verifyRXDP_random(n, k, 100);
            }
        }
        printf("\n");
    }

    if (argc == 1 || atoi(argv[1]) <= 3) {
        printf("Verifying XDS/Rn/Tn formulas exhaustively for n <= 8\n");
        for(uint n = 2; n <= 8; n++){ // 8
            verifyDiff_exhaustive(n);
        }
        printf("\n");
    }

    if (argc == 1 || atoi(argv[1]) <= 4) {
        printf("Verifying XDS/Rn/Tn formulas for random 1000 input differences for n <= 12\n");
        for(uint n = 2; n <= 12; n++){ // 12
            verifyDiff_random(n, 100);
        }
        printf("\n");
    }

    if (argc == 1 || atoi(argv[1]) <= 5) {
        // note: number of inputs is 2^2n so we limite n<=30 to fit results in uint64_t
        printf("Verifying XDS/Rn/Tn propositions for random 1 000 000 input differences for n <= 30\n");
        for(uint n = 2; n <= 30; n++){ // 30
            verifyTheory_random(n, 1000000);
        }
        printf("\n");
    }

    if (argc == 1 || atoi(argv[1]) <= 6) {
        printf("Verifying TnC formulas (Theorem 3) for random 1000 constants/differences for n <= 12\n");
        for(uint n = 2; n <= 12; n++){ // 12/10k
            verifyTnC_random(n, 1000);
        }
        printf("\n");
    }

     if (argc == 1 || atoi(argv[1]) <= 7) {
        printf("Verifying RXDP const formula exhaustively for n <= 8\n");
        for(uint n = 2; n <= 8; n++){ // 8
            for(uint k = 1; k < n; k++){
                verifyRXDPconst_exhaustive(n, k);
            }
        }
        printf("\n");
    }

    if (argc == 1 || atoi(argv[1]) <= 8) {
        printf("Verifying RXDP const formula for random 50 input differences for n <= 12\n");
        for(uint n = 2; n <= 12; n++){ // 12+
            #pragma omp parallel for
            for(uint k = 1; k < n; k++){
                verifyRXDPconst_random(n, k, 100);
            }
        }
        printf("\n");
    }


    return 0;
}