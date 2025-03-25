#include <cstdint>
#include <iostream>
#include <algorithm> 
#include <set>
#include <vector>
#include <fstream>
#include <string>
#include <assert.h>

#include "rxdp.hpp"

using namespace std;


bool checkCriteria(uint64_t const a,
                   uint64_t const b,
                   uint64_t const d,
                   uint const n,
                   uint const k){

    uint64_t mask = (1ULL << n) - 1;
    uint64_t u = (a^b^d) ^ shl(a^b^d,1,n);
    uint64_t v = shl( (a^d) | (b^d) ,1,n);

    if(k == 0){
        return ((u&v) == u);
    }
    else{
        uint64_t mask1 = (~(1ULL))&mask;
        uint64_t maskgamma = (~(1ULL << k))&mask;
        u = u&mask1&maskgamma;
        v = v&mask1&maskgamma;
        return ((u&v) == u);
    }
}


uint64_t Tcounter(uint64_t const alpha,
                  uint64_t const beta,
                  uint64_t const delta,
                  uint const n,
                  uint const v){

    // cout << "Tcounter " << alpha << " " << beta << " " << delta << " " << n << " " << v << endl;
    uint64_t mask = (1ULL << n) - 1;
    uint64_t tmp = (((alpha ^ delta) | (beta ^ delta)) << 1)&mask;
    uint d = __builtin_popcountll(tmp);
    uint64_t xorall = alpha^beta^delta;

    // cout << "tmp " << tmp << endl;
    // cout << "d " << d << endl;
    // cout << "xorall " << xorall << endl;

    if(xorall == 0 || xorall == mask){
        if(v == 0){ 
            return (1ULL << (2*n-d-1)) + (1ULL << (n-1));
        }
        else{
            return (1ULL << (2*n-d-1)) - (1ULL << (n-1));
        }
    }
    else{
        //Number of solutions = 2^(2n-d-1)
        return (1ULL << (2*n-d-1));
    }
}

uint64_t getRXDiffCount(uint64_t const alpha,
                        uint64_t const beta,
                        uint64_t const delta,
                        uint const n,
                        uint const k){

    if(!checkCriteria(alpha,beta,delta,n,k))
        return 0;

    if(k == 0){
        uint64_t mask = (1ULL << n) - 1;
        uint64_t tmp = (((alpha ^ delta) | (beta ^ delta)) << 1)&mask;
        uint d = __builtin_popcountll(tmp);
        return (1ULL << (2*n-d));
    }
    else{

        uint64_t mask_k = (1ULL << k)-1;
        uint64_t mask_nk = (1ULL << (n-k))-1;

        uint64_t alpha_L = (alpha >> k)&mask_nk;
        uint64_t  beta_L = (beta  >> k)&mask_nk;
        uint64_t delta_L = (delta >> k)&mask_nk;

        uint64_t alpha_R = alpha&mask_k;
        uint64_t  beta_R = beta&mask_k;
        uint64_t delta_R = delta&mask_k;

        uint v  = (alpha_L^beta_L^delta_L)&1;
        uint vp = (alpha_R^beta_R^delta_R)&1;

        uint64_t factornk = Tcounter(alpha_L,beta_L,delta_L,n-k,vp);
        uint64_t factork = Tcounter(alpha_R,beta_R,delta_R,k,v);

        return factornk*factork;
    }
}


const uint64_t MATRIX_RHO[8][2][2] = {
    {{2, 1}, {0, 1}},
    {{0, 0}, {0, 0}},
    {{1, 1}, {0, 0}},
    {{1, 0}, {0, 1}},
    {{2, 0}, {0, 0}},
    {{0, 1}, {0, 1}},
    {{1, 1}, {0, 0}},
    {{1, 0}, {0, 1}},
};


uint64_t RXDPconstTcounter_i(uint64_t a, uint64_t b, uint64_t delta, uint64_t Delta, uint64_t C, int n, int i) {
    if (i == 0) {
        return 1 - C;
    }
    if (i == n) {
        return (
            RXDPconstTcounter_i(a, b, delta, Delta, C, n+1, i) +
            RXDPconstTcounter_i(a, b, delta, Delta ^ (1ull << n), C, n+1, i)
        );
    }
    // 1 <= i <= n-1
    uint64_t prev0 = RXDPconstTcounter_i(a, b, delta, Delta, 0, n, i-1);
    uint64_t prev1 = RXDPconstTcounter_i(a, b, delta, Delta, 1, n, i-1);
    uint64_t bi1 = ((b >> (i-1)) & 1);
    
    if (bi1) {
        swap(prev0, prev1);
    }

    uint64_t out_c_diff = ((a ^ b ^ delta ^ Delta) >> i) & 1;
    uint64_t rho0 = ((a ^ b ^ delta) >> (i-1)) & 1;
    uint64_t rho1 = ((delta ^ Delta) >> (i-1)) & 1;
    uint64_t rho2 = (((a ^ b) >> (i-1)) & 1) ^ out_c_diff;
    uint64_t rho = (rho0 << 2) | (rho1 << 1) | (rho2);
    assert(rho < 8);

    auto mat = MATRIX_RHO[rho];
    uint64_t cur0 = mat[0][0] * prev0 + mat[0][1] * prev1;
    uint64_t cur1 = mat[1][0] * prev0 + mat[1][1] * prev1;

    if (bi1) {
        swap(cur0, cur1);
    }
    return (C == 0) ? cur0 : cur1;
}


uint64_t getRXDiffConstCount(uint64_t a, uint64_t delta, uint64_t Delta, int n, int k) {
    if (k == 0) {
        assert(0 and "not implemented");
    }
    assert(0 < k && k < n);

    if (2*k > n) { // k > n/2
        delta = CSHR(delta, k, n);
        Delta = CSHR(Delta, k, n);
        k = n-k;
        assert (2*k < n);
    }

    uint64_t mask_k = (1ULL << k)-1;
    uint64_t mask_nk = (1ULL << (n-k))-1;
    // uint64_t mask_mid = (1ULL << (n-2*k))-1;

    uint64_t a_L = (a >> (n-k)) & mask_k;
    // uint64_t a_M = (a >> k) & mask_mid;
    uint64_t a_LM = (a >> k) & mask_nk;
    uint64_t a_Rp = a & mask_k;
    uint64_t a_MRp = a & mask_nk;
    
    uint64_t delta_Lp = (delta >> k) & mask_nk;
    uint64_t delta_Rp = delta & mask_k;
    uint64_t Delta_Lp = (Delta >> k) & mask_nk;
    uint64_t Delta_Rp = Delta & mask_k;

    // uint64_t cR = (delta_Lp ^ a_LM ^ a_Rp ^ Delta_Lp) & 1;
    uint64_t cR = (delta_Lp ^ a_LM ^ a_Rp ^ Delta_Lp) & 1;
    uint64_t cL = (delta_Rp ^ a_Rp ^ a_L ^ Delta_Rp) & 1;

    // RXDPconstTcounter_i(a, b, delta, Delta, C, n, i)
    uint64_t factorL = RXDPconstTcounter_i(
        (a_LM + cR) & mask_nk, // a
        a_MRp, // b
        delta_Lp, // delta
        Delta_Lp, // Delta
        cL, // C
        n-k, n-k
    );
    uint64_t factorR = RXDPconstTcounter_i(
        (a_L + cL) & mask_k, // a
        a_Rp, // b
        delta_Rp, // delta
        Delta_Rp, // Delta
        cR, // C
        k, k
    );
    // printf("cntL = %lu; cntR = %lu\n", factorL, factorR);
    return factorL*factorR;
}