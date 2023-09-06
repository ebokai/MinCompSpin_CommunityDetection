#include "header.h"

// map<__uint128_t, unsigned int> build_pdata(map<__uint128_t, unsigned int> &data, __uint128_t community) {

map<__uint128_t, unsigned int> build_pdata(vector<pair<__uint128_t, unsigned int>> &data, __uint128_t community) {

    map<__uint128_t, unsigned int> pdata;
    // map<__uint128_t, unsigned int>::iterator it;

    // unsigned int ks;
    // __uint128_t state, mask_state;
    __uint128_t mask_state;

    // for (it = data.begin(); it != data.end(); it++) {

    //     state = it -> first;
    //     ks = it -> second;

    //     mask_state = state & community;
    //     pdata[mask_state] += ks;
    // }

    for (auto const &it : data) {
        mask_state = ((it).first) & community;
        pdata[mask_state] += ((it).second);
    }

    return pdata;
}

double icc_evidence(__uint128_t community, Partition &p_struct){

    double logE = 0;
    double pf;
    unsigned int k;

    map<__uint128_t, unsigned int> pdata = build_pdata(p_struct.data, community);
    map<__uint128_t, unsigned int>::iterator it;

    unsigned int rank = bit_count(community);
    

    if (rank > 32) {
        pf = -((double) rank - 1.) * p_struct.N * log(2);
    } else {
        double rank_pow = (double) (ONE << (rank - 1));
        pf = lgamma(rank_pow) - lgamma(p_struct.N + rank_pow);
    }

    logE += pf;
    for (it = pdata.begin(); it != pdata.end(); it ++){
            k = it -> second;
            logE += lgamma(k + 0.5) - lgamma(0.5);
        }

    return logE;
}