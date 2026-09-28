// VHTNP amplitude cache for POWHEG-BOX-V2/ggxy_ggZH (gg -> ZH at NLO, ggxy + Recola).
//
// Why it exists.  POWHEG's rwl reweighting (pwhg_main.f: rwl_compute_new_weight ->
// setrandom -> gen_btilderw) recomputes B-tilde for every weight in the list with the
// random state restored, i.e. at bitwise-identical kinematics.  Per weight only
// pdf_ndns1/2, st_facfact and st_renfact change (rwl_setup_param_weights.f).  In
// ggxy_ggZH the Born and the real emission are Recola one-loop amplitudes
// (calc_matsq_rec, with a quad-precision rescue) and the virtual is the ggxy two-loop
// amplitude; none of them depends on the PDF or on mu_F, and at fixed order each is an
// exact power of alpha_s.  With the CMS 884-weight list the amplitudes were therefore
// evaluated 884 times per event for nothing - 2 min/event on the 200-event validation.
//
// What it does.  A process-wide hash map from (flavours, momenta, top mass[, mu_R])
// to the alpha_s-stripped squared amplitude.  The key is the raw bit pattern of the
// doubles, so only exactly repeated kinematics hit; a hit reproduces the uncached
// result to the last bit apart from the alpha_s power, which the Fortran caller
// re-applies with the current st_alpha.
//
// Controls (environment):  GGXY_AMPCACHE=0 disables the cache (the binary then behaves
// exactly as without the patch, which is how the A/B test is done);  GGXY_AMPCACHE_MAX
// caps the number of entries (default 400000, ~100 MB); when full the table is cleared.
// A one-line hit/miss summary goes to stderr at exit.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

struct AmpCache {
    std::unordered_map<std::string, double> table;
    std::unordered_map<std::string, std::vector<double>> vtable;   // Born value + spin-correlated tensor
    unsigned long hits = 0, misses = 0, stores = 0, clears = 0;
    std::size_t cap = 400000;
    bool enabled = true;

    AmpCache() {
        const char* e = std::getenv("GGXY_AMPCACHE");
        if (e && e[0] == '0') enabled = false;
        const char* c = std::getenv("GGXY_AMPCACHE_MAX");
        if (c) {
            unsigned long v = std::strtoul(c, nullptr, 10);
            if (v > 0) cap = v;
        }
        if (enabled) table.reserve(1 << 17);
    }
    ~AmpCache() {
        std::fprintf(stderr,
                     "ggxy ampcache: enabled=%d hits=%lu misses=%lu stores=%lu clears=%lu size=%zu vsize=%zu\n",
                     enabled ? 1 : 0, hits, misses, stores, clears, table.size(), vtable.size());
    }
    static std::string key(int nd, const double* d, int ni, const int* i) {
        std::string k(static_cast<std::size_t>(nd) * sizeof(double) +
                          static_cast<std::size_t>(ni) * sizeof(int),
                      '\0');
        std::memcpy(&k[0], d, static_cast<std::size_t>(nd) * sizeof(double));
        std::memcpy(&k[static_cast<std::size_t>(nd) * sizeof(double)], i,
                    static_cast<std::size_t>(ni) * sizeof(int));
        return k;
    }
};

AmpCache& cache() {
    static AmpCache c;
    return c;
}

}  // namespace

// gfortran name mangling: lower case + trailing underscore, all arguments by reference.
extern "C" {

void ampcache_lookup_(const int* nd, const double* d, const int* ni, const int* i,
                      double* val, int* found) {
    AmpCache& c = cache();
    *found = 0;
    if (!c.enabled) return;
    auto it = c.table.find(AmpCache::key(*nd, d, *ni, i));
    if (it == c.table.end()) {
        ++c.misses;
        return;
    }
    *val = it->second;
    *found = 1;
    ++c.hits;
}

void ampcache_store_(const int* nd, const double* d, const int* ni, const int* i,
                     const double* val) {
    AmpCache& c = cache();
    if (!c.enabled) return;
    if (c.table.size() >= c.cap) {
        c.table.clear();
        ++c.clears;
    }
    c.table[AmpCache::key(*nd, d, *ni, i)] = *val;
    ++c.stores;
}

void ampcache_lookupv_(const int* nd, const double* d, const int* ni, const int* i,
                       const int* nv, double* vals, int* found) {
    AmpCache& c = cache();
    *found = 0;
    if (!c.enabled) return;
    auto it = c.vtable.find(AmpCache::key(*nd, d, *ni, i));
    if (it == c.vtable.end() || static_cast<int>(it->second.size()) != *nv) {
        ++c.misses;
        return;
    }
    std::memcpy(vals, it->second.data(), static_cast<std::size_t>(*nv) * sizeof(double));
    *found = 1;
    ++c.hits;
}

void ampcache_storev_(const int* nd, const double* d, const int* ni, const int* i,
                      const int* nv, const double* vals) {
    AmpCache& c = cache();
    if (!c.enabled) return;
    if (c.vtable.size() >= c.cap) {
        c.vtable.clear();
        ++c.clears;
    }
    c.vtable[AmpCache::key(*nd, d, *ni, i)] = std::vector<double>(vals, vals + *nv);
    ++c.stores;
}

void ampcache_stats_(long* hits, long* misses, long* stores) {
    AmpCache& c = cache();
    *hits = static_cast<long>(c.hits);
    *misses = static_cast<long>(c.misses);
    *stores = static_cast<long>(c.stores);
}

}  // extern "C"
