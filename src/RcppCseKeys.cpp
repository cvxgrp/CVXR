// Structural-key store for the CSE (common-subexpression elimination) pass in
// Dcp2Cone -- the C++ half of rsrc_tree/reductions/subexpr_cache.R.
//
// WHY THIS IS IN C++ (measured, see notes/session_handoff_2026-08-12_perf_regression.md):
// the pass computes a structural key for EVERY node of EVERY problem, so its
// cost is per-node. The R implementation built those keys as strings -- paste0
// for composite atoms (3.41us) and sprintf("%a") per element for constants
// (13.1us for 50 doubles) -- then used environments, whose keys must be
// character, so integers were converted back to strings to be looked up. On the
// qp bench cell that bookkeeping was 370us per compile, about +6% of total
// compile time, on problems with nothing to share.
//
// R offers no integer-keyed hash map, which is what this actually needs.
// std::unordered_map does, and hashing an integer vector never builds a string
// at all. Measured against the R implementation, same operations:
//
//   composite node, first visit   6.86us -> 1.35us   (5.1x)
//   constant keyed by value (50)  13.10us -> 0.78us  (16.8x)
//   affine-above memo lookup       1.56us -> 0.41us  (3.8x)
//
// (C++ figures are via direct .Call, per D_PERF.6.)
//
// CORRECTNESS -- the invariants this file must never break:
//
//   1. A FALSE MERGE IS A WRONG ANSWER; a missed merge only costs time. Every
//      ambiguous case therefore resolves to "uncacheable".
//   2. NA_INTEGER is INT_MIN, an ordinary value on the wire, and R's
//      `.UNCACHEABLE` IS NA_integer_. A leaked uncacheable child key would
//      otherwise arrive looking like a legitimate component. Any NA in any
//      component makes the whole node uncacheable, and such a node is never
//      memoised. The R side still returns early on `.is_uncacheable()`; this is
//      a backstop, not a substitute.
//   3. Runs are LENGTH-PREFIXED, never sentinel-separated. A separator is only
//      unambiguous if no component can equal it, which nothing guarantees:
//      (shape = c(1,2), children = 3) must not alias (shape = 1,
//      children = c(2,3)) whatever values they hold.
//   4. Doubles are keyed by their exact bit pattern, which is the property
//      sprintf("%a") supplied: equal doubles always collide, unequal ones never
//      do. NaNs are canonicalized before keying (see push_double): NA_real_ is
//      one value and every other NaN is another, matching R and removing the
//      platform-dependence of raw NaN payloads (an earlier version keyed raw
//      bits and CI caught K(NaN) != K(0/0) on x86). These keys are per-apply
//      and must never be persisted or compared across processes.
//   5. unordered_map compares keys EXACTLY after hashing, so a hash collision
//      selects a slower bucket and can never merge unlike subtrees.
//   6. The store is PER-APPLY. `aa_memo` depends on `quad_obj`, so sharing a
//      store across applies with different `quad_obj` would be wrong.

#include <Rcpp.h>
#include <unordered_map>
#include <string>
#include <vector>
#include <cstring>
#include <cstdint>

namespace {

struct VecHash {
  std::size_t operator()(const std::vector<int>& v) const {
    // FNV-1a over the words. Only a bucket selector -- equality is exact.
    std::size_t h = 1469598103934665603ULL;
    for (int x : v) {
      h ^= static_cast<std::size_t>(static_cast<uint32_t>(x));
      h *= 1099511628211ULL;
    }
    return h;
  }
};

struct CseStore {
  std::unordered_map<int, int> expr_memo;                 // expression id -> key
  std::unordered_map<int, int> aa_memo;                   // expression id -> 0/1
  std::unordered_map<std::vector<int>, int, VecHash> sig; // signature -> key
  std::unordered_map<std::string, int> str_codes;         // class name / string -> code
  std::unordered_map<int, Rcpp::RObject> results;         // key -> canonicalized node
  int n_sig = 0;
};

inline CseStore* store_of(SEXP s) {
  Rcpp::XPtr<CseStore> p(s);
  CseStore* st = p.get();
  if (st == nullptr) Rcpp::stop("CSE key store is NULL (was it carried across a session?)");
  return st;
}

// Intern a fully-built signature, optionally memoising it under an expression
// id. `id == NA_INTEGER` means "payload, do not memoise".
inline int intern(CseStore* st, std::vector<int>& v, int id) {
  int k;
  auto f = st->sig.find(v);
  if (f != st->sig.end()) {
    k = f->second;
  } else {
    k = st->n_sig++;
    st->sig.emplace(std::move(v), k);
  }
  if (id != NA_INTEGER) st->expr_memo[id] = k;
  return k;
}

// Fold a double into the signature as its exact 64 bits, low word first.
// NaNs are canonicalized FIRST: raw NaN bit patterns are platform-dependent
// (0/0 is the negative quiet NaN on x86, the positive one on ARM), so keying
// raw bits made key equality platform-dependent -- caught by CI on Linux and
// Windows, where K(NaN) != K(0/0) while both R and the legacy sprintf("%a")
// key (R normalizes NaN formatting to the string "NaN") consider them equal.
// R's semantics are the spec: NA_real_ is distinct from every other NaN
// (R_IsNA tests its payload), and all other NaNs are ONE value.
inline void push_double(std::vector<int>& v, double x) {
  double canon = x;
  if (ISNAN(x)) canon = R_IsNA(x) ? NA_REAL : R_NaN;
  uint64_t bits;
  std::memcpy(&bits, &canon, sizeof(bits));
  v.push_back(static_cast<int>(bits & 0xffffffffULL));
  v.push_back(static_cast<int>(bits >> 32));
}

inline bool has_na(const Rcpp::IntegerVector& v) {
  for (int x : v) if (x == NA_INTEGER) return true;
  return false;
}

}  // namespace

// Create a per-apply key store.
//
// @return an external pointer (Rcpp::XPtr) to a CseStore instance.
// [[Rcpp::export(.CseKeys__new)]]
SEXP CseKeys__new() {
  return Rcpp::XPtr<CseStore>(new CseStore, true);
}

// Intern a string (a class name, or a character payload element) to a code.
// [[Rcpp::export(.CseKeys__class_code)]]
int CseKeys__class_code(SEXP s, std::string cls) {
  CseStore* st = store_of(s);
  auto it = st->str_codes.find(cls);
  if (it != st->str_codes.end()) return it->second;
  int k = static_cast<int>(st->str_codes.size());
  st->str_codes.emplace(std::move(cls), k);
  return k;
}

// Memo lookup by expression id. NA_INTEGER means MISS -- uncacheable nodes are
// never memoised, so a hit is always a valid key (>= 0).
// [[Rcpp::export(.CseKeys__memo_get)]]
int CseKeys__memo_get(SEXP s, int id) {
  CseStore* st = store_of(s);
  auto it = st->expr_memo.find(id);
  return it == st->expr_memo.end() ? NA_INTEGER : it->second;
}

// The hot path: memo check, signature build, intern and memo store, in ONE
// crossing. Composite atoms key by (class, shape, get_data() payload, child
// keys) -- subexpr_cache.py lines 81-134.
//
// Types are NOT coerced here: callers use .Call directly (D_PERF.6), so `shape`
// and `child_keys` must already be INTEGER vectors and `id`, `class_code`,
// `data_key` plain integers.
// [[Rcpp::export(.CseKeys__key_node)]]
int CseKeys__key_node(SEXP s, int id, int class_code,
                      Rcpp::IntegerVector shape, int data_key,
                      Rcpp::IntegerVector child_keys) {
  if (id == NA_INTEGER || class_code == NA_INTEGER || data_key == NA_INTEGER)
    return NA_INTEGER;
  if (has_na(shape) || has_na(child_keys)) return NA_INTEGER;

  CseStore* st = store_of(s);
  auto hit = st->expr_memo.find(id);
  if (hit != st->expr_memo.end()) return hit->second;

  std::vector<int> v;
  v.reserve(shape.size() + child_keys.size() + 4);
  v.push_back(class_code);
  v.push_back(data_key);
  v.push_back(static_cast<int>(shape.size()));
  for (int d : shape) v.push_back(d);
  v.push_back(static_cast<int>(child_keys.size()));
  for (int c : child_keys) v.push_back(c);
  return intern(st, v, id);
}

// Integer-vector signature: leaves (tag + id), logical payloads, and composed
// list payloads. `id = NA_INTEGER` -> intern without memoising.
// [[Rcpp::export(.CseKeys__intern_ints)]]
int CseKeys__intern_ints(SEXP s, int id, int tag, Rcpp::IntegerVector x) {
  if (tag == NA_INTEGER || has_na(x)) return NA_INTEGER;
  CseStore* st = store_of(s);
  std::vector<int> v;
  v.reserve(x.size() + 2);
  v.push_back(tag);
  v.push_back(static_cast<int>(x.size()));
  for (int e : x) v.push_back(e);
  return intern(st, v, id);
}

// Numeric values keyed by exact bits, with their dims so a 1x3 and a 3x1
// holding the same numbers never key alike.
// [[Rcpp::export(.CseKeys__intern_doubles)]]
int CseKeys__intern_doubles(SEXP s, int id, int tag,
                            Rcpp::NumericVector x, Rcpp::IntegerVector dims) {
  if (tag == NA_INTEGER || has_na(dims)) return NA_INTEGER;
  CseStore* st = store_of(s);
  std::vector<int> v;
  v.reserve(2 * x.size() + dims.size() + 3);
  v.push_back(tag);
  v.push_back(static_cast<int>(dims.size()));
  for (int d : dims) v.push_back(d);
  v.push_back(static_cast<int>(x.size()));
  for (double e : x) push_double(v, e);
  return intern(st, v, id);
}

// Sparse values keyed from sparse storage, never densified. `i`, `j` and `x`
// must already be ordered canonically by the caller.
// [[Rcpp::export(.CseKeys__intern_sparse)]]
int CseKeys__intern_sparse(SEXP s, int id, int tag,
                           Rcpp::IntegerVector i, Rcpp::IntegerVector j,
                           Rcpp::NumericVector x, Rcpp::IntegerVector dims) {
  if (tag == NA_INTEGER || has_na(dims) || has_na(i) || has_na(j))
    return NA_INTEGER;
  if (i.size() != j.size() || i.size() != x.size())
    Rcpp::stop("CSE sparse key: i, j and x must have equal length");
  CseStore* st = store_of(s);
  std::vector<int> v;
  v.reserve(2 * i.size() + 2 * x.size() + dims.size() + 4);
  v.push_back(tag);
  v.push_back(static_cast<int>(dims.size()));
  for (int d : dims) v.push_back(d);
  v.push_back(static_cast<int>(i.size()));
  for (R_xlen_t n = 0; n < i.size(); ++n) { v.push_back(i[n]); v.push_back(j[n]); }
  for (double e : x) push_double(v, e);
  return intern(st, v, id);
}

// Character payloads. NA_STRING makes the payload uncacheable: a character NA
// is legitimate data and there is no code that could represent it without
// risking collision with the literal string "NA".
// [[Rcpp::export(.CseKeys__intern_strings)]]
int CseKeys__intern_strings(SEXP s, int id, int tag, Rcpp::CharacterVector x) {
  if (tag == NA_INTEGER) return NA_INTEGER;
  CseStore* st = store_of(s);
  std::vector<int> v;
  v.reserve(x.size() + 2);
  v.push_back(tag);
  v.push_back(static_cast<int>(x.size()));
  for (R_xlen_t n = 0; n < x.size(); ++n) {
    if (Rcpp::CharacterVector::is_na(x[n])) return NA_INTEGER;
    std::string e = Rcpp::as<std::string>(x[n]);
    auto it = st->str_codes.find(e);
    if (it != st->str_codes.end()) {
      v.push_back(it->second);
    } else {
      int c = static_cast<int>(st->str_codes.size());
      st->str_codes.emplace(std::move(e), c);
      v.push_back(c);
    }
  }
  return intern(st, v, id);
}

// Introspection for tests. The R implementation exposed its maps as
// environments, so a test could count them with ls(); these keep that
// observable available now that the maps are C++ containers.
// [[Rcpp::export(.CseKeys__memo_size)]]
int CseKeys__memo_size(SEXP s) {
  return static_cast<int>(store_of(s)->expr_memo.size());
}

// [[Rcpp::export(.CseKeys__n_signatures)]]
int CseKeys__n_signatures(SEXP s) {
  return store_of(s)->n_sig;
}

// affine_above relevance memo (dcp2cone.py:275-302). NA_INTEGER = miss,
// 0 = FALSE, 1 = TRUE.
// [[Rcpp::export(.CseKeys__aa_get)]]
int CseKeys__aa_get(SEXP s, int id) {
  CseStore* st = store_of(s);
  auto it = st->aa_memo.find(id);
  return it == st->aa_memo.end() ? NA_INTEGER : it->second;
}

// [[Rcpp::export(.CseKeys__aa_set)]]
void CseKeys__aa_set(SEXP s, int id, int val) {
  if (id == NA_INTEGER || val == NA_INTEGER) return;
  store_of(s)->aa_memo[id] = val;
}

// Canonicalized-result cache, keyed by the integer cache key. Absent -> NULL.
// A stored result is always an Expression, never NULL, so NULL is unambiguous.
// [[Rcpp::export(.CseKeys__result_get)]]
SEXP CseKeys__result_get(SEXP s, int key) {
  if (key == NA_INTEGER) return R_NilValue;
  CseStore* st = store_of(s);
  auto it = st->results.find(key);
  return it == st->results.end() ? R_NilValue : static_cast<SEXP>(it->second);
}

// [[Rcpp::export(.CseKeys__result_set)]]
void CseKeys__result_set(SEXP s, int key, SEXP value) {
  if (key == NA_INTEGER) return;
  store_of(s)->results[key] = Rcpp::RObject(value);
}
