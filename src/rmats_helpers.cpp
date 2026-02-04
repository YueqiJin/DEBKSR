#include <Rcpp.h>
using namespace Rcpp;

inline double logit_scalar(double x) {
    if (NumericVector::is_na(x)) return NA_REAL;
    if (x < 0.01) x = 0.01;
    if (x > 0.99) x = 0.99;
    return std::log(x / (1.0 - x));
}

// [[Rcpp::export]]
NumericVector logit_clamped_cpp(NumericVector x) {
    NumericVector out(x.size());
    for (R_xlen_t i = 0; i < x.size(); ++i) {
        double v = x[i];
        if (NumericVector::is_na(v)) {
            out[i] = NA_REAL;
        } else {
            if (v < 0.01) v = 0.01;
            if (v > 0.99) v = 0.99;
            out[i] = std::log(v / (1.0 - v));
        }
    }
    return out;
}

// [[Rcpp::export]]
NumericVector vec2psi_cpp(NumericVector inc, NumericVector skp,
                          double effective_inclusion_length,
                          double effective_skipping_length) {
    R_xlen_t n = inc.size();
    NumericVector psi(n);
    for (R_xlen_t i = 0; i < n; ++i) {
        double incv = inc[i];
        double skpv = skp[i];
        if (NumericVector::is_na(incv) || NumericVector::is_na(skpv) || incv + skpv == 0.0) {
            psi[i] = 0.5;
        } else {
            double inc_adj = incv / effective_inclusion_length;
            double skp_adj = skpv / effective_skipping_length;
            psi[i] = inc_adj / (inc_adj + skp_adj);
        }
    }
    return psi;
}

// [[Rcpp::export]]
double pop_var_cpp(NumericVector x) {
    double sum = 0.0;
    double cnt = 0.0;
    for (R_xlen_t i = 0; i < x.size(); ++i) {
        double v = x[i];
        if (NumericVector::is_na(v)) continue;
        sum += v;
        cnt += 1.0;
    }
    if (cnt == 0.0) return NA_REAL;
    double mean = sum / cnt;
    double acc = 0.0;
    for (R_xlen_t i = 0; i < x.size(); ++i) {
        double v = x[i];
        if (NumericVector::is_na(v)) continue;
        double d = v - mean;
        acc += d * d;
    }
    return acc / cnt;
}

// [[Rcpp::export]]
double myfunc_individual_cpp(double x, double I, double S, double beta, double var,
                             double inclusion_length, double skipping_length) {
    double denom = inclusion_length * x + skipping_length * (1.0 - x);
    double new_psi = inclusion_length * x / denom;
    double logit_x = logit_scalar(x);
    double logit_beta = logit_scalar(beta);
    return -1.0 * (I * std::log(new_psi) + S * std::log(1.0 - new_psi) -
                   (logit_x - logit_beta) * (logit_x - logit_beta) / (2.0 * var) -
                   std::log(x) - std::log(1.0 - x));
}

// [[Rcpp::export]]
double myfunc_individual_der_cpp(double x, double I, double S, double beta, double var,
                                 double inclusion_length, double skipping_length) {
    double denom = inclusion_length * x + skipping_length * (1.0 - x);
    double new_psi = inclusion_length * x / denom;
    double new_psi_der = inclusion_length * skipping_length / (denom * denom);
    double term = logit_scalar(x) - logit_scalar(beta);
    return -1.0 * (I / new_psi * new_psi_der - S / (1.0 - new_psi) * new_psi_der -
                   term / var / x / (1.0 - x) - 1.0 / x + 1.0 / (1.0 - x));
}

struct PsiOptRes {
    double x;
    double f;
};

inline PsiOptRes golden_min(double lower, double upper, double tol, int max_iter,
                            double I, double S, double beta, double var,
                            double inclusion_length, double skipping_length) {
    const double invphi = 0.6180339887498949;      // 1/phi
    const double invphi2 = 1.0 - invphi;           // 1/phi^2
    double a = lower;
    double b = upper;
    double h = b - a;
    if (h <= tol) {
        double mid = 0.5 * (a + b);
        return {mid, myfunc_individual_cpp(mid, I, S, beta, var, inclusion_length, skipping_length)};
    }
    int n = static_cast<int>(std::ceil(std::log(tol / h) / std::log(invphi)));
    double c = a + invphi2 * h;
    double d = a + invphi * h;
    double fc = myfunc_individual_cpp(c, I, S, beta, var, inclusion_length, skipping_length);
    double fd = myfunc_individual_cpp(d, I, S, beta, var, inclusion_length, skipping_length);
    for (int k = 0; k < n && k < max_iter; ++k) {
        if (fc < fd) {
            b = d;
            d = c;
            fd = fc;
            h = invphi * h;
            c = a + invphi2 * h;
            fc = myfunc_individual_cpp(c, I, S, beta, var, inclusion_length, skipping_length);
        } else {
            a = c;
            c = d;
            fc = fd;
            h = invphi * h;
            d = a + invphi * h;
            fd = myfunc_individual_cpp(d, I, S, beta, var, inclusion_length, skipping_length);
        }
        if (std::abs(b - a) < tol) break;
    }
    if (fc < fd) return {c, fc};
    return {d, fd};
}

// [[Rcpp::export]]
List psi_optimize_golden_cpp(NumericVector I, NumericVector S, double beta, double var,
                             double inclusion_length, double skipping_length,
                             double lower = 0.01, double upper = 0.99,
                             int max_iter = 64, double tol = 1e-8) {
    R_xlen_t n = I.size();
    NumericVector psi(n);
    NumericVector value(n);
    for (R_xlen_t idx = 0; idx < n; ++idx) {
        double i = I[idx];
        double s = S[idx];
        PsiOptRes res = golden_min(lower, upper, tol, max_iter, i, s, beta, var, inclusion_length, skipping_length);
        psi[idx] = res.x;
        value[idx] = res.f;
    }
    return List::create(_["psi"] = psi, _["value"] = value);
}
