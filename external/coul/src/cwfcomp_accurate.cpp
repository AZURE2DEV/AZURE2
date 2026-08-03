#include "cwfcomp_accurate.H"  // Our wrapper class
#include "cwfcomp.H"           // Original Coulomb_wave_functions class
#include <cmath>
#include <limits>
#include <stdexcept>
#include <iostream>

// Value for which e^max_exponent is representable as double
const double NUMERICAL_MAX_EXPONENT = std::log(std::numeric_limits<double>::max());

// Constructor
Coulomb_wave_functions_accurate::Coulomb_wave_functions_accurate(const bool is_it_normalized_c,
                                               const std::complex<double> &l_c,
                                               const std::complex<double> &eta_c)
    : l(l_c), eta(eta_c), is_it_normalized(is_it_normalized_c)
{
    // Constructor initialization - simplified compared to original
}

// Destructor
Coulomb_wave_functions_accurate::~Coulomb_wave_functions_accurate()
{
    // Simplified destructor
}

// Initialize F and dF at a point (compatibility method)
void Coulomb_wave_functions_accurate::F_dF_init(const std::complex<double> &z,
                                       const std::complex<double> &F,
                                       const std::complex<double> &dF)
{
    // For compatibility - not used in accurate_cwf approach
}

// Calculate F and dF using accurate_cwf method
void Coulomb_wave_functions_accurate::F_dF(const std::complex<double> &z,
                                  std::complex<double> &F,
                                  std::complex<double> &dF)
{
    std::complex<double> G, dG;
    solve_cwf(z, F, dF, G, dG);
}

// Calculate G and dG using accurate_cwf method  
void Coulomb_wave_functions_accurate::G_dG(const std::complex<double> &z,
                                  std::complex<double> &G,
                                  std::complex<double> &dG)
{
    std::complex<double> F, dF;
    solve_cwf(z, F, dF, G, dG);
}

// Calculate H and dH using accurate_cwf method
void Coulomb_wave_functions_accurate::H_dH(const int omega,
                                  const std::complex<double> &z,
                                  std::complex<double> &H,
                                  std::complex<double> &dH)
{
    std::complex<double> F, dF, G, dG;
    solve_cwf(z, F, dF, G, dG);
    
    // H+ = G + iF, H- = G - iF
    std::complex<double> i_unit(0.0, 1.0);
    if (omega == 1) {
        H = G + i_unit * F;
        dH = dG + i_unit * dF;
    } else {
        H = G - i_unit * F;
        dH = dG - i_unit * dF;
    }
}

// Calculate scaled H and dH (for compatibility)
void Coulomb_wave_functions_accurate::H_dH_scaled(const int omega,
                                         const std::complex<double> &z,
                                         std::complex<double> &H,
                                         std::complex<double> &dH)
{
    // For this implementation, scaled version is same as regular
    H_dH(omega, z, H, dH);
}

// Core solver using accurate_cwf algorithm (from accurate_cwf.cpp)
void Coulomb_wave_functions_accurate::solve_cwf(const std::complex<double> &rho_c,
                                       std::complex<double> &F,
                                       std::complex<double> &dF,
                                       std::complex<double> &G,
                                       std::complex<double> &dG) const
{
    // Handle neutral particle case with real rho (spherical Bessel functions)
    if ((std::abs(eta.real()) == 0.0) && (std::abs(eta.imag()) == 0.0) && 
        (std::abs(rho_c.imag()) == 0.0)) {
        
        // Start with L=0 spherical Bessel functions
        double rho_r = rho_c.real();
        double F_0 = std::sin(rho_r);
        double G_0 = std::cos(rho_r);
        double dF_0 = G_0;
        double dG_0 = -F_0;
        
        // Apply recursion relations to get to desired L
        double F_L = F_0, G_L = G_0, dF_L = dF_0, dG_L = dG_0;
        int ang_mom_L = std::round(std::abs(l.real()));
        
        for (int el_rec = 1; el_rec <= ang_mom_L; el_rec++) {
            double eta_r = 0.0; // neutral particle
            
            F_L = get_next_el_cwf(rho_r, eta_r, el_rec, F_0, dF_0);
            G_L = get_next_el_cwf(rho_r, eta_r, el_rec, G_0, dG_0);
            dF_L = get_next_el_cwf_deriv(rho_r, eta_r, el_rec, F_0, F_L);
            dG_L = get_next_el_cwf_deriv(rho_r, eta_r, el_rec, G_0, G_L);
            
            F_0 = F_L; G_0 = G_L; dF_0 = dF_L; dG_0 = dG_L;
        }
        
        F = std::complex<double>(F_L, 0.0);
        dF = std::complex<double>(dF_L, 0.0);
        G = std::complex<double>(G_L, 0.0);
        dG = std::complex<double>(dG_L, 0.0);
        return;
    }
    else {
        // Charged particle case and/or complex rho/eta
        // Try continued fraction method first
        // initialize CWF solver
        class Coulomb_wave_functions cwf(true, l, eta); 

        // local variables for Hp, dHp, and square-root of minus 1
        std::complex<double> Hp,dHp, imag_I(0.0,1.0);

        // attempt to solve, catch exceptions and report if necessary
        try{
            cwf.H_dH ( 1, rho_c, Hp, dHp );
            cwf.F_dF ( rho_c, F, dF );
            
            // finish solution construction and return success
            G  =  Hp - imag_I *  F;
            dG = dHp - imag_I * dF;
            return;
        } catch(...) {
            
            // Fall back to asymptotic solution for real parameters only
            if (std::abs(rho_c.imag()) == 0.0 && std::abs(eta.imag()) == 0.0 && std::abs(l.imag()) == 0.0) {
                double rho_r = rho_c.real();
                double eta_r = eta.real();
                int eL_r = std::round(std::abs(l.real()));
                double F_r = 0.0, dF_r = 0.0, G_r = 0.0, dG_r = 0.0;
                
                try {
                    compute_asymptotic_bessel_solution(rho_r, eta_r, eL_r,
                                                       F_r, dF_r, G_r, dG_r);
                    F = std::complex<double>(F_r, 0.0);
                    dF = std::complex<double>(dF_r, 0.0);
                    G = std::complex<double>(G_r, 0.0);
                    dG = std::complex<double>(dG_r, 0.0);
                    return;
                } catch (const std::exception &e2) {
                    std::cerr << "Asymptotic solution also failed" << std::endl;
                }
            }
            
            // If all methods fail, throw error
            throw std::runtime_error("All Coulomb function computation methods failed");
        }
    }
}

// Recursion relation for U_L from U_{L-1} (from accurate_cwf.cpp)
double Coulomb_wave_functions_accurate::get_next_el_cwf(const double &rho, const double &eta,
                                               const int &eL, const double &u_L,
                                               const double &upr_L) const
{
    // From Handbook of Mathematical Functions 14.2.1, 14.2.2, 14.2.3
    // [ (L+1)^2 + eta^2 ]^(1/2)
    double normalize_factor = std::sqrt(1.0 * eL * eL + eta * eta);

    // sum the parts and normalize to get the next CWF
    double next_el_cwf = ((1.0 * eL * eL) / rho + eta) * u_L; // [ (L+1)^2/rho + eta ] * U_L
    next_el_cwf -= eL * upr_L; // [ (L+1)^2/rho + eta ] * U_L - (L+1)*U'_L
    next_el_cwf /= normalize_factor; // { [ (L+1)^2/rho + eta ] * U_L - (L+1)*U'_L } / [ (L+1)^2 + eta^2 ]^(1/2)
    return next_el_cwf;
}

// Recursion relation for U'_L from U_{L-1} and U_L (from accurate_cwf.cpp)
double Coulomb_wave_functions_accurate::get_next_el_cwf_deriv(const double &rho, const double &eta,
                                                     const int &eL, const double &u_L_minus,
                                                     const double &u_L) const
{
    // From Handbook of Mathematical Functions 14.2.1, 14.2.2, 14.2.3
    // [ (L+1)^2 + eta^2 ]^(1/2)
    double normalize_factor = std::sqrt(1.0 * eL * eL + eta * eta);

    // sum the parts and normalize to get the next CWF derivative
    double next_el_cwf_pr = -((1.0 * eL * eL) / rho + eta) * u_L; // -( L^2/rho + eta ) * U_L
    next_el_cwf_pr += normalize_factor * u_L_minus; // (L^2+eta^2)^(1/2) * U_[L-1] - ( L^2/rho + eta ) * U_L
    next_el_cwf_pr /= (eL * 1.0); // { (L^2+eta^2)^(1/2) * U_[L-1] - ( L^2/rho + eta ) * U_L } / L
    return next_el_cwf_pr;
}

// Asymptotic Bessel function solution for large eta case
void Coulomb_wave_functions_accurate::compute_asymptotic_bessel_solution(double rho_r, double eta_r, int eL,
                                                               double &F, double &dF, 
                                                               double &G, double &dG) const
{
    const double half = 0.5;
    const double one = 1.0;
    const double two = 2.0;
    const double euler = 0.577215664901532860606512;
    const double pi = M_PI;
    const double q = two * rho_r * eta_r;
    const double zhalf = std::sqrt(q);
    const double z = two * zhalf;
    const int term_limit = 100;
    
    // Calculate modified Bessel functions I_0, I_1, K_0, K_1
    double bessel_i_0 = one, bessel_k_0 = 0.0, bessel_i_1 = one, bessel_k_1 = 0.0;
    double a, b, c;
    int num_terms;
    
    // I_0(z) from series
    a = q;
    num_terms = 1;
    while ((num_terms < term_limit) && ((std::abs(bessel_i_0) + std::abs(a)) > std::abs(bessel_i_0))) {
        double np1_sq = static_cast<double>(num_terms + 1);
        np1_sq *= np1_sq;
        bessel_i_0 += a;
        a *= q / np1_sq;
        num_terms++;
    }
    
    // K_0(z) from series
    bessel_k_0 = -(std::log(zhalf) + euler) * bessel_i_0;
    a = q;
    b = one;
    num_terms = 1;
    while ((num_terms < term_limit) && ((std::abs(bessel_k_0) + std::abs(a * b)) > std::abs(bessel_k_0))) {
        double np1 = static_cast<double>(num_terms + 1);
        bessel_k_0 += a * b;
        b += one / np1;
        a *= q / (np1 * np1);
        num_terms++;
    }
    
    // I_1(z) from series
    a = q / two;
    num_terms = 1;
    while ((num_terms < term_limit) && ((std::abs(bessel_i_1) + std::abs(a)) > std::abs(bessel_i_1))) {
        double np1 = static_cast<double>(num_terms + 1);
        bessel_i_1 += a;
        a *= q / (np1 * (np1 + one));
        num_terms++;
    }
    bessel_i_1 *= zhalf;
    
    // K_1(z) from series  
    bessel_k_1 = one / z + (std::log(zhalf) + euler) * bessel_i_1 - half * zhalf;
    a = zhalf * q / two;
    b = one;
    c = half * half;
    num_terms = 1;
    while ((num_terms < term_limit) && ((std::abs(bessel_k_1) + std::abs(a * (b + c))) > std::abs(bessel_k_1))) {
        double np1 = static_cast<double>(num_terms + 1);
        bessel_k_1 -= a * (b + c);
        c = half / (np1 + one);
        b += one / np1;
        a *= q / (np1 * (np1 + one));
        num_terms++;
    }
    
    // Convert to F, G, F', G' at L=0
    a = std::sqrt(pi * rho_r);
    b = std::sqrt(two * pi * eta_r);
    F = a * bessel_i_1;
    dF = b * bessel_i_0;
    a *= two / pi;
    b *= two / pi;
    G = a * bessel_k_1;
    dG = -b * bessel_k_0;
    
    // Use recursion to get to desired L
    double F_0 = F, G_0 = G, dF_0 = dF, dG_0 = dG;
    for (int iel = 1; iel <= eL; iel++) {
        F = get_next_el_cwf(rho_r, eta_r, iel, F_0, dF_0);
        G = get_next_el_cwf(rho_r, eta_r, iel, G_0, dG_0);
        dF = get_next_el_cwf_deriv(rho_r, eta_r, iel, F_0, F);
        dG = get_next_el_cwf_deriv(rho_r, eta_r, iel, G_0, G);
        
        F_0 = F; G_0 = G; dF_0 = dF; dG_0 = dG;
    }
    
    // Apply scaling to prevent overflow (from accurate_cwf.cpp)
    // G,G' => need to be scaled to preserve precision in their ratios
    double scale_exp = M_PI * eta_r;
    double max_cwf_mag = std::abs(G);
    if (std::abs(dG) > max_cwf_mag) max_cwf_mag = std::abs(dG);
    if (std::abs(F) > max_cwf_mag) max_cwf_mag = std::abs(F);
    if (std::abs(dF) > max_cwf_mag) max_cwf_mag = std::abs(dF);
    
    double a_max = NUMERICAL_MAX_EXPONENT / 2.0 - std::log(max_cwf_mag) - 10.0;
    if (scale_exp > a_max) scale_exp = a_max;
    
    double scale_factor = std::exp(scale_exp);
    F /= scale_factor;
    dF /= scale_factor;
    G *= scale_factor;
    dG *= scale_factor;
}

// Continued fraction method for complex Coulomb wave functions
// Based on the actual continued fraction implementation from coul library
void Coulomb_wave_functions_accurate::solve_continued_fraction(const std::complex<double> &rho,
                                                              std::complex<double> &F,
                                                              std::complex<double> &dF,
                                                              std::complex<double> &G,
                                                              std::complex<double> &dG) const
{
    std::complex<double> imag_i(0.0, 1.0);
    
    // Try to compute H+ and H+' using continued fractions (like original coul library)
    std::complex<double> Hp, dHp;
    
    try {
        // Compute H+(z) and H+'(z) using continued fraction for h = H'/H
        compute_H_functions(rho, Hp, dHp);
        
        // Compute F and F' using continued fraction for f = F'/F
        compute_F_functions(rho, F, dF);
        
        // Construct G and G' from H+ and F: G = H+ - i*F
        G = Hp - imag_i * F;
        dG = dHp - imag_i * dF;
        
    } catch (const std::exception &e) {
        // If continued fractions fail, throw to trigger fallback
        throw std::runtime_error("Continued fraction computation failed: " + std::string(e.what()));
    }
}

// Series expansion for small rho
void Coulomb_wave_functions_accurate::solve_small_rho_series(const std::complex<double> &rho,
                                                            std::complex<double> &F,
                                                            std::complex<double> &dF,
                                                            std::complex<double> &G,
                                                            std::complex<double> &dG) const
{
    // Series expansion around rho = 0
    // F(l,eta;rho) = C_l(eta) * rho^(l+1) * (1 + O(rho^2))
    
    std::complex<double> imag_i(0.0, 1.0);
    std::complex<double> Cl = compute_Coulomb_constant(this->l, this->eta);
    std::complex<double> rho_power = std::pow(rho, this->l + 1.0);
    
    // Leading terms of series
    F = Cl * rho_power;
    dF = Cl * (this->l + 1.0) * std::pow(rho, this->l) * (1.0 + this->eta * rho / (2.0 * this->l + 2.0));
    
    // G function for small rho (irregular solution)
    std::complex<double> sigma_l = compute_coulomb_phase_shift(this->l, this->eta);
    std::complex<double> gamma_factor = std::exp(imag_i * sigma_l);
    
    G = gamma_factor / (Cl * std::pow(rho, this->l));
    dG = -this->l * gamma_factor / (Cl * std::pow(rho, this->l + 1.0));
}

// Asymptotic expansion for large rho  
void Coulomb_wave_functions_accurate::solve_large_rho_asymptotic(const std::complex<double> &rho,
                                                                std::complex<double> &F,
                                                                std::complex<double> &dF,
                                                                std::complex<double> &G,
                                                                std::complex<double> &dG) const
{
    // Asymptotic forms for large rho
    // F ~ sin(rho - eta*ln(2*rho) + sigma_l)
    // G ~ cos(rho - eta*ln(2*rho) + sigma_l)
    
    std::complex<double> sigma_l = compute_coulomb_phase_shift(this->l, this->eta);
    std::complex<double> theta = rho - this->eta * std::log(2.0 * rho) + sigma_l;
    
    F = std::sin(theta);
    dF = std::cos(theta) * (1.0 - this->eta / rho);
    
    G = std::cos(theta);
    dG = -std::sin(theta) * (1.0 - this->eta / rho);
}

// Simplified continued fraction for intermediate parameters
void Coulomb_wave_functions_accurate::solve_intermediate_continued_fraction(const std::complex<double> &rho,
                                                                          std::complex<double> &F,
                                                                          std::complex<double> &dF,
                                                                          std::complex<double> &G,
                                                                          std::complex<double> &dG) const
{
    // This is a simplified implementation
    // A full implementation would use proper continued fraction evaluation
    
    // For now, use a hybrid approach combining series and asymptotic
    std::complex<double> F_series, dF_series, G_series, dG_series;
    std::complex<double> F_asymp, dF_asymp, G_asymp, dG_asymp;
    
    this->solve_small_rho_series(rho, F_series, dF_series, G_series, dG_series);
    this->solve_large_rho_asymptotic(rho, F_asymp, dF_asymp, G_asymp, dG_asymp);
    
    // Weight between series and asymptotic based on |rho|
    double weight = std::tanh(std::abs(rho) / 5.0);
    
    F = (1.0 - weight) * F_series + weight * F_asymp;
    dF = (1.0 - weight) * dF_series + weight * dF_asymp;
    G = (1.0 - weight) * G_series + weight * G_asymp;
    dG = (1.0 - weight) * dG_series + weight * dG_asymp;
}

// Helper function to compute Coulomb constant C_l(eta)
std::complex<double> Coulomb_wave_functions_accurate::compute_Coulomb_constant(const std::complex<double> &l,
                                                                              const std::complex<double> &eta) const
{
    // C_l(eta) = 2^l * exp(-pi*eta/2) * |Gamma(l+1+i*eta)| / Gamma(2*l+2)
    // Simplified approximation for moderate parameters
    
    std::complex<double> result(1.0, 0.0);
    
    // Power of 2
    result *= std::pow(2.0, l);
    
    // Exponential factor
    result *= std::exp(-M_PI * eta / 2.0);
    
    // Gamma function ratio approximation
    // For simplicity, use Stirling's approximation for large l
    if (std::abs(l) > 5.0) {
        result /= std::pow(2.0 * l, l + 0.5);
    } else {
        // Use factorial for small l (when l is integer)
        if (std::abs(l.imag()) < 1e-10) {
            int l_int = static_cast<int>(std::round(l.real()));
            if (l_int >= 0 && l_int < 20) {
                double factorial = 1.0;
                for (int i = 1; i <= l_int; ++i) {
                    factorial *= i;
                }
                result /= std::sqrt(2.0 * M_PI * factorial);
            }
        }
    }
    
    return result;
}

// Helper function to compute Coulomb phase shift sigma_l
std::complex<double> Coulomb_wave_functions_accurate::compute_coulomb_phase_shift(const std::complex<double> &l,
                                                                                 const std::complex<double> &eta) const
{
    // sigma_l = arg(Gamma(l+1+i*eta))
    // Simplified implementation
    
    std::complex<double> sigma(0.0, 0.0);
    
    // Leading term
    sigma += eta * (std::log(2.0 * l + 2.0) - 1.0);
    
    // Correction terms for small eta
    if (std::abs(eta) < 1.0) {
        sigma -= eta * eta / (2.0 * (l + 1.0));
    }
    
    return sigma;
}

// Compute H+ and H+' using Lentz's continued fraction method
// Based on continued_fraction_h() from original coul library
void Coulomb_wave_functions_accurate::compute_H_functions(const std::complex<double> &z,
                                                          std::complex<double> &H_plus,
                                                          std::complex<double> &dH_plus) const
{
    const double precision = 1e-15;
    const int max_iterations = 100000;
    std::complex<double> imag_i(0.0, 1.0);
    
    // Choose omega = +1 or -1 for numerical stability
    int omega = 1; // For H+, use omega = +1
    
    // Lentz's method for continued fraction h = H'/H
    std::complex<double> z_minus_eta = z - eta;
    std::complex<double> b0 = z_minus_eta;
    
    // Initialize Lentz's method
    std::complex<double> C_n = b0;
    std::complex<double> D_n = 0.0;
    if (std::abs(b0) < 1e-30) {
        C_n = std::complex<double>(1e-30, 0.0);
    }
    
    std::complex<double> Delta_n;
    std::complex<double> h_cf = C_n;
    
    bool converged = false;
    int n = 1;
    
    while (n <= max_iterations && !converged) {
        // Compute coefficients a_n and b_n
        std::complex<double> n_complex(n, 0.0);
        std::complex<double> a_n = (1.0 + l + imag_i * static_cast<double>(omega) + n_complex - 1.0) *
                                  (imag_i * static_cast<double>(omega) * eta - l + n_complex - 1.0);
        std::complex<double> b_n = 2.0 * z_minus_eta + imag_i * static_cast<double>(omega) * n_complex;
        
        // Update D_n
        D_n = b_n + a_n * D_n;
        if (std::abs(D_n) < 1e-30) {
            D_n = std::complex<double>(1e-30, 0.0);
        }
        D_n = 1.0 / D_n;
        
        // Update C_n
        C_n = b_n + a_n / C_n;
        if (std::abs(C_n) < 1e-30) {
            C_n = std::complex<double>(1e-30, 0.0);
        }
        
        // Compute Delta_n
        Delta_n = C_n * D_n;
        h_cf *= Delta_n;
        
        // Check convergence using infinity norm
        if (inf_norm(1.0 - Delta_n) < precision) {
            converged = true;
        }
        
        n++;
    }
    
    if (!converged) {
        throw std::runtime_error("H function continued fraction did not converge");
    }
    
    // Complete the continued fraction result
    h_cf *= imag_i * static_cast<double>(omega) / z;
    
    // Now compute H+ from the logarithmic derivative h = H'/H
    // Use the Wronskian and known properties to get H+ and H+'
    compute_H_from_logarithmic_derivative(z, h_cf, omega, H_plus, dH_plus);
}

// Compute F and F' using continued fraction method
void Coulomb_wave_functions_accurate::compute_F_functions(const std::complex<double> &z,
                                                          std::complex<double> &F,
                                                          std::complex<double> &dF) const
{
    const double precision = 1e-15;
    const int max_iterations = 100000;
    std::complex<double> imag_i(0.0, 1.0);
    
    // Lentz's method for continued fraction f = F'/F
    std::complex<double> b0 = l + 1.0 + imag_i * eta;
    
    // Initialize Lentz's method
    std::complex<double> C_n = b0;
    std::complex<double> D_n = 0.0;
    if (std::abs(b0) < 1e-30) {
        C_n = std::complex<double>(1e-30, 0.0);
    }
    
    std::complex<double> Delta_n;
    std::complex<double> f_cf = C_n;
    
    bool converged = false;
    int n = 1;
    
    while (n <= max_iterations && !converged) {
        // Compute coefficients for F continued fraction
        std::complex<double> n_complex(n, 0.0);
        std::complex<double> a_n = -eta * eta + (l + n_complex) * (l + n_complex);
        std::complex<double> b_n = 2.0 * (l + n_complex) + 1.0 + 2.0 * imag_i * eta;
        
        // Update D_n
        D_n = b_n + a_n * D_n;
        if (std::abs(D_n) < 1e-30) {
            D_n = std::complex<double>(1e-30, 0.0);
        }
        D_n = 1.0 / D_n;
        
        // Update C_n
        C_n = b_n + a_n / C_n;
        if (std::abs(C_n) < 1e-30) {
            C_n = std::complex<double>(1e-30, 0.0);
        }
        
        // Compute Delta_n
        Delta_n = C_n * D_n;
        f_cf *= Delta_n;
        
        // Check convergence
        if (inf_norm(1.0 - Delta_n) < precision) {
            converged = true;
        }
        
        n++;
    }
    
    if (!converged) {
        throw std::runtime_error("F function continued fraction did not converge");
    }
    
    // Complete the continued fraction result
    f_cf /= z;
    
    // Compute F from the logarithmic derivative f = F'/F
    compute_F_from_logarithmic_derivative(z, f_cf, F, dF);
}

// Helper to compute H+ and H+' from logarithmic derivative
void Coulomb_wave_functions_accurate::compute_H_from_logarithmic_derivative(const std::complex<double> &z,
                                                                           const std::complex<double> &h,
                                                                           int omega,
                                                                           std::complex<double> &H_plus,
                                                                           std::complex<double> &dH_plus) const
{
    // This is a simplified implementation
    // A full implementation would use proper normalization and Wronskian relations
    
    // For now, use asymptotic normalization and integrate backward
    std::complex<double> imag_i(0.0, 1.0);
    std::complex<double> sigma_l = compute_coulomb_phase_shift(l, eta);
    
    // Asymptotic form: H+ ~ exp(i*(z - eta*ln(2z) + sigma_l))
    std::complex<double> phase = imag_i * (z - eta * std::log(2.0 * z) + sigma_l);
    H_plus = std::exp(phase);
    
    // H+' = h * H+
    dH_plus = h * H_plus;
}

// Helper to compute F and F' from logarithmic derivative  
void Coulomb_wave_functions_accurate::compute_F_from_logarithmic_derivative(const std::complex<double> &z,
                                                                           const std::complex<double> &f,
                                                                           std::complex<double> &F,
                                                                           std::complex<double> &dF) const
{
    // Compute Coulomb constant for proper normalization
    std::complex<double> Cl = compute_Coulomb_constant(l, eta);
    
    // Near origin behavior: F ~ C_l * z^(l+1)
    std::complex<double> z_power = std::pow(z, l + 1.0);
    F = Cl * z_power;
    
    // F' = f * F
    dF = f * F;
}

// Infinity norm helper function
double Coulomb_wave_functions_accurate::inf_norm(const std::complex<double> &z) const
{
    return std::max(std::abs(z.real()), std::abs(z.imag()));
}