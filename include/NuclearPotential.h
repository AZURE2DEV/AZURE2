#ifndef NUCLEAR_POTENTIAL_H
#define NUCLEAR_POTENTIAL_H

#include <functional>
#include <string>
#include <stdexcept>

/**
 * @brief Base class for nuclear potentials
 */
class NuclearPotential {
public:
    virtual ~NuclearPotential() = default;

    /**
     * @brief Evaluate the potential at radius r (in fm)
     * @param r Radius in fm
     * @return Potential value in MeV
     */
    virtual double operator()(double r) const = 0;

    /**
     * @brief Get potential type name
     * @return String describing the potential type
     */
    virtual std::string name() const = 0;
};

/**
 * @brief Woods-Saxon potential: V(r) = -V0 / (1 + exp((r - R) / a))
 */
class WoodsSaxonPotential : public NuclearPotential {
public:
    /**
     * @brief Constructor
     * @param V0 Depth of the potential well (MeV)
     * @param R Nuclear radius (fm)
     * @param a Surface diffuseness (fm)
     */
    WoodsSaxonPotential(double V0, double R, double a)
        : V0_(V0), R_(R), a_(a) {
        if (V0 <= 0.0 || R <= 0.0 || a <= 0.0) {
            throw std::invalid_argument("WoodsSaxon: V0, R, a must be positive");
        }
    }

    double operator()(double r) const override {
        return -V0_ / (1.0 + std::exp((r - R_) / a_));
    }

    std::string name() const override { return "Woods-Saxon"; }

    double get_V0() const { return V0_; }
    double get_R() const { return R_; }
    double get_a() const { return a_; }

private:
    double V0_;  // Depth (MeV)
    double R_;   // Radius (fm)
    double a_;   // Diffuseness (fm)
};

/**
 * @brief Gaussian potential: V(r) = -V0 * exp(-(r / r0)^2)
 */
class GaussianPotential : public NuclearPotential {
public:
    /**
     * @brief Constructor
     * @param V0 Depth of the potential well (MeV)
     * @param r0 Gaussian width parameter (fm)
     */
    GaussianPotential(double V0, double r0)
        : V0_(V0), r0_(r0) {
        if (V0 <= 0.0 || r0 <= 0.0) {
            throw std::invalid_argument("Gaussian: V0, r0 must be positive");
        }
    }

    double operator()(double r) const override {
        double x = r / r0_;
        return -V0_ * std::exp(-x * x);
    }

    std::string name() const override { return "Gaussian"; }

    double get_V0() const { return V0_; }
    double get_r0() const { return r0_; }

private:
    double V0_;   // Depth (MeV)
    double r0_;   // Width (fm)
};

/**
 * @brief Centrifugal potential: V_l(r) = hbar^2/(2*mu) * l(l+1)/r^2
 */
class CentrifugalPotential : public NuclearPotential {
public:
    /**
     * @brief Constructor
     * @param l Orbital angular momentum quantum number
     * @param redmass_MeV Reduced mass in MeV/c^2
     * @param hbarc hbar*c in MeV*fm (default 197.3269804)
     */
    CentrifugalPotential(int l, double redmass_MeV, double hbarc = 197.3269804)
        : l_(l), redmass_MeV_(redmass_MeV), hbarc_(hbarc) {
        if (l < 0) {
            throw std::invalid_argument("Angular momentum l must be non-negative");
        }
        if (redmass_MeV <= 0.0) {
            throw std::invalid_argument("Reduced mass must be positive");
        }
        if (hbarc <= 0.0) {
            throw std::invalid_argument("hbarc must be positive");
        }
    }

    double operator()(double r) const override {
        if (r <= 0.0) return 0.0;
        double hbar2_over_2mu = (hbarc_ * hbarc_) / (2.0 * redmass_MeV_);
        return hbar2_over_2mu * l_ * (l_ + 1) / (r * r);
    }

    std::string name() const override {
        return "Centrifugal(l=" + std::to_string(l_) + ")";
    }

    int get_l() const { return l_; }
    double get_redmass_MeV() const { return redmass_MeV_; }
    double get_hbarc() const { return hbarc_; }

private:
    int l_;
    double redmass_MeV_;
    double hbarc_;
};

/**
 * @brief Combined potential: Coulomb + Nuclear + Centrifugal
 * Coulomb: uniform sphere for r < Rc, point charge for r >= Rc
 * Nuclear: user-supplied potential
 * Centrifugal: l(l+1)/(2*mu*r^2)
 */
class CoulombNuclearPotential : public NuclearPotential {
public:
    /**
     * @brief Constructor
     * @param Z1 Projectile charge
     * @param Z2 Target charge
     * @param nuclear_potential Nuclear potential (Woods-Saxon, Gaussian, etc.)
     * @param l Orbital angular momentum
     * @param redmass_MeV Reduced mass in MeV/c^2
     * @param Rc Coulomb radius - switch from uniform sphere to point charge (fm)
     * @param hbarc hbar*c in MeV*fm (default 197.3269804)
     */
    CoulombNuclearPotential(int Z1, int Z2,
                            const NuclearPotential& nuclear_potential,
                            int l = 0,
                            double redmass_MeV = 931.494,
                            double Rc = 3.5,
                            double hbarc = 197.3269804)
        : Z1_(Z1), Z2_(Z2), nuclear_(nuclear_potential),
          l_(l), redmass_MeV_(redmass_MeV), Rc_(Rc), hbarc_(hbarc) {
        if (Rc <= 0.0) {
            throw std::invalid_argument("Coulomb radius Rc must be positive");
        }
        if (l < 0) {
            throw std::invalid_argument("Angular momentum l must be non-negative");
        }
        if (redmass_MeV <= 0.0) {
            throw std::invalid_argument("Reduced mass must be positive");
        }
    }

    double operator()(double r) const override {
        // Coulomb potential: uniform sphere (r < Rc) or point charge (r >= Rc)
        constexpr double e2 = 1.4399764;  // MeV*fm
        double Vc;

        if (r < Rc_) {
            // Uniform charged sphere: Vc = Z1*Z2*e^2/(2*Rc) * (3 - (r/Rc)^2)
            double x = r / Rc_;
            Vc = Z1_ * Z2_ * e2 / (2.0 * Rc_) * (3.0 - x * x);
        } else {
            // Point charge: Vc = Z1*Z2*e^2/r
            Vc = Z1_ * Z2_ * e2 / r;
        }

        // Add nuclear potential
        double V_nuclear = nuclear_(r);

        // Add centrifugal potential
        double V_centrifugal = 0.0;
        if (r > 0.0) {
            double hbar2_over_2mu = (hbarc_ * hbarc_) / (2.0 * redmass_MeV_);
            V_centrifugal = hbar2_over_2mu * l_ * (l_ + 1) / (r * r);
        }

        return Vc + V_nuclear + V_centrifugal;
    }

    std::string name() const override {
        return "Coulomb + " + nuclear_.name() + " + Centrifugal(l=" + std::to_string(l_) + ")";
    }

    double get_Rc() const { return Rc_; }
    int get_l() const { return l_; }
    double get_redmass_MeV() const { return redmass_MeV_; }

private:
    int Z1_, Z2_;
    const NuclearPotential& nuclear_;
    int l_;
    double redmass_MeV_;
    double Rc_;
    double hbarc_;
};

#endif // NUCLEAR_POTENTIAL_H
