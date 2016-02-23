#include "enthalpyConverter.hh"

namespace glint2 {
namespace gpism {

class GLINT2EnthalpyConverter : public pism::EnthalpyConverter {
public:

  double L_from_Tm(double T_m) const
    { return L + (c_w - c_i) * (T_m - 273.15); }

  double L_from_p(double p) const
    { return L_from_Tm(getMeltingTemp(p)); }

protected:
  double c_w;

public:

  GLINT2EnthalpyConverter(const pism::ConfigI &config);

  PetscErrorCode viewConstants(PetscViewer viewer) const;

  PetscErrorCode getEnthalpyInterval(
                         double p, double &E_s, double &E_l) const;

  PetscErrorCode getWaterFraction(double E, double p, double &omega) const;

  PetscErrorCode getEnth(
                    double T, double omega, double p, double &E) const;

  PetscErrorCode getEnthAtWaterFraction(
                          double omega, double p, double &E) const;

};


}}	// namespace glint2::gpism
