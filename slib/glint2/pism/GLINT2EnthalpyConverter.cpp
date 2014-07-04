#include <petsc.h>  // for PetscErrorPrintf, etc.
#include "pism_const.hh"
#include <glint2/pism/GLINT2EnthalpyConverter.hpp>
#include "ConfigI.hh"

namespace glint2 {
namespace gpism {

GLINT2EnthalpyConverter::GLINT2EnthalpyConverter(const pism::ConfigI &config) :
	EnthalpyConverter(config)
{
  c_w   = config.get("water_specific_heat_capacity");            // J kg-1 K-1
}

//! Simple view of state of EnthalpyConverter.  viewer==NULL sends to stdout.
PetscErrorCode GLINT2EnthalpyConverter::viewConstants(PetscViewer viewer) const {
  PetscErrorCode ierr;

  ierr = pism::EnthalpyConverter::viewConstants(viewer); CHKERRQ(ierr);

  ierr = PetscViewerASCIIPrintf(viewer,
    "\n<showing GLINT2EnthalpyConverter constants:\n"); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
      "   c_w   = %12.5f (J kg-1 K-1)\n",c_w); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
      ">\n"); CHKERRQ(ierr);

  return 0;

}

//! Get enthalpies E_s(p) and E_l(p) (endpoints of temperate ice enthalpy range) from pressure p.
/*! Ice at enthalpy \f$E\f$ is temperate if \f$E_s(p) < E < E_l(p)\f$:
     \f[ E_s(p) = c_i (T_m(p) - T_0), \f]
     \f[ E_l(p) = E_s(p) + L. \f]
 */
PetscErrorCode GLINT2EnthalpyConverter::getEnthalpyInterval(
                       double p, double &E_s, double &E_l) const {
  E_s = getEnthalpyCTS(p);
  E_l = E_s + L_from_p(p);
  return 0;
}

//! Get liquid water fraction from enthalpy and pressure.
/*!
From \ref AschwandenBuelerKhroulevBlatter,
   \f[ \omega(E,p) = \begin{cases}  0.0,            & E \le E_s(p), \\
                                    (E-E_s(p)) / L, & E_s(p) < E < E_l(p).
                     \end{cases} \f]

We do not allow liquid water (i.e. water fraction \f$\omega=1.0\f$) so we return
error code 1 if \f$E \ge E_l(p)\f$, but we still compute \f$\omega=1.0\f$.
 */
PetscErrorCode GLINT2EnthalpyConverter::getWaterFraction(double E, double p, double &omega) const {
  double E_s, E_l;
  PetscErrorCode ierr = getEnthalpyInterval(p, E_s, E_l);
#if PISM_DEBUG==1
  CHKERRQ(ierr);
#endif
  if (E >= E_l) {
    omega = 1.0;
    return 1;
  }
  if (E <= E_s) {
    omega = 0.0;
  } else {
    omega = (E - E_s) / L_from_p(p);
  }
  return 0;
}

//! Compute enthalpy from absolute temperature, liquid water fraction, and pressure.
/*! This is an inverse function to the functions \f$T(E,p)\f$ and
\f$\omega(E,p)\f$ [\ref AschwandenBuelerKhroulevBlatter].  It returns:
  \f[E(T,\omega,p) =
       \begin{cases}
         c_i (T - T_0),     & T < T_m(p) \quad\text{and}\quad \omega = 0, \\
         E_s(p) + \omega L, & T = T_m(p) \quad\text{and}\quad \omega \ge 0.
       \end{cases} \f]
Certain cases are not allowed and return errors:
- \f$T<=0\f$ (error code 1)
- \f$\omega < 0\f$ or \f$\omega > 1\f$ (error code 2)
- \f$T>T_m(p)\f$ (error code 3)
- \f$T<T_m(p)\f$ and \f$\omega > 0\f$ (error code 4)
These inequalities may be violated in the sixth digit or so, however.
 */
PetscErrorCode GLINT2EnthalpyConverter::getEnth(
                  double T, double omega, double p, double &E) const {
  const double T_m = getMeltingTemp(p);
#if (PISM_DEBUG==1)
  if (T <= 0.0) {
    SETERRQ1(PETSC_COMM_SELF, 1,"\n\nT = %f <= 0 is not a valid absolute temperature\n\n",T);
  }
  if ((omega < 0.0 - 1.0e-6) || (1.0 + 1.0e-6 < omega)) {
    SETERRQ1(PETSC_COMM_SELF, 2,"\n\nwater fraction omega=%f not in range [0,1]\n\n",omega);
  }
  if (T > T_m + 1.0e-6) {
    SETERRQ2(PETSC_COMM_SELF, 3,"T=%f exceeds T_m=%f; not allowed\n\n",T,T_m);
  }
  if ((T < T_m - 1.0e-6) && (omega > 0.0 + 1.0e-6)) {
    SETERRQ3(PETSC_COMM_SELF, 4,"T < T_m AND omega > 0 is contradictory\n\n",T,T_m,omega);
  }
#endif
  if (T < T_m) {
    E = c_i * (T - T_0);
  } else {
    E = getEnthalpyCTS(p) + omega * L_from_Tm(T_m);
  }
  return 0;
}

//! Returns enthalpy for temperate ice with a given liquid fraction.
/*! Computes
  \f[E = E_s(p) + \omega L.\f]
Only the following case returns an error:
- \f$\omega < 0\f$ or \f$\omega > 1\f$ (error code 2)
These inequalities may be violated in the sixth digit or so, however.
 */
PetscErrorCode GLINT2EnthalpyConverter::getEnthAtWaterFraction(
                        double omega, double p, double &E) const {
#if (PISM_DEBUG==1)
  if ((omega < 0.0 - 1.0e-6) || (1.0 + 1.0e-6 < omega)) {
    SETERRQ1(PETSC_COMM_SELF, 2,"\n\nwater fraction omega=%f not in range [0,1]\n\n",omega);
  }
#endif
  E = getEnthalpyCTS(p) + omega * L_from_p(p);
  return 0;
}

}}	// namespace glint2::gpism
