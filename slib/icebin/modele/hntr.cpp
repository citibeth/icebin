#include <cmath>
#include <icebin/error.hpp>
#include <icebin/modele/hntr.hpp>

using namespace blitz;
using namespace spsparse;
using namespace ibmisc;

namespace icebin {
namespace modele {

HntrGrid::HntrGrid(HntrGrid const &other)
    : spec(other.spec),
    dxyp(blitz::Range(other.dxyp.lbound(0), other.dxyp.ubound(0)))
{
    dxyp = other.dxyp;
}


HntrGrid::HntrGrid(HntrSpec const &_spec) :
    spec(_spec)
{
    init();
}


void HntrGrid::init()
{
    auto const jm(spec.jm);

    dxyp.reference(blitz::Array<double,1>(blitz::Range(1,jm)));

    // Calculate the sperical area of grid cells
    double dLON = (2.*M_PI) / spec.im;
    double dLAT = M_PI / jm;
    for (int j=1; j<=spec.jm; ++j) {
        double SINS = sin(dLAT*(j-jm/2-1));
        double SINN = sin(dLAT*(j-jm/2));
        dxyp(j) = dLON * (SINN - SINS);
    }
}

void HntrGrid::ncio(ibmisc::NcIO &ncio, std::string const &vname)
{
    spec.ncio(ncio, vname);
    init();
}

// =============================================================

Hntr::Hntr(double yp17, HntrSpec const &_Bgrid, HntrSpec const &_Agrid, double _DATMIS)
    : Agrid(HntrGrid(_Agrid)), Bgrid(HntrGrid(_Bgrid)),
      SINA(Range(0, Agrid.spec.jm)),
      SINB(Range(0, Bgrid.spec.jm)),
      FMIN(Range(1, Bgrid.spec.im)),
      FMAX(Range(1, Bgrid.spec.im)),
      IMIN(Range(1, Bgrid.spec.im)),
      IMAX(Range(1, Bgrid.spec.im)),
      GMIN(Range(1, Bgrid.spec.jm)),
      GMAX(Range(1, Bgrid.spec.jm)),
      JMIN(Range(1, Bgrid.spec.jm)),
      JMAX(Range(1, Bgrid.spec.jm)),
      DATMIS(_DATMIS)
{

    partition_east_west();
    partition_north_south();
}

// Partitions in east-west (I) direction
// Domain, around the globe, is scaled to fit from 0 to IMA*Bgrid.spec.im
void Hntr::partition_east_west()
{
    double DIA = Bgrid.spec.im;  //  width of single A grid cell in scaled domain
    int IA = 1;
    double RIA = (IA+Agrid.spec.offi - Agrid.spec.im)*Bgrid.spec.im;  //  scaled longitude of eastern edge
    int IB  = Bgrid.spec.im;
    for (int IBp1=1; IBp1 <= Bgrid.spec.im; ++IBp1) {
        double RIB = (IBp1-1+Bgrid.spec.offi)*Agrid.spec.im;    //  scaled longitude of eastern edge
        while (RIA < RIB) {
            IA  += 1;
            RIA += DIA;
        }

        if (RIA == RIB) {
        // Eastern edges of cells IA of grid A and IB of grid B coincide
            IMAX(IB) = IA;
            FMAX(IB) = 0;
            IA += 1;
            RIA += DIA;
            IMIN(IBp1) = IA;
            FMIN(IBp1) = 0;
        } else {
        // Cell IA of grid A contains western edge of cell IB of grid B
            IMAX(IB) = IA;
            FMAX(IB) = (RIA-RIB)/DIA;
            IMIN(IBp1) = IA;
            FMIN(IBp1) = 1-FMAX(IB);
        }
        IB = IBp1;
    }
    IMAX(Bgrid.spec.im) += Agrid.spec.im;
}


void Hntr::partition_north_south()
{
    // Convert minutes to radians
    const double MIN_TO_RAD = (2. * M_PI) / (360*60);

    // ------------------------------------------------
    // Partitions in the north-south (J) direction
    // Domain is measured in minutes (1/60-th of a degree)
    double FJEQA = .5*(1+Agrid.spec.jm);
    for (int JA=1; JA <= Agrid.spec.jm-1; ++JA) {
        double RJA = (JA + .5-FJEQA) * Agrid.spec.dlat;  //  latitude in minutes of northern edge
        SINA(JA) = sin(RJA * MIN_TO_RAD);
    }
    SINA(0) = -1;
    SINA(Agrid.spec.jm)=  1;

    // -----------
    double FJEQB = .5*(1+Bgrid.spec.jm);
    for (int JB=1; JB <= Bgrid.spec.jm-1; ++JB) {
        double RJB = (JB+.5-FJEQB)*Bgrid.spec.dlat;  //  latitude in minutes of northern edge
        SINB(JB) = sin(RJB * MIN_TO_RAD);
    }
    SINB(0) = -1;
    SINB(Bgrid.spec.jm) = 1;

    // -----------
    JMIN(1) = 1;
    GMIN(1) = 0;
    int JA = 1;
    for (int JB=1; JB <= Bgrid.spec.jm-1; ++JB) {
        while (SINA(JA) < SINB(JB)) ++JA;

        if (SINA(JA) == SINB(JB)) {
            //  Northern edges of cells JA of grid A and JB of grid B coincide
            JMAX(JB) = JA;
            GMAX(JB) = 0;
            JA += 1;
            JMIN(JB+1) = JA;
            GMIN(JB+1) = 0;
        } else {
            // Cell JA of grid A contains northern edge of cell JB of grid B
            JMAX(JB) = JA;
            GMAX(JB) = SINA(JA) - SINB(JB);
            JMIN(JB+1) = JA;
            GMIN(JB+1) = SINB(JB) - SINA(JA-1);
        }
    }
    JMAX(Bgrid.spec.jm) = Agrid.spec.jm;
    GMAX(Bgrid.spec.jm) = 0;

}

void Hntr::regrid1(
    blitz::Array<double,1> const &WTA,
    blitz::Array<double,1> const &A,
    blitz::Array<double,1> &B,
    bool mean_polar) const
{
    // Check array dimensions
    if ((WTA.extent(0) != Agrid.spec.size()) ||
        (A.extent(0) != Agrid.spec.size()) ||
        (B.extent(0) != Bgrid.spec.size()))
    {
        (*icebin_error)(-1, "Error in dimensions: (%d, %d, %d) vs. (%d, %d)\n",
            WTA.extent(0), A.extent(0), B.extent(0),
            Agrid.spec.size(), Bgrid.spec.size());
    }

    // ------------------
    // Interpolate the A grid onto the B grid

    for (int JB=1; JB <= Bgrid.spec.jm; ++JB) {
        int JAMIN = JMIN(JB);
        int JAMAX = JMAX(JB);

        for (int IB=1; IB <= Bgrid.spec.im; ++IB) {
            int const IJB = IB + Bgrid.spec.im * (JB-1);
            double WEIGHT= 0;
            double VALUE = 0;
            int const IAMIN = IMIN(IB);
            int const IAMAX = IMAX(IB);
            for (int JA=JAMIN; JA <= JAMAX; ++JA) {
                double G = SINA(JA) - SINA(JA-1);
                if (JA==JAMIN) G -= GMIN(JB);
                if (JA==JAMAX) G -= GMAX(JB);

                for (int IAREV=IAMIN; IAREV <= IAMAX; ++IAREV) {
                    int const IA  = 1 + ((IAREV-1) % Agrid.spec.im);
                    int const IJA = IA + Agrid.spec.im * (JA-1);
                    double F = 1;
                    if (IAREV==IAMIN) F -= FMIN(IB);
                    if (IAREV==IAMAX) F -= FMAX(IB);

                    double const wt = F*G*WTA(IJA);
                    WEIGHT += wt;
                    VALUE  += wt*A(IJA);
                }
            }
            B(IJB) = (WEIGHT == 0 ? DATMIS : VALUE / WEIGHT);
        }
    }

    if (mean_polar) {
        // Replace individual values near the poles by longitudinal mean
        for (int JB=1; JB <= Bgrid.spec.jm; JB += Bgrid.spec.jm-1) {
            double BMEAN  = DATMIS;
            double WEIGHT = 0;
            double VALUE  = 0;
            for (int IB=1; ; ++IB) {
                if (IB > Bgrid.spec.im) {
                    if (WEIGHT != 0) BMEAN = VALUE / WEIGHT;
                    break;
                }
                int IJB = IB + Bgrid.spec.im * (JB-1);
                if (B(IJB) == DATMIS) break;
                WEIGHT += 1;
                VALUE  += B(IJB);
            }
            for (int IB=1; IB <= Bgrid.spec.im; ++IB) {
                int IJB = IB + Bgrid.spec.im * (JB-1);
                B(IJB) = BMEAN;
            }
        }
    }
}


#if 0
/** Produce an almost-identity matrix that takes the mean of the polar grid cells */
void Hntr::mean_polar_matrix()
{
    if (mean_polar) {
        // Replace individual values near the poles by longitudinal mean
        for (int JB=1; JB <= Bgrid.spec.jm; JB += Bgrid.spec.jm-1) {
            double BMEAN  = DATMIS;
            double WEIGHT = 0;
            double VALUE  = 0;
            for (int IB=1; ; ++IB) {
                if (IB > Bgrid.spec.im) {
                    if (WEIGHT != 0) BMEAN = VALUE / WEIGHT;
                    break;
                }
                int IJB = IB + Bgrid.spec.im * (JB-1);
                if (B(IJB) == DATMIS) break;
                WEIGHT += 1;
                VALUE  += B(IJB);
            }
            for (int IB=1; IB <= Bgrid.spec.im; ++IB) {
                int IJB = IB + Bgrid.spec.im * (JB-1);
                B(IJB) = BMEAN;
            }
        }
    }
}
#endif




// -------------------------------------------------------------


}}



// =================================================================
// Original Fortran Code follows
// (author: Gary Russell)
#if 0
!**** HNTR4.F90   Horizontal Interpolation Program Real*4   2014/12/24
!****
      Subroutine HNTR40 (IMA,JMA,OFFIA,DLATA, IMB,JMB,OFFIB,DLATB, DATMIS)
!****
!**** HNTR40 fills in the common block HNTRCB with coordinate
!**** parameters that will be used by subsequent calls to HNTR4.
!**** The 5 Real input values are expected to be Real*4.
!****
!**** Input: IMA = number of cells in east-west direction of grid A
!****        JMA = number of cells in north-south direction of grid A
!****      OFFIA = number of cells of grid A in east-west direction
!****              from IDL (180) to western edge of cell IA=1
!****      DLATA = minutes of latitude for non-polar cells on grid A
!****        IMB = number of cells in east-west direction of grid B
!****        JMB = number of cells in north-south direction of grid B
!****      OFFIB = number of cells of grid B in east-west direction
!****              from IDL (180) to western edge of cell IB=1
!****      DLATB = minutes of latitude for non-polar cells on grid B
!****     DATMIS = missing data value inserted in output array B when
!****              cell (IB,JB) has integrated value 0 of WTA
!****
!**** Output: common block /HNTRCB/
!**** SINA(JA) = sine of latitude of northern edge of cell JA on grid A
!**** SINB(JB) = sine of latitude of northern edge of cell JB on grid B
!**** FMIN(IB) = fraction of cell IMIN(IB) on grid A west of cell IB
!**** FMAX(IB) = fraction of cell IMAX(IB) on grid A east of cell IB
!**** GMIN(JB) = fraction of cell JMIN(JB) on grid A south of cell JB
!**** GMAX(JB) = fraction of cell JMAX(JB) on grid A north of cell JB
!**** IMIN(IB) = western most cell of grid A that intersects cell IB
!**** IMAX(IB) = eastern most cell of grid A that intersects cell IB
!**** JMIN(JB) = southern most cell of grid A that intersects cell JB
!**** JMAX(JB) = northern most cell of grid A that intersects cell JB
!****
      Implicit None
      Real*8,Parameter :: TWOPI=6.283185307179586477d0
      Integer :: IMA,JMA,IMB,JMB
      Real*4  :: OFFIA,DLATA, OFFIB,DLATB, DATMIS
!**** Local variables
      Integer :: INA,JNA,INB,JNB, IMIN,IMAX,JMIN,JMAX, &
                 IA,IB,JA,JB, IBp1
      Real*4  :: DATMCB
      Real*8  :: SINA,SINB, FMIN,FMAX,GMIN,GMAX, &
                 DIA,DIB, RIA,RIB,RJA,RJB, FJEQA,FJEQB
      Common /HNTRCB/ SINA(0:5401),SINB(0:5401), &
                      FMIN(10800),FMAX(10800),GMIN(5401),GMAX(5401), &
                      IMIN(10800),IMAX(10800),JMIN(5401),JMAX(5401), &
                      DATMCB, INA,JNA,INB,JNB
!****
      INA = IMA  ;  JNA = JMA  ;  INB = IMB  ;  JNB = JMB
      DATMCB = DATMIS
      If (IMA<1 .or. IMA>10800 .or. JMA<1 .or. JMA>5401 .or. IMB<1 .or. IMB>10800 .or. JMB<1 .or. JMB>5401)  GoTo 400
!****
!**** Partitions in east-west (I) direction
!**** Domain, around the globe, is scaled to fit from 0 to IMA*IMB
!****
      DIA = IMB  !  width of single A grid cell in scaled domain
      DIB = IMA  !  width of single B grid cell in scaled domain
      IA  = 1
      RIA = (IA+OFFIA - IMA)*IMB  !  scaled longitude of eastern edge
      IB  = IMB
      Do 150 IBp1=1,IMB
      RIB = (IBp1-1+OFFIB)*IMA    !  scaled longitude of eastern edge
  110 If (RIA-RIB)  120,130,140
  120 IA  = IA  + 1
      RIA = RIA + DIA
      GoTo 110
!**** Eastern edges of cells IA of grid A and IB of grid B coincide
  130 IMAX(IB) = IA
      FMAX(IB) = 0
      IA  = IA  + 1
      RIA = RIA + DIA
      IMIN(IBp1) = IA
      FMIN(IBp1) = 0
      GoTo 150
!**** Cell IA of grid A contains western edge of cell IB of grid B
  140 IMAX(IB) = IA
      FMAX(IB) = (RIA-RIB)/DIA
      IMIN(IBp1) = IA
      FMIN(IBp1) = 1-FMAX(IB)
  150 IB = IBp1
      IMAX(IMB) = IMAX(IMB) + IMA
!        WRITE (0,915) 'IMIN=',IMIN(1:IMB)
!        WRITE (0,915) 'IMAX=',IMAX(1:IMB)
!        WRITE (0,916) 'FMIN=',FMIN(1:IMB)
!        WRITE (0,916) 'FMAX=',FMAX(1:IMB)
!****
!**** Partitions in the north-south (J) direction
!**** Domain is measured in minutes (1/60-th of a degree)
!****
      FJEQA = .5*(1+JMA)
      Do 210 JA=1,JMA-1
      RJA = (JA+.5-FJEQA)*DLATA  !  latitude in minutes of northern edge
  210 SINA(JA) = Sin (RJA*TWOPI/(360*60))
      SINA(0)  = -1
      SINA(JMA)=  1
!****
      FJEQB = .5*(1+JMB)
      Do 220 JB=1,JMB-1
      RJB = (JB+.5-FJEQB)*DLATB  !  latitude in minutes of northern edge
  220 SINB(JB) = Sin (RJB*TWOPI/(360*60))
      SINB(0)  = -1
      SINB(JMB)=  1
!****
      JMIN(1) = 1
      GMIN(1) = 0
      JA = 1
      Do 350 JB=1,JMB-1
  310 If (SINA(JA)-SINB(JB))  320,330,340
  320 JA = JA + 1
      GoTo 310
!**** Northern edges of cells JA of grid A and JB of grid B coincide
  330 JMAX(JB) = JA
      GMAX(JB) = 0
      JA = JA + 1
      JMIN(JB+1) = JA
      GMIN(JB+1) = 0
      GoTo 350
!**** Cell JA of grid A contains northern edge of cell JB of grid B
  340 JMAX(JB) = JA
      GMAX(JB) = SINA(JA) - SINB(JB)
      JMIN(JB+1) = JA
      GMIN(JB+1) = SINB(JB) - SINA(JA-1)
  350 Continue
      JMAX(JMB) = JMA
      GMAX(JMB) = 0
!        WRITE (0,915) 'JMIN=',JMIN(1:JMB)
!        WRITE (0,915) 'JMAX=',JMAX(1:JMB)
!        WRITE (0,916) 'GMIN=',GMIN(1:JMB)
!        WRITE (0,916) 'GMAX=',GMAX(1:JMB)
      Return
!****
!**** Invalid parameters or dimensions out of range
!****
  400 Write (0,940) IMA,JMA,OFFIA,DLATA, IMB,JMB,OFFIB,DLATB, DATMIS
      Stop 400
!****
! 915 Format (/ 1X,A5 / (20I6))
! 916 Format (/ 1X,A5 / (20F6.2))
  940 Format (/ 'Arguments received by HNTRP0 in order:'/ &
              2I12,' = IMA,JMA = array dimensions for A grid'/ &
             E24.8,' = OFFIA   = fractional number of grid cells from IDL to western edge of grid cell I=1'/ &
             E24.8,' = DLATA   = minutes of latitude for interior grid cell'/ &
              2I12,' = IMB,JMB = array dimensions for B grid'/ &
             E24.8,' = OFFIB   = fractional number of grid cells from IDL to western edge of grid cell I=1'/ &
             E24.8,' = DLATB   = minute of latitude for interior grid cell'/ &
             E24.8,' = DATMIS  = missing data value to be put in B array when integrated WTA = 0'// &
             'These arguments are invalid or out of range.')
      EndSubroutine HNTR40


      Subroutine HNTR4 (WTA,A,B)
!****
!**** HNTR4 performs a horizontal interpolation of per unit area or per
!**** unit mass quantities defined on grid A, calculating the quantity
!**** on grid B.  B grid values that cannot be calculated because the
!**** covering A grid boxes have WTA = 0, are set to the value DATMIS.
!**** The area weighted integral of the quantity is conserved.
!**** The 3 Real input values are expected to be Real*4.
!****
!**** Input: WTA = weighting array for values on the A grid
!****          A = per unit area or per unit mass quantity
!**** Output:  B = horizontally interpolated quantity on B grid
!****
      Implicit None
      Real*4  :: WTA(*), A(*), B(*)
      Integer :: IMA,JMA,IMB,JMB, IMIN,IMAX,JMIN,JMAX
      Real*4  :: DATMIS
      Real*8  :: SINA,SINB, FMIN,FMAX,GMIN,GMAX
      Common /HNTRCB/ SINA(0:5401),SINB(0:5401), &
                      FMIN(10800),FMAX(10800),GMIN(5401),GMAX(5401), &
                      IMIN(10800),IMAX(10800),JMIN(5401),JMAX(5401), &
                      DATMIS, IMA,JMA,IMB,JMB
!**** Local variables
      Integer :: IA,JA,IJA, IB,JB,IJB, IAMIN,IAMAX,JAMIN,JAMAX, IAREV
      Real*8  :: WEIGHT,VALUE,F,G
!****
!**** Interpolate the A grid onto the B grid
!****
      Do 20 JB=1,JMB
      JAMIN = JMIN(JB)
      JAMAX = JMAX(JB)
      Do 20 IB=1,IMB
      IJB  = IB + IMB*(JB-1)
      WEIGHT= 0
      VALUE = 0
      IAMIN = IMIN(IB)
      IAMAX = IMAX(IB)
      Do 10 JA=JAMIN,JAMAX
      G = SINA(JA)-SINA(JA-1)
      If (JA==JAMIN)  G = G - GMIN(JB)
      If (JA==JAMAX)  G = G - GMAX(JB)
      Do 10 IAREV=IAMIN,IAMAX
      IA  = 1 + Mod(IAREV-1,IMA)
      IJA = IA + IMA*(JA-1)
      F   = 1
      If (IAREV==IAMIN)  F = F - FMIN(IB)
      If (IAREV==IAMAX)  F = F - FMAX(IB)
      WEIGHT = WEIGHT + F*G*WTA(IJA)
   10 VALUE  = VALUE  + F*G*WTA(IJA)*A(IJA)
      B(IJB) = DATMIS
      If (WEIGHT /= 0)  B(IJB) = VALUE/WEIGHT
   20 Continue
      Return
      EndSubroutine HNTR4


      Subroutine HNTR4P (WTA,A,B)
!****
!**** HNTR4P is similar to HNTR4 but polar values are replaced by their longitudinal mean.
!**** The 3 Real input values are expected to be Real*4.
!****
      Implicit None
      Real*4  :: WTA(*), A(*), B(*) 
      Integer :: IMA,JMA,IMB,JMB, IMIN,IMAX,JMIN,JMAX
      Real*4  :: DATMIS       
      Real*8  :: SINA,SINB, FMIN,FMAX,GMIN,GMAX
      Common /HNTRCB/ SINA(0:5401),SINB(0:5401), &
                      FMIN(10800),FMAX(10800),GMIN(5401),GMAX(5401), &
                      IMIN(10800),IMAX(10800),JMIN(5401),JMAX(5401), &
                      DATMIS, IMA,JMA,IMB,JMB
!**** Local variables
      Integer :: IB,JB,IJB
      Real*8  :: BMEAN,WEIGHT,VALUE
!****
      Call HNTR4 (WTA,A,B)
!****
!**** Replace individual values near the poles by longitudinal mean
!****
      Do 40 JB=1,JMB,JMB-1
      BMEAN  = DATMIS
      WEIGHT = 0
      VALUE  = 0
      Do 10 IB=1,IMB
      IJB  = IB + IMB*(JB-1)
      If (B(IJB) == DATMIS)  GoTo 20
      WEIGHT = WEIGHT + 1
      VALUE  = VALUE  + B(IJB)
   10 Continue
      If (WEIGHT.ne.0)  BMEAN = VALUE/WEIGHT
   20 Do 30 IB=1,IMB
      IJB  = IB + IMB*(JB-1)
   30 B(IJB) = BMEAN
   40 Continue
      Return
      EndSubroutine HNTR4P
#endif
