#include <ibmisc/rlarray.hpp>
#include <icebin/RLWeightedSparse.hpp>

using namespace ibmisc;

namespace icebin {


void RLWeightedSparse::ncio(NcIO &ncio, std::string const &vname)
{
    auto info_v = get_or_add_var(ncio, vname + ".info", "int", {});
    get_or_put_att(info_v, ncio.rw, "conservative", conservative);

    wM.ncio(ncio, vname + ".wM", "int", "int", "double");
    M.ncio(ncio, vname + ".M", "int", "int", "double");
    Mw.ncio(ncio, vname + ".Mw", "int", "int", "double");
}

#if 0
void RLWeightedSparse::apply(
    // this = BvA
    blitz::Array<double,2> const &A,      //  IN: A{nj} one row per variable
    blitz::Array<double,2> &B,            // OUT: B{ni} one row per variable
    SparseFillType fill_type,
    bool force_conservation) const
{
    // Clear the output
    switch(fill_type) {
        case SparseFillType::zero :
            B = 0;
        break;
        case SparseFillType::nan :
            B = std::numeric_limits<double>::quiet_NaN();
        break;
    }

    // Multiply
    for (auto gen(M.generator()); ++gen; ) {
        auto const iB(gen.index(0));
        auto const iA(gen.index(1));
        double const val = A(iA) * gen.value();


        switch(fill_type) {
            case SparseFillType::zero :
                B(iB) += val;
            break;
            case SparseFillType::nan :
                if (std::isnan(B(iB))) B(iB) = val;
                else B(iB) += val;
            break;
        }
    }

    // Adjust for conservation
    if (force_conservation && !conservative) {

        // Determine original mass
        for (auto gen(Mw.generator()); ++gen; ) {
            auto const iA(gen.index(0));
            auto const Mw_iA(gen.value());
            mass_A += Mw_iA * A(iA);
        }

        // Determine final mass
        for (auto gen(wM.generator()); ++gen; ) {
            auto const iB(gen.index(0));
            auto const wM_iB(gen.value());
            mass_B += wM_iB * B(iB);
        }

        // Determine adjustment factor for B
        double const factor = mass_A / mass_B;

        // Multiply by the adjustment factor
        // Determine final mass
        for (auto gen(wM.generator()); ++gen; ) {
            auto const iB(gen.index(0));
            B(iB) *= factor;
        }
    }
}
#endif

}    // namespace icebin
