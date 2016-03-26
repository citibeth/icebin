#pragma once

namespace icebin {

/** Base class for per-ice-sheet parameters provided by the GCM, for
the coupler.  The TYPE of these parameters depends on the GCM, and
there is one INSTANCE per ice sheet. */
class GCMPerIceSheetParams {
public:
    virtual ~GCMPerIceSheetParams() {}
};

}
