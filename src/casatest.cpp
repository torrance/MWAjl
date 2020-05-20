#include <math.h>
#include <casacore/casa/Quanta/MVPosition.h>
#include <casacore/measures/Measures/MPosition.h>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/casa/Quanta.h>
#include <casacore/casa/Quanta/Quantum.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/casa/Quanta/MVEpoch.h>
#include <casacore/casa/Quanta/MVTime.h>
#include <casacore/measures/Measures/MeasConvert.h>
#include <casacore/measures/Measures/MCEpoch.h>
#include <casacore/measures/Measures/MCDirection.h>
#include <iostream>


int main() {
    // Need position to calculate LST
    // MWAPOS = EarthLocation.from_geodetic(lon="116:40:14.93", lat="-26:42:11.95", height=377.8)
    // arrayPos = casacore::MPosition(casacore::MVPosition(-2.55952e+06, 5.09585e+06, -2.84899e+06)); // pos of tile 011
    casacore::MPosition pos(
        casacore::MVPosition(
            casacore::Quantity(1, "km"),
            (116 + 40.0 / 60.0 + 14.93 / 3600.0) * (M_PI / 180),
            -(26 + 42.0 / 60.0 + 11.95 / 3600.0) * (M_PI / 180)
        ),
        casacore::MPosition::ITRF
    );
    std::cout << pos << "\n";

    casacore::MEpoch time(
        casacore::MVEpoch(
            casacore::MVTime(2018, 5, 17, (8 + 18.0 / 60.0) / 24.0)
        ),
        casacore::MEpoch::UTC
    );
    std::cout << time << "\n";

    casacore::MeasFrame frame(pos, time);

    // casacore::MEpoch::Ref sidref(casacore::MEpoch::LAST, frame);
    casacore::MDirection::Ref azelref(casacore::MDirection::AZEL, frame);

    // casacore::MEpoch::Convert tosid(time, sidref);
    // casacore::MEpoch sidtime = tosid();
    // std::cout << sidtime << "\n";

    casacore::MDirection radec(casacore::Quantity(11, "deg"), casacore::Quantity(-30, "deg"));
    casacore::MDirection azel = casacore::MDirection::Convert(radec, azelref)();
    std::cout << "Apparent coordinates: " << azel.toString() << "\n";
    std::cout << "Azimuth : " << azel.getAngle().getValue("rad")[0] << " Altitude: " << azel.getAngle().getValue("rad")[1] << "\n";

}