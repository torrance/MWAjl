#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/IPosition.h>
#include <casacore/casa/Arrays/Slicer.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/ms/MeasurementSets.h>
#include <casacore/tables/Tables.h>
#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/Tables/ColumnDesc.h>
#include <casacore/tables/Tables/ScalarColumn.h>
#include <casacore/tables/TaQL/TableParse.h>

#include <iostream>
#include <vector>
#include <array>
#include <chrono>

int main() {
    std::cout << "Starting...\nOpening table...";
    fflush(stdout);
    casacore::Table tbl("~/scratch/MWAjl/1248714872/1248714872.ms");
    std::cout << "Done.\nCreating array column...";
    fflush(stdout);
    casacore::ArrayColumn<casacore::Complex> arrcol(tbl, "DATA");
    std::cout << "Done.\n\n";

    size_t npol = 4;
    size_t nchan = 200;

    casacore::IPosition start(2, 0, 0), end(2, npol - 1, nchan - 1);
    casacore::Slicer slice(start, end, casacore::Slicer::LengthOrLast::endIsLast);

    std::chrono::steady_clock::time_point begint = std::chrono::steady_clock::now();
    std::cout << "Done.\nGetting array...";
    fflush(stdout);
    auto arr = new casacore::Array<casacore::Complex>;
    arrcol.getColumn(slice, *arr);
    std::cout << "Done.\n Get storage...";
    fflush(stdout);
    casacore::Bool deleteIt;
    casacore::Complex* ptr = arr->getStorage(deleteIt);

    if (deleteIt) {
        delete arr;
    }
    std::cout << "Done\nFirst values...";
    fflush(stdout);
    for (size_t i = 0; i < 10;  i++) {
        std::cout << " " << ptr[i];
    }
    std::cout << " Done\n";
    fflush(stdout);
    std::chrono::steady_clock::time_point endt = std::chrono::steady_clock::now();
    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::milliseconds>(endt - begint).count() << "[ms]" << std::endl;



    begint = std::chrono::steady_clock::now();
    std::cout << "Allocate memory...";
    fflush(stdout);
    auto nrow = arrcol.nrow();
    auto ptr1 = (casacore::Complex*) malloc(npol * nchan * nrow * sizeof(casacore::Complex));
    std::cout << "Done.\n Reading and copying...";
    fflush(stdout);
    casacore::Array<casacore::Complex> arr1;
    for (size_t row = 0; row < nrow; row++) {
        arrcol.getSlice(row, slice, arr1);
        for (size_t chan = 0; chan < nchan; chan++) {
            for (size_t pol = 0; pol < npol; pol++) {
                ptr1[chan * nrow * npol + row * npol + pol] = arr1.data()[chan * npol + pol];
            }
        }

    }
    std::cout << "Done\nFirst values...";
    fflush(stdout);
    for (size_t i = 0; i < 10;  i++) {
        std::cout << " " << ptr1[i];
    }
    std::cout << " Done\n";
    endt = std::chrono::steady_clock::now();
    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::milliseconds>(endt - begint).count() << "[ms]" << std::endl;

}