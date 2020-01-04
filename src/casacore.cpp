// g++ -shared -fPIC -lcasa_casa -lcasa_tables -lcasa_measures casacore.cpp -o libcasacorejl.so

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


enum errors {
    OK,
    TableNoFile,
    ArraySlicerError,
    TableError,
};

template <typename T>
T* getColumn(casacore::Table* tbl, char* name, int* ndim, size_t** shape, int sliceLength, size_t* blc, size_t* trc, int* const error) {
    try {
        casacore::Array<T> array;

        auto table_desc = tbl->tableDesc();
        auto column_desc = table_desc.columnDesc(name);
        if (column_desc.isScalar()) {
            casacore::ScalarColumn<T> column(*tbl, name);
            array = column.getColumn();
        }
        else {
            casacore::ArrayColumn<T> column(*tbl, name);

            // Only array types can be sliced
            if (sliceLength) {
                casacore::IPosition start(sliceLength), end(sliceLength);
                for (uint i = 0; i < sliceLength; i++) {
                    start[i] = blc[i];
                    end[i] = trc[i];
                }
                casacore::Slicer slice(start, end, casacore::Slicer::LengthOrLast::endIsLast);
                array = column.getColumn(slice);
            }
            else {
                array = column.getColumn();
            }
        }

        // Populate shape
        *ndim = array.ndim();
        *shape = (size_t*) malloc(*ndim * sizeof(size_t));
        memcpy(*shape, array.shape().storage(), *ndim * sizeof(size_t));

        // Get raw storage data
        casacore::Bool deleteIt;
        T* data = array.getStorage(deleteIt);
        size_t nelements = array.nelements();
        // Unless deleteIt is true, we need to copy the data
        if (deleteIt) {
            return data;
        }
        else {
            T* mydata = (T*) malloc(nelements * sizeof(T));
            memcpy(mydata, data, nelements * sizeof(T));
            return mydata;
        }
    }
    catch (casacore::ArraySlicerError e) {
        *error = ArraySlicerError;
        return 0;
    }
}

extern "C" {
    casacore::Table* table_open(char* path, int* const error) {
        *error = OK;
        try {
            return new casacore::Table(path);
        }
        catch (const casacore::TableNoFile& e) {
            *error = TableNoFile;
            return 0;
        }
    }

    void table_close(casacore::Table* tbl) {
        delete tbl;
    }

    casacore::Table* taql(casacore::Table* tbl, char* query) {
        auto res = tableCommand(query, *tbl).table();
        return new casacore::Table(res);
    }

    bool column_exists(casacore::Table* tbl, char* name) {
        return tbl->tableDesc().isColumn(name);
    }

    int column_type(casacore::Table* tbl, char* name) {
        casacore::ROTableColumn col(*tbl, name);
        return col.columnDesc().dataType();
    }

    bool column_is_fixed_shape(casacore::Table* tbl, char* name) {
        casacore::ROTableColumn col(*tbl, name);
        return (col.columnDesc().options() & casacore::ColumnDesc::FixedShape) == casacore::ColumnDesc::FixedShape;
    }

    size_t* column_info(casacore::Table* tbl, char* name, int* element_type, int* dimension) {
        casacore::ROTableColumn col(*tbl, name);
        *element_type = col.columnDesc().dataType();
        if (col.columnDesc().isScalar()) {
            *dimension = 1;
            size_t* shape = new size_t[1];
            shape[0] = tbl->nrow();
            return shape;
        }
        else {
            if (column_is_fixed_shape(tbl, name)) {
                // for fixed shape columns we can use col.shapeColumn() to get the shape
                auto colshape = col.shapeColumn();
                *dimension = colshape.size() + 1;
                size_t* shape = new size_t[*dimension];
                for (uint i = 0; i < colshape.size(); ++i) {
                    shape[i] = colshape[i];
                }
                shape[*dimension - 1] = tbl->nrow();
                return shape;
            }
            else {
                // if the column doesn't have a fixed shape, we will rely on the shape of the array
                // in the first row to get the shape of the entire column, but we need to be sure
                // that the first row actually has an array in it (ie. it's not left undefined).
                if (col.isDefined(0)) {
                    *dimension = col.ndim(0) + 1;
                    size_t* shape = new size_t[*dimension];
                    auto colshape0 = col.shape(0);
                    for (uint i = 0; i < colshape0.size(); ++i) {
                        shape[i] = colshape0[i];
                    }
                    shape[*dimension - 1] = tbl->nrow();
                    return shape;
                }
                else {
                    *dimension = 1;
                    size_t* shape = new size_t[1];
                    shape[0] = tbl->nrow();
                    return shape;
                }
            }
        }
    }

    casacore::Bool* get_column_boolean(casacore::Table* tbl, char* name, int* ndim, size_t** shape, int sliceLength, size_t* start, size_t* end, int* const error) {
        return getColumn<casacore::Bool>(tbl, name, ndim, shape, sliceLength, start, end, error);
    }

    casacore::Int* get_column_int(casacore::Table* tbl, char* name, int* ndim, size_t** shape, int sliceLength, size_t* start, size_t* end, int* const error) {
        return getColumn<casacore::Int>(tbl, name, ndim, shape, sliceLength, start, end, error);
    }

    casacore::Float* get_column_float(casacore::Table* tbl, char* name, int* ndim, size_t** shape, int sliceLength, size_t* start, size_t* end, int* const error) {
        return getColumn<casacore::Float>(tbl, name, ndim, shape, sliceLength, start, end, error);
    }

    casacore::Double* get_column_double(casacore::Table* tbl, char* name, int* ndim, size_t** shape, int sliceLength, size_t* start, size_t* end, int* const error) {
        return getColumn<casacore::Double>(tbl, name, ndim, shape, sliceLength, start, end, error);
    }

    casacore::Complex* get_column_complex(casacore::Table* tbl, char* name, int* ndim, size_t** shape, int sliceLength, size_t* start, size_t* end, int* const error) {
        return getColumn<casacore::Complex>(tbl, name, ndim, shape, sliceLength, start, end, error);
    }

    // mjd :: Modified Julian Date (UTC)
    // lon, lat :: longitude and lattitude (radians) of observing site (ITRF)
    casacore::MeasFrame* frame_new(double mjd, double lon, double lat) {
        casacore::MPosition pos(
            casacore::MVPosition(casacore::Quantity(1, "km"), lon, lat),
            casacore::MPosition::ITRF
        );
        casacore::MEpoch time(
            casacore::MVEpoch(casacore::MVEpoch(casacore::Quantity(mjd, "d"))),
            casacore::MEpoch::UTC
        );
        return new casacore::MeasFrame(pos, time);
    }

    void frame_del(casacore::MeasFrame* frame) {
        delete frame;
    }

    void radec_to_altaz(double ra, double dec, casacore::MeasFrame* frame, double* alt, double* az) {
        casacore::MDirection::Ref azelref(casacore::MDirection::AZELNE, *frame);
        casacore::MDirection radec(casacore::Quantity(ra, "rad"), casacore::Quantity(dec, "rad"));
        casacore::MDirection azel = casacore::MDirection::Convert(radec, azelref)();
        auto vals = azel.getAngle().getValue("rad"); // [az, el]
        *alt = vals[1];
        *az = vals[0];
    }
}