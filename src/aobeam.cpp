// g++ -O3 -shared -fPIC -I ~/Downloads/mwa-reduce/beam -I/usr/include/hdf5/serial src/aobeam.cpp ~/Downloads/mwa-reduce/beam/beam2016implementation.cpp -lhdf5_cpp -lboost_filesystem -o src/libaobeam.so

#include <beam2016implementation.h>
#include <complex>
#include <math.h>

extern "C" {
    Beam2016Implementation* beam_new(double* delays, double* amps, char* path) {
        return new Beam2016Implementation(delays, amps, path);
    }

    void beam_del(Beam2016Implementation* beam) {
        delete beam;
    }

    void beamjones(Beam2016Implementation* beam, int freq, size_t N, double* alts, double* azs, std::complex<double>* jones) {
        for (size_t i = 0; i < N; i++) {
            double az_deg = azs[i] * (180./M_PI);
            double za_deg = 90. - alts[i] * (180./M_PI);
            auto j = beam->CalcJones(az_deg, za_deg, freq, true);

            jones[4 * i + 0] = j.j00;
            jones[4 * i + 1] = j.j01;
            jones[4 * i + 2] = j.j10;
            jones[4 * i + 3] = j.j11;
        }
    }

    int find_closest_freq(Beam2016Implementation* beam, int freq) {
        return beam->find_closest_freq(freq);
    }
}