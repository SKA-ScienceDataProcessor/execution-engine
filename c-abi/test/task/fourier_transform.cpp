
#include <complex>
#include <fftw3.h>

#include "c-abi.h"
#include "c-abi-helpers.h"
#include "tasks.h"

using namespace std;

// Domains
const int VIS = 0;
const int BLCOO = 1;
const int CHAN = 1;

const int GRIDV = 0;
const int GRIDU = 1;

// Coordinates
const int BLCOO_U = 0;
const int BLCOO_V = 1;
const int BLCOO_W = 2;
const int BLCOO_COUNT = 3;

struct fft_data
{
    fftw_plan plan;
    fftw_complex *grid;
};

//
//  task_fourier_transform()
//
//  CJS: 11/08/2015
//
//  Make a dirty image by inverse FFTing the gridded visibilites.

int task_fourier_transform(int stage,
                           task_instance inst,
                           int argc,
                           const task_param* argv)
{
    const task_param *grid = task_input_output(argc, argv, "grid", TASK_PARAM_TYPE_COMPLEX_DOUBLE, 2);

    // Check grid size (must be quadratic)
    int pSize = grid->regions[GRIDV].extent;
    if (strcmp(grid->regions[GRIDV].name, "v") || strcmp(grid->regions[GRIDU].name, "u") ||
        grid->regions[GRIDU].extent != pSize) {

        fprintf(stderr, "grid dimensions mismatch!\n");
        return 1;
    }

    // Check grid packing. FFTW seems to require that there is no padding?
    if (grid->regions[GRIDU].stride != 1 || grid->regions[GRIDV].stride != pSize) {

        fprintf(stderr, "grid stride mismatch!\n");
        return 1;
    }

    // Prepare plan
    fft_data *user = reinterpret_cast<fft_data *>(inst.user);
    if (stage & TASK_STAGE_PREPARE) {

        // Allocate user data
        if (!user) {
            user = new fft_data();
            inst.user = reinterpret_cast<void *>(user);
        }

        // Make plan
        user->grid = (fftw_complex *) fftw_malloc( pSize * pSize * sizeof( fftw_complex ) );
        user->plan = fftw_plan_dft_2d( pSize, pSize, user->grid, user->grid, FFTW_BACKWARD, FFTW_MEASURE );
    }

    if (stage & TASK_STAGE_START) {

        // Copy the image from the original memory location to the
        // FFTW memory location. execute the fft.
        complex<double> *pGrid = &task_c64(grid, 0, 0);
        memcpy( user->grid, grid->data, pSize * pSize * sizeof( fftw_complex ) );
        fftw_execute( user->plan );

        // Do FFT shift.
        for ( int i = 0; i < pSize / 2; i++ )
            for ( int j = 0; j < pSize / 2; j++ )
            {
                memcpy( &pGrid[ (j * pSize) + i ],
                        &user->grid[ ((j + (pSize / 2)) * pSize) + i + (pSize / 2) ],
                        sizeof( complex<double> ) );
                memcpy( &pGrid[ (j * pSize) + i + (pSize / 2) ],
                        &user->grid[ ((j + (pSize / 2)) * pSize) + i ],
                        sizeof( complex<double> ) );
                memcpy( &pGrid[ ((j + (pSize / 2)) * pSize) + i ],
                        &user->grid[ (j * pSize) + i + (pSize / 2) ],
                        sizeof( complex<double> ) );
                memcpy( &pGrid[ ((j + (pSize / 2)) * pSize) + i + (pSize / 2) ],
                        &user->grid[ (j * pSize) + i ],
                        sizeof( complex<double> ) );
            }

    }

    if (stage & TASK_STAGE_FREE) {

        // destroy the FFT plan.
        fftw_destroy_plan( user->plan );
        fftw_cleanup();

        // Free FFTW memory
        fftw_free( user->grid );

        delete user;
        inst.user = 0;
    }

    return 0;
}
