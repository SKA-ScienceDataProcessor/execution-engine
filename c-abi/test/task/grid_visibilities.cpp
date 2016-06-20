
#include <algorithm>
#include <math.h>
#include <complex>

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

//
//  task_grid_visibilities()
//
//  CJS: 10/08/2015
//
//  Produce gridded visibilities by convolving the complex visibilities with the kernel function.

int task_grid_visibilities(int stage,
                           task_instance inst,
                           int argc,
                           const task_param* argv)
{

    const task_param *wavelengths = task_input(argc, argv, "wavelength", TASK_PARAM_TYPE_DOUBLE, 1);
    const task_param *uvw = task_input(argc, argv, "baselines", TASK_PARAM_TYPE_DOUBLE, 2);
    const task_param *vis = task_input(argc, argv, "visibilities", TASK_PARAM_TYPE_COMPLEX_DOUBLE, 2);
    const task_param *w_planes = task_input(argc, argv, "w_planes", TASK_PARAM_TYPE_DOUBLE, 1);
    const task_param *kernels = task_input(argc, argv, "kernels", TASK_PARAM_TYPE_COMPLEX_DOUBLE, K_COUNT);

    const task_param *grid = task_input_output(argc, argv, "grid", TASK_PARAM_TYPE_COMPLEX_DOUBLE, 2);

    if (!wavelengths || !uvw || !vis || !w_planes || !kernels || !grid) {
        fprintf(stderr, "parameter not found!\n");
        return 1;
    }

    if (strcmp(kernels->regions[K_PLANE].name, "w_plane") ||
        strcmp(kernels->regions[K_OVER_V].name, "over_v") ||
        strcmp(kernels->regions[K_OVER_U].name, "over_u") ||
        strcmp(kernels->regions[K_V].name, "v")) {
        strcmp(kernels->regions[K_U].name, "u") ||

        fprintf(stderr, "Unexpected kernel domains!\n");
        return 1;
    }

    // Get kernel size, must be quadratic and odd
    int64_t pKernelSize = kernels->regions[K_U].extent;
    if (kernels->regions[K_V].extent != pKernelSize ||
        (pKernelSize % 2) != 1) {

        fprintf(stderr, "Unacceptable kernel dimensions!\n");
        return 1;
    }
    int64_t pSupport = (pKernelSize - 1) / 2;

    // Get oversample factor, must be quadratic as well
    int64_t pOversample = kernels->regions[K_OVER_U].extent;
    if (kernels->regions[K_OVER_V].extent != pOversample) {

        fprintf(stderr, "Unacceptable kernel oversampling dimensions!\n");
        return 1;
    }

    // Check number of w-planes
    int64_t _wPlanes = w_planes->regions[K_PLANE].extent;
    if (kernels->regions[K_PLANE].extent != _wPlanes) {

        fprintf(stderr, "Unacceptable kernel plane count!\n");
        return 1;
    }

    // Get number of visibilities
    if (strcmp(uvw->regions[VIS].name, "vis") || strcmp(vis->regions[VIS].name, "vis") ||
        uvw->regions[VIS].extent != vis->regions[VIS].extent) {
        fprintf(stderr, "visibility region mismatch!\n");
        return 1;
    }
    int pNumSamples = uvw->regions[VIS].extent;

    // Check baseline coordinates
    if (strcmp(uvw->regions[BLCOO].name, "bl_coordinate") ||
        uvw->regions[BLCOO].extent < BLCOO_COUNT) {
        fprintf(stderr, "baseline coordinate region mismatch!\n");
        return 1;
    }

    // Check number of channels
    if (strcmp(vis->regions[CHAN].name, "channel") || strcmp(wavelengths->regions[0].name, "channel") ||
        vis->regions[CHAN].extent != wavelengths->regions[0].extent) {

        fprintf(stderr, "channel region mismatch!\n");
        return 1;
    }
    int pNumChannels = vis->regions[CHAN].extent;

    // Check grid size (must be quadratic)
    int _uvPixels = grid->regions[GRIDV].extent;
    if (strcmp(grid->regions[GRIDV].name, "v") || strcmp(grid->regions[GRIDU].name, "u") ||
        grid->regions[GRIDU].extent != _uvPixels) {

        fprintf(stderr, "grid dimensions mismatch!\n");
        return 1;
    }

    // Calculate uv cell size
    double _uvCellSize = (grid->regions[GRIDU].max - grid->regions[GRIDU].min) / _uvPixels;
    double _uvCellSizeV = (grid->regions[GRIDV].max - grid->regions[GRIDV].min) / _uvPixels;
    if (fabs(_uvCellSize / _uvCellSizeV - 1.0) > 1e-15) {

        fprintf(stderr, "grid resolution mismatch!\n");
        return 1;
    }

    if (stage & TASK_STAGE_START) {

        const double INTERVAL = 0.05; // progress is updated at these % intervals.

        // loop through the list of visilibities.  use the w-coordinate to
        // pick an appropriate w-plane use the u and v coordinates to pick
        // an appropriate oversampled kernel.  loop through the kernel
        // pixels, multiplying each one by the complex visibility and then
        // adding to the grid.

        // loop through the samples.
        int fraction = -1;
        for ( long int sample = 0; sample < pNumSamples; sample++ )
        {

            // display progress.
            if ((double)sample / (double)pNumSamples >= ((double)(fraction + 1) * INTERVAL))
            {
                fraction = fraction + 1;
                printf( "%i%%.", (int)(fraction * INTERVAL * 100) );
                fflush( stdout );
            }

            // loop through the channels.
            for ( long int channel = 0; channel < pNumChannels; channel++ )
            {

                // convert uvw coordinates to units of lambda.
                double uvwSampleU = task_f64(uvw, sample, BLCOO_U) / task_f64(wavelengths, channel);
                double uvwSampleV = task_f64(uvw, sample, BLCOO_V) / task_f64(wavelengths, channel);
                // double uvwSampleW = task_f64(uvw, sample, BLCOO_W) / task_f64(wavelengths, channel);

                // get the exact u and v positions within the grid.
                double uExact = uvwSampleU / _uvCellSize;
                double vExact = uvwSampleV / _uvCellSize;

                // find the nearest grid point.
                int uGrid = floor( uExact );
                int vGrid = floor( vExact );

                // choose which oversampled kernel to use in the u and v diretions.
                int uOversample = (int)((uExact - (double)uGrid) * pOversample);
                int vOversample = (int)((vExact - (double)vGrid) * pOversample);

                // translate the u and v coordinates into the centre of the image.
                uGrid = uGrid + (_uvPixels / 2.0);
                vGrid = vGrid + (_uvPixels / 2.0);

                // get the w-plane index. we currently use 0 since all the kernels are the same anyway.
                int wPlane = 0;

                // get the complex visibility.
                complex<double> visibility = task_c64(vis, sample, channel);

                // loop over the kernel in the u and v directions.
                for (int j = max<int64_t>(0, pSupport - vGrid); j < min(pKernelSize, _uvPixels + pSupport - vGrid); j++ )
                    for ( int i = max<int64_t>(0, pSupport - uGrid); i < min(pKernelSize, _uvPixels + pSupport - uGrid); i++ )
                    {

                        // get exact grid coordinates,
                        int iGrid = uGrid + i - pSupport;
                        int jGrid = vGrid + j - pSupport;

                        // get kernel value.
                        complex<double> kernel = task_c64(kernels, wPlane, vOversample, uOversample, j, i);

                        // update the grid.
                        task_c64(grid, jGrid, iGrid) += visibility * kernel;

                    }

            }

        }
        printf("100%%\n\n");

    }

    // return success/failure.
    return 0;
}
