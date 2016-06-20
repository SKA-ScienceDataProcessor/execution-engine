
#include "c-abi.h"
#include "c-abi-helpers.h"
#include "tasks.h"

#include <complex.h>

//
//	task_generate_w_kernel()
//
//	CJS: 16/12/2015
//	PMW: 16/06/2016
//
//	Generates the convolution function. Currently just a single pixel at the centre of the kernel.
//

int task_generate_w_kernel(int stage,
                           task_instance inst,
                           int argc,
                           const task_param* argv)
{

    // Locate and check parameters
    const task_param *far_field   = task_value(argc, argv, "far_field",    TASK_PARAM_TYPE_INT_64);
    const task_param *theta       = task_value(argc, argv, "theta",        TASK_PARAM_TYPE_DOUBLE);

    const task_param *w_planes    = task_input(argc, argv, "w_planes",     TASK_PARAM_TYPE_DOUBLE, 1);

    const task_param *kernels     = task_output(argc, argv, "kernels",     TASK_PARAM_TYPE_COMPLEX_DOUBLE, K_COUNT);

    if (!far_field || !theta || !w_planes || !kernels) {
        fprintf(stderr, "Required parameter not found!\n");
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
    int64_t _kernelSize = kernels->regions[K_U].extent;
    if (kernels->regions[K_V].extent != _kernelSize ||
        (_kernelSize % 2) != 1) {

        fprintf(stderr, "Unacceptable kernel dimensions!\n");
        return 1;
    }
    int64_t _support = (_kernelSize - 1) / 2;

    // Kernels must have packed inner array dimensions
    if (kernels->regions[K_V].stride != _kernelSize ||
        kernels->regions[K_U].stride != 1) {

        fprintf(stderr, "Unacceptable inner kernel strides!\n");
        return 1;
    }

    // Get oversample factor, must be quadratic as well
    int64_t _oversample = kernels->regions[K_OVER_U].extent;
    if (kernels->regions[K_OVER_V].extent != _oversample) {

        fprintf(stderr, "Unacceptable kernel oversampling dimensions!\n");
        return 1;
    }

    // Check number of w-planes
    int64_t _wPlanes = w_planes->regions[K_PLANE].extent;
    if (kernels->regions[K_PLANE].extent != _wPlanes) {

        fprintf(stderr, "Unacceptable kernel plane count!\n");
        return 1;
    }

    if (stage & TASK_STAGE_START) {

        // calculate separate kernels for each w-value.
        for ( int w_plane = 0; w_plane < _wPlanes; w_plane++ )

            // calculate separate kernels for each (oversampled) intermediate grid position.
            for ( int oversampleI = 0; oversampleI < _oversample; oversampleI++ )
                for ( int oversampleJ = 0; oversampleJ < _oversample; oversampleJ++ )
                {

                    // clear aa-kernel.
                    memset( &task_c64(kernels, w_plane, oversampleI, oversampleJ, 0, 0), 0,
                            _kernelSize * _kernelSize * sizeof( std::complex<double> ) );

                    // update the central pixel to 1.
                    task_c64(kernels, w_plane, oversampleI, oversampleJ, _support, _support)
                        = std::complex<double>( 1, 0 );

                }

    }

    return 0;
}
