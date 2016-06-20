
#include <complex>

#include "c-abi.h"
#include "c-abi-helpers.h"
#include "tasks.h"

using namespace std;

// Domains
const int IMAGE_L = 1;
const int IMAGE_M = 0;

static double grid_stats(int64_t _uvPixels, const task_param *grid, const char *name);

//  task_hogbom_clean()
//
//  CJS: 03/11/2015
//
//  Perform a Hogbom clean on our dirty image.

int task_hogbom_clean(const int stage,
                      task_instance inst,
                      int argc,
                      const task_param* argv)
{

    const task_param *minor_cycles = task_value(argc, argv, "minor_cycles", TASK_PARAM_TYPE_INT_64);
    const task_param *loop_gain = task_value(argc, argv, "loop_gain", TASK_PARAM_TYPE_DOUBLE);
    const task_param *psf_x = task_value(argc, argv, "psf_x", TASK_PARAM_TYPE_INT_64);
    const task_param *psf_y = task_value(argc, argv, "psf_y", TASK_PARAM_TYPE_INT_64);

    // TODO: Neither of these should actually be complex...
    const task_param *dirty_beam  = task_input(argc, argv, "dirty_beam", TASK_PARAM_TYPE_COMPLEX_DOUBLE, 2);
    const task_param *dirty_image = task_input(argc, argv, "dirty_image", TASK_PARAM_TYPE_COMPLEX_DOUBLE, 2);
    const task_param *clean_beam  = task_input(argc, argv, "clean_beam", TASK_PARAM_TYPE_COMPLEX_DOUBLE, 2);
    const task_param *clean_image = task_output(argc, argv, "clean_image", TASK_PARAM_TYPE_COMPLEX_DOUBLE, 2);

    if (!minor_cycles || !loop_gain || !dirty_beam || !dirty_image || !clean_beam || !clean_image) {
        fprintf(stderr, "parameter not found!\n");
        return 1;
    }

    // Get image size. Must be the same quadratic region on all of them.
    int _uvPixels = clean_image->regions[IMAGE_M].extent;
    if (strcmp(dirty_beam->regions[IMAGE_L].name, "l") ||
        strcmp(dirty_beam->regions[IMAGE_M].name, "m") ||
        strcmp(dirty_image->regions[IMAGE_L].name, "l") ||
        strcmp(dirty_image->regions[IMAGE_M].name, "m") ||
        strcmp(clean_beam->regions[IMAGE_L].name, "l") ||
        strcmp(clean_beam->regions[IMAGE_M].name, "m") ||
        strcmp(clean_image->regions[IMAGE_L].name, "l") ||
        strcmp(clean_image->regions[IMAGE_M].name, "m") ||
        dirty_beam->regions[IMAGE_L].extent != _uvPixels ||
        dirty_beam->regions[IMAGE_M].extent != _uvPixels ||
        dirty_image->regions[IMAGE_L].extent != _uvPixels ||
        dirty_image->regions[IMAGE_M].extent != _uvPixels ||
        clean_beam->regions[IMAGE_L].extent != _uvPixels ||
        clean_beam->regions[IMAGE_M].extent != _uvPixels ||
        clean_image->regions[IMAGE_L].extent != _uvPixels) {

        fprintf(stderr, "Unexpected image dimensions!\n");
        return 1;
    }

    // Clean image needs to be packed
    if (clean_image->regions[IMAGE_L].stride != 1 ||
        clean_image->regions[IMAGE_M].stride != _uvPixels) {

        fprintf(stderr, "Unexpected clean image strides!\n");
        return 1;
    }
    complex<double> *_cleanImage = &task_c64(clean_image, 0, 0);

    // Get parameters
    int64_t _minorCycles = minor_cycles->i64;
    double _loopGain = loop_gain->f64;
    int64_t _psfX = psf_x->i64;
    int64_t _psfY = psf_y->i64;

    if (stage & TASK_STAGE_START) {

        // Initialise clean image
        memset( _cleanImage, 0, _uvPixels * _uvPixels * sizeof( complex<double> ) );

        printf( "\ncleaning....." );

        double noiseLevel = grid_stats(_uvPixels, dirty_image, "mean");

        // loop over each minor cycle.
        for ( int cycle = 0; cycle < _minorCycles; cycle++ )
        {

            printf( "minor cycle %i.....", cycle );

            // identify the pixel with the largest value.
            double maxPixel = -1;
            int maxX = 0, maxY = 0;
            for ( int i = 0; i < _uvPixels; i++ )
                for ( int j = 0; j < _uvPixels; j++ )
                    if ( real( task_c64(dirty_image, i, j) ) > maxPixel )
                    {
                        maxPixel = real( task_c64(dirty_image, i, j) );
                        maxX = j;
                        maxY = i;
                    }
            printf( "position (%i, %i), value %f\n", maxX, maxY, maxPixel );

            // add the clean beam to the model image at the required
            // position and intensity, and subtract from the dirty
            // image. here we loop over all the pixels in the clean
            // beam, and add them one at a time to the clean image.
            for ( int i = 0; i < _uvPixels; i++ )
                for ( int j = 0; j < _uvPixels; j++ )
                {

                    // calculate position in clean image (x,y).
                    int x = maxX - _psfX + j;
                    int y = maxY - _psfY + i;

                    // are we within the image bounds ?
                    if (x >= 0 && x < _uvPixels && y >= 0 && y < _uvPixels)
                    {

                        // add/subtract the psf (scaled).
                        task_c64(clean_image, y, x) += maxPixel * _loopGain * task_c64(clean_beam, i, j);
                        task_c64(dirty_image, y, x) -= maxPixel * _loopGain * task_c64(dirty_beam, i, j);

                    }

                }

            // check to see if we've reached noise level (2-sigma). if so, exit loop.
            if (maxPixel <= noiseLevel)
            {
                printf( "noise level reached. cleaning will stop.\n" );
                break;
            }
        }

        // add the residuals to the clean image.
        for ( int i = 0; i < _uvPixels; i++ )
            for ( int j = 0; j < _uvPixels; j++ )
                task_c64(clean_image, i, j) += task_c64(dirty_image, i, j);

        // show grid stats for residuals
        grid_stats(_uvPixels, dirty_image, "mean of residuals");

    }
    return 0;
}

double grid_stats(int64_t _uvPixels, const task_param *image, const char *name)
{

    // calculate the average pixel value and standard deviation
    double total = 0, total_sqr = 0;
    for ( int i = 0; i < _uvPixels; i++ )
        for ( int j = 0; i < _uvPixels; i++ ) {
            double r = real( task_c64(image, i, j) );
            total += r;
            total_sqr += r*r;
        }
    double pixels = _uvPixels * _uvPixels;
    double averageValue = total / pixels;
    double standardDeviation = sqrt( total_sqr / pixels - averageValue*averageValue );

    printf( "%s %f, s.d. %f\n", name, averageValue, standardDeviation );

    // Return noise level (2-sigma)
    return averageValue + 2 * standardDeviation;
}
