
#include <complex>

#include "c-abi.h"
#include "c-abi-helpers.h"
#include "tasks.h"

using namespace std;

// Domains
const int IMAGE_L = 1;
const int IMAGE_M = 0;

struct Params
{
    double angle;
    double r1;
    double r2;
    double x;
    double y;
};

//
//	gaussian2D()
//
//	CJS: 05/11/2015
//
//	Create an elliptical 2D Gaussian at position (pXPos, pYPos), with
//	long and short axes pR1 and pR2, rotated at pAngle.
//

static double gaussian2D( double pX, double pY, double pAngle, double pR1, double pR2, double pXPos, double pYPos )
{

	double iOffset = pX - pXPos;
	double jOffset = pY - pYPos;
	double rLong = ((jOffset * cos( pAngle )) + (iOffset * sin( pAngle ))) / pR1;
	double rShort = ((iOffset * cos( pAngle )) - (jOffset * sin( pAngle ))) / pR2;

	return exp( -pow( rLong, 2 ) - pow( rShort, 2 ) );

} // gaussian2D

//  task_hogbom_clean()
//
//	CJS: 04/11/2015
//
//	Generate the clean beam by fitting elliptical Gaussian to the dirty beam.

int task_generate_clean_beam(const int stage,
                             task_instance inst,
                             int argc,
                             const task_param* argv)
{
    const task_param *max_size_for_psf_fitting = task_value(argc, argv, "max_size_for_psf_fitting", TASK_PARAM_TYPE_INT_64);
    const task_param *dirty_beam  = task_input(argc, argv, "dirty_beam", TASK_PARAM_TYPE_COMPLEX_DOUBLE, 2);

    const task_param *clean_beam  = task_output(argc, argv, "clean_beam", TASK_PARAM_TYPE_COMPLEX_DOUBLE, 2);
    const task_param *psf_x = task_output(argc, argv, "psf_x", TASK_PARAM_TYPE_INT_64, 0);
    const task_param *psf_y = task_output(argc, argv, "psf_y", TASK_PARAM_TYPE_INT_64, 0);

    if (!dirty_beam || !clean_beam || !psf_x || !psf_y) {
        fprintf(stderr, "parameter not found!\n");
        return 1;
    }

    // Get image size. Must be the same quadratic region
    int _uvPixels = clean_beam->regions[IMAGE_M].extent;
    if (strcmp(dirty_beam->regions[IMAGE_L].name, "l") ||
        strcmp(dirty_beam->regions[IMAGE_M].name, "m") ||
        strcmp(clean_beam->regions[IMAGE_L].name, "l") ||
        strcmp(clean_beam->regions[IMAGE_M].name, "m") ||
        dirty_beam->regions[IMAGE_L].extent != _uvPixels ||
        dirty_beam->regions[IMAGE_M].extent != _uvPixels ||
        clean_beam->regions[IMAGE_L].extent != _uvPixels) {

        fprintf(stderr, "Unexpected beam image dimensions!\n");
        return 1;
    }

    // get parameter values
    int64_t MAX_SIZE_FOR_PSF_FITTING = max_size_for_psf_fitting->i64;

    if (stage & TASK_STAGE_START) {

        struct Params bestFit;
        bestFit.angle = 0;
        bestFit.r1 = 0;
        bestFit.r2 = 0;
        bestFit.x = ((double)_uvPixels / 2);
        bestFit.y = ((double)_uvPixels / 2) - 1;

        struct Params testFit;
        double r = 0;
        double bestError = -1;

        // initialise random seed.
        srand( time( NULL ) );

        // find the size of the region we need to fit.
        for ( int i = 0; i < _uvPixels; i++ )
            for ( int j = 0; j < _uvPixels; j++ )
                if ( abs( task_c64(dirty_beam, i, j) ) >= 0.1 )
                {
                    double rTemp = pow((double)j - bestFit.x, 2) + pow((double)i - bestFit.y, 2);
                    if (rTemp > r)
                        r = rTemp;
                }
        r = sqrt( r );
        bestFit.r1 = r / 2;
        bestFit.r2 = r / 2;

        // ensure psf fitting region isn't huge.
        if (r > MAX_SIZE_FOR_PSF_FITTING)
            r = MAX_SIZE_FOR_PSF_FITTING;

        printf("size of ellipse fitting region - %f\n", r);

        bool updatedAnything;
        do
        {

            // reset flag.
            updatedAnything = false;

            // the scale loop will attempt to fit the parameters at four different scales. the scale of each change
            // divides by two with each iteration of the loop.
            struct Params scaleFit;
            scaleFit.angle = 1;
            scaleFit.r1 = r / 2;
            scaleFit.r2 = r / 2;
            scaleFit.x = 1;
            scaleFit.y = 1;
            for ( int scale = 0; scale < 6; scale++ )
            {

                // the param loop will change each of the 5 fit parameters in turn.
                // param 0 = r1, param 1 = angle, param 2 = r2, param 3 = x, param 4 = y.
                for ( int param = 0; param < 5; param++ )

                    // attempt 100 random changes to this parameter.
                    for ( int test = 0; test < 40; test++ )
                    {

                        // randomly add or subtract some amount from this parameter, and re-assess the error.
                        // here we produce a random double between -0.5 and 0.5.
                        double change = (static_cast <double> (rand()) / static_cast <double> (RAND_MAX)) - 0.5;

                        testFit = bestFit;
                        switch (param)
                        {
                            case 0:
                                testFit.r1 = bestFit.r1 + (change * scaleFit.r1);
                                if (testFit.r1 < 0)
                                    testFit.r1 = 0;
                                break;
                            case 1:
                                testFit.angle = bestFit.angle + (change * scaleFit.angle);
                                break;
                            case 2:
                                testFit.r2 = bestFit.r2 + (change * scaleFit.r2);
                                if (testFit.r2 < 0)
                                    testFit.r2 = 0;
                                break;
                            case 3:
                                testFit.x = bestFit.x + (change * scaleFit.x);
                                break;
                            case 4:
                                testFit.y = bestFit.y + (change * scaleFit.y);
                                break;
                        }

                        // get error between the actual dirty beam and the fitted ellipse.
                        double error = 0;
                        for ( int i = (int)floor( testFit.x - r ); i < (int)ceil( testFit.x + r ); i++ )
                            for ( int j = testFit.y - ceil( r ); j < testFit.y + ceil( r ); j++ )
                                error = error + pow( abs( task_c64(dirty_beam, j, i) ) -
                                        gaussian2D( (double)i, (double)j, testFit.angle, testFit.r1, testFit.r2, testFit.x, testFit.y ), 2 );

                        // are these parameters an improvement?
                        if (error < bestError || bestError == -1)
                        {
                            bestFit = testFit;
                            bestError = error;
                            updatedAnything = true;
                        }

                    }

                // divide all the scales by 2 so that we are fitting finer values.
                scaleFit.angle = scaleFit.angle / 2;
                scaleFit.r1 = scaleFit.r1 / 2;
                scaleFit.r2 = scaleFit.r2 / 2;
                scaleFit.x = scaleFit.x / 2;
                scaleFit.y = scaleFit.y / 2;

            }

        } while (updatedAnything == true);

        // TEMP: fix the beam for LOFAR.
        //bestFit.angle = 0;
        //bestFit.r1 = 35;
        //bestFit.r2 = 35;

        printf( "parameters fitted to psf: angle %f, r1 %f, r2 %f, x %f, y %f\n", bestFit.angle, bestFit.r1, bestFit.r2, bestFit.x, bestFit.y );

        // populate Gaussian image.
        for ( int i = 0; i < _uvPixels; i++ )
            for ( int j = 0; j < _uvPixels; j++ )
                task_c64(clean_beam, i, j) =
                    gaussian2D( (double)j, (double)i, bestFit.angle, bestFit.r1,
                                bestFit.r2, bestFit.x, bestFit.y );

        // make a note of the centre position of the psf - we'll need to for cleaning.
        task_f64(psf_x) = (int)floor( bestFit.x + 0.5 );
        task_f64(psf_y) = (int)floor( bestFit.y + 0.5 );

    }

    return 0;
}
