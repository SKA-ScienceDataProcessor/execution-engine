
#include <math.h>

#include "c-abi.h"
#include "c-abi-helpers.h"
#include "tasks.h"

// Domains
const int VIS = 0;
const int BLCOO = 1;
const int CHAN = 1;

// Coordinates
const int BLCOO_U = 0;
const int BLCOO_V = 1;
const int BLCOO_W = 2;
const int BLCOO_COUNT = 3;

// Chris' phase correction code.
#include "phasecorrection.h"

int task_phase_correction(int stage,
                          task_instance inst,
                          int argc,
                          const task_param* argv)
{

    // Locate and check parameters
    const task_param *in_long  = task_value(argc, argv, "in_longitude",    TASK_PARAM_TYPE_DOUBLE);
    const task_param *in_lat   = task_value(argc, argv, "in_latitude",     TASK_PARAM_TYPE_DOUBLE);
    const task_param *out_long = task_value(argc, argv, "out_longitude",   TASK_PARAM_TYPE_DOUBLE);
    const task_param *out_lat  = task_value(argc, argv, "out_latitude",    TASK_PARAM_TYPE_DOUBLE);

    const task_param *wavelengths = task_input(argc, argv, "wavelength",   TASK_PARAM_TYPE_DOUBLE, 1);

    const task_param *uvw = task_input_output(argc, argv, "baselines",    TASK_PARAM_TYPE_DOUBLE, 2);
    const task_param *vis = task_input_output(argc, argv, "visibilities", TASK_PARAM_TYPE_COMPLEX_DOUBLE, 2);

    if (!in_long || !in_lat || !out_long || !out_lat || !wavelengths || !uvw || !vis) {
        fprintf(stderr, "parameter not found!\n");
        return 1;
    }

    // Get number of visibilities
    if (strcmp(uvw->regions[VIS].name, "vis") || strcmp(vis->regions[VIS].name, "vis") ||
        uvw->regions[VIS].extent != vis->regions[VIS].extent) {
        fprintf(stderr, "visibility region mismatch!\n");
        return 1;
    }
    int _numSamples = uvw->regions[VIS].extent;

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
    int _numChannels = vis->regions[CHAN].extent;

    // Prepare our matrices
    PhaseCorrection *phaseCorrection = reinterpret_cast<PhaseCorrection *>(inst.user);
    if (stage & TASK_STAGE_PREPARE) {

        // Initialise phase correction
        delete phaseCorrection;
        phaseCorrection = new PhaseCorrection();
        inst.user = reinterpret_cast<void *>(phaseCorrection);

        // set up the coordinate system.
        const char J2000[] = "J2000";

        phaseCorrection->inCoords.longitude = in_long->f64;
        phaseCorrection->inCoords.latitude = in_lat->f64;
        strcpy( phaseCorrection->inCoords.epoch, J2000 );

        phaseCorrection->outCoords.longitude = out_long->f64;
        phaseCorrection->outCoords.latitude = out_lat->f64;
        strcpy( phaseCorrection->outCoords.epoch, J2000 );

        phaseCorrection->uvProjection = false;

        // initialise phase rotation matrices.
        phaseCorrection->init();

    }

    if (stage & TASK_STAGE_START) {
        if (!phaseCorrection) { return 1; }

        const double INTERVAL = 0.05; // progress is updated at these % intervals.
        int fraction = -1;
        printf( "\nphase rotating visibilities.....\n" );

        // loop through all the visibilities.
        for ( int i = 0; i < _numSamples; i++ )
        {

            // Rotate uvw coordinate.
            phaseCorrection->uvwIn.x = task_f64(uvw, i, BLCOO_U);
            phaseCorrection->uvwIn.y = task_f64(uvw, i, BLCOO_V);
            phaseCorrection->uvwIn.z = task_f64(uvw, i, BLCOO_W);
            phaseCorrection->rotate();
            task_f64(uvw, i, BLCOO_U) = phaseCorrection->uvwOut.x;
            task_f64(uvw, i, BLCOO_V) = phaseCorrection->uvwOut.y;
            task_f64(uvw, i, BLCOO_W) = phaseCorrection->uvwOut.z;

            // Phase correction to apply to visibilities
            double phase = phaseCorrection->phase;

            // loop through all the channels.
            for ( int j = 0; j < _numChannels; j++ )
            {

                double wavelength = task_f64(wavelengths, j);

                // calculate the phasor.
                const double PI = 3.141592654;
                std::complex<double> phasor(
                    cos( 2 * PI * phase / wavelength ),
                    sin( 2 * PI * phase / wavelength ) );

                // multiply phasor by visibility.
                task_c64(vis, i, j) *= phasor;

            }

            // display progress.
            if ((double)i / (double)_numSamples >= ((double)(fraction + 1) * INTERVAL))
            {
                fraction = fraction + 1;
                printf( "%i%%.", (int)(fraction * INTERVAL * 100) );
                fflush( stdout );
            }

        }
        printf( "100%%\n" );

    }

    if (stage & TASK_STAGE_FREE) {
        delete phaseCorrection;
        inst.user = 0;
    }

    return 0;
}
