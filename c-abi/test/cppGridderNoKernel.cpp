#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <complex>
#include <math.h>
#include <string.h>
#include <fftw3.h>
#include <time.h>

#include "c-abi.h"
#include "c-abi-helpers.h"
#include "tasks.h"

using namespace std;

//
//  FORWARDS DECLARATIONS
//

bool saveBitmap( char * pOutputFilename );
void performFFT( complex<double> * pGrid, int pSize, bool pInverse, bool pConvertToAbsolute, bool pFFTShift );

//
//  CONSTANTS
//

// the input parameters from file gridder-params.
const char CELL_SIZE[] = "cell_size:";
const char PIXELS_UV[] = "pixels_uv:";
const char W_PLANES[] = "w_planes:";
const char OVERSAMPLE[] = "oversample:";
const char MINOR_CYCLES[] = "minor_cycles:";
const char LOOP_GAIN[] = "loop_gain:";
const char INPUT_RA[] = "input-ra:";
const char INPUT_DEC[] = "input-dec:";
const char OUTPUT_RA[] = "output-ra:";
const char OUTPUT_DEC[] = "output-dec:";

// speed of light.
const long int CONST_C = 299792458;
const double PI = 3.141592654;

// bitmap file header positions.
const int BIT_CONST = 0x00;
const int MAP_CONST = 0x01;
const int IMAGE_SIZE = 0x02;
const int RESERVED = 0x06;
const int FILE_HEADER_SIZE = 0x0A;
const int BITMAP_INFO_HEADER = 0x0E;
const int IMAGE_WIDTH = 0x12;
const int IMAGE_HEIGHT = 0x16;
const int COLOUR_PLANES = 0x1A;
const int BIT_COUNT = 0x1C;
const int COMPRESSION_TYPE = 0x1E;
const int COLOURS_USED = 0x2E;
const int SIGNIFICANT_COLOURS = 0x32;

// other constants.
const int MAX_SIZE_FOR_PSF_FITTING = 40;

//
//  STRUCTURES
//

// vector with floats. Can be used either as a 2 or 3 element vector.
struct VectorF
{
    double u;
    double v;
    double w;
};
typedef struct VectorF VectorF;

//
//  GLOBAL VARIABLES
//

// grid parameters.
double _cellSize = 0;       // the angular size of each output pixel
double _uvCellSize = 0;     // in units of lambda
int _uvPixels = 0;
int _wPlanes = 0;
double _inputRA = 0;
double _inputDEC = 0;
double _outputRA = 0;
double _outputDEC = 0;

// samples.
int _numSamples = 0;
VectorF * _sample = NULL;

// channels.
int _numChannels = 0;
double * _wavelength = NULL;

// w-plane details.
double * _wPlaneMean = NULL;
double * _wPlaneMax = NULL;

// hogbom parameters.
int _minorCycles = 10;
double _loopGain = 0.1;

// visibilities.
complex<double> * _visibility = NULL;

// kernel parameters.
double _oversample = 0;
int _support = 0;
int _kernelSize = 0;

// anti-aliasing kernel parameters.
int _aaSupport = 0;
int _aaKernelSize = 0;

// w-kernel parameters.
int _wSupport = 0;
int _wKernelSize = 0;

// kernel data.
complex<double> * _kernel = NULL;

// gridded data.
complex<double> * _grid = NULL;

// psf.
complex<double> * _psf = NULL;
int64_t _psfX = 0, _psfY = 0;

// clean beam.
complex<double> * _cleanBeam = NULL;

// clean image.
complex<double> * _cleanImage = NULL;

//
//  GENERAL FUNCTIONS
//

//
//  DEVICE FUNCTIONS
//

//
//  HOST FUNCTIONS
//

//
//  gaussian2D()
//
//  CJS: 05/11/2015
//
//  Create an elliptical 2D Gaussian at position (pXPos, pYPos), with long and short axes pR1 and pR2, rotated at pAngle.
//

double gaussian2D( double pX, double pY, double pAngle, double pR1, double pR2, double pXPos, double pYPos )
{

    double iOffset = pX - pXPos;
    double jOffset = pY - pYPos;
    double rLong = ((jOffset * cos( pAngle )) + (iOffset * sin( pAngle ))) / pR1;
    double rShort = ((iOffset * cos( pAngle )) - (jOffset * sin( pAngle ))) / pR2;

    return exp( -pow( rLong, 2 ) - pow( rShort, 2 ) );

} // gaussian2D

//
//  calculateKernelSize()
//
//  CJS: 07/08/2015
//
//  Calculate some values such as support, kernel size and w-cell size. We need to work out the maximum baseline length first. The support is given
//  by the square root of the number of uv cells that fit into the maximum baseline, multiplied by 1.5. The w-cell size is given by the maximum
//  baseline length divided by the number of w-planes, multiplied by two.
//

void calculateKernelSize()
{

    // calculate the size of each uv pixel.
    _uvCellSize = (1 / (_uvPixels * (_cellSize / 3600) * (PI / 180)));

    // set the properties of the anti-aliasing kernel.
    _aaSupport = 3;
    _aaKernelSize = (2 * _aaSupport) + 1;

    // set the properties of the w-kernel.
    if (_wPlanes > 1) _wSupport = 30;
    _wKernelSize = (2 * _wSupport) + 1;

    // calculate the support, which should be either the w-support or the aa-support, whichever is larger.
    if (_aaSupport > _wSupport)
        _support = _aaSupport;
    else
        _support = _wSupport;

    // calculate the kernel size.
    _kernelSize = (2 * _support) + 1;

    printf( "_support = %i\n", _support );
    printf( "_aaSupport = %i\n", _aaSupport );
    printf( "_wSupport = %i\n", _wSupport );
    printf( "_cellSize = %f arcsec, %1.12f rad\n", _cellSize, (_cellSize / 3600) * (PI / 180) );
    printf( "_uvCellSize = %f\n", _uvCellSize );

} // calculateKernelSize

//
//  runTask()
//
//  PMW: 16/06/2016
//
//  Run a task using default configuration
//

void runTask(const char *name, task_function fun, std::vector<task_param> params)
{

    // Populate task instance. Resource amounts are made-up for the moment.
    task_instance inst;
    memset(&inst, 0, sizeof(inst));
    inst.name = name;
    inst.resources[TASK_RESOURCE_CORE].amount = 1;
    inst.resources[TASK_RESOURCE_MEMORY].amount = 1000000; // ?

    // Call function
    int ret = fun(TASK_STAGE_ALL, inst, params.size(), &params[0]);

    // Free parameter regions
    for (size_t i = 0; i < params.size(); i++) {
        task_param_free(&params[i]);
    }

    // Check return result
    if (ret) {
        fprintf(stderr, "Task '%s' failed with code %d.\n", name, ret);
        exit(1);
    }

}


//
//  doPhaseCorrection()
//
//  CJS: 12/08/2015
//
//  Phase correct all the visibilities using the PhaseCorrection class.
//

void doPhaseCorrection()
{

    runTask("phase correction", task_phase_correction, {
        make_task_value("in_longitude", _inputRA),
        make_task_value("in_latitude", _inputDEC),
        make_task_value("out_longitude", _outputRA),
        make_task_value("out_latitude", _outputDEC),

        make_task_input("wavelength", _wavelength, "channel", _numChannels),

        make_task_input_output("baselines", &_sample->u,
                               "vis", _numSamples, "bl_coordinate", 3),
        make_task_input_output("visibilities", &_visibility[0],
                               "vis", _numSamples, "channel", _numChannels),
    });

} // doPhaseCorrection

//
//  getTime()
//
//  CJS: 16/11/2015
//
//  Get the elapsed time.
//

float getTime( struct timespec start, struct timespec end )
{

    return ((float)(end.tv_sec - start.tv_sec) * 1000.0) + ((float)(end.tv_nsec - start.tv_nsec) / 1000000.0);

} // getTime

//
//  generateKernel()
//
//  CJS: 16/12/2015
//
//  Generates the convolution function. Currently just a single pixel at the centre of the kernel.
//

bool generateKernel()
{

    // reserve some memory for the kernel.
    _kernel = (complex<double> *) malloc( _kernelSize * _kernelSize *
                                          _oversample * _oversample *
                                          _wPlanes * sizeof( complex<double> ) );

    // Generate w-planes. Unused so far.
    double *w_planes = (double *) malloc( sizeof(double) * _wPlanes );
    memset(w_planes, 0, sizeof(double) * _wPlanes);

    runTask("generate w kernel", task_generate_w_kernel, {
        make_task_value("far_field", int64_t(0)), // unused
        make_task_value("theta", double(0)), // unused

        make_task_input("w_planes", w_planes,
                        "w_plane", _wPlanes),

        make_task_output("kernels", _kernel,
                         "w_plane", _wPlanes,
                         "over_v", _oversample, "over_u", _oversample,
                         "v", _kernelSize, "u", _kernelSize),
    });

    delete [] w_planes;

    return true;

} // generateKernel

//
//  getParameters()
//
//  CJS: 07/08/2015
//
//  Load the following parameters from the parameter file gridder-params: uv cell size, uv grid size, # w-planes, oversample.
//

void getParameters()
{

    char params[80], line[1024], par[80];

    // Open the parameter file and get all lines.
    FILE *fr = fopen( "gridder-params", "rt" );
    while (fgets( line, 80, fr ) != NULL)
    {

        sscanf( line, "%s %s", par, params );
        if (strcmp( par, CELL_SIZE ) == 0)
            _cellSize = atof( params );
        else if (strcmp( par, PIXELS_UV ) == 0)
            _uvPixels = atoi( params );
        else if (strcmp( par, W_PLANES ) == 0)
            _wPlanes = atoi( params );
        else if (strcmp( par, OVERSAMPLE ) == 0)
            _oversample = atof( params );
        else if (strcmp( par, MINOR_CYCLES ) == 0)
            _minorCycles = atoi( params );
        else if (strcmp( par, LOOP_GAIN ) == 0)
            _loopGain = atof( params );
        else if (strcmp( par, INPUT_RA ) == 0)
            _inputRA = atof( params );
        else if (strcmp( par, INPUT_DEC ) == 0)
            _inputDEC = atof( params );
        else if (strcmp( par, OUTPUT_RA ) == 0)
            _outputRA = atof( params );
        else if (strcmp( par, OUTPUT_DEC ) == 0)
            _outputDEC = atof( params );

    }
    fclose( fr );

} // getParameters

//
//  gridVisibilities()
//
//  CJS: 10/08/2015
//
//  Produce an image of gridding visibilities by convolving the complex visibilities with the kernel function.
//

bool gridVisibilities( complex<double> * pGrid, complex<double> * pVisibility,
                       int pOversample, int pKernelSize, int pSupport,
                       complex<double> * pKernel, int pWPlanes, VectorF * pSample,
                       int pNumSamples, int pNumChannels )
{

    // clear the memory.
    memset( pGrid, 0, _uvPixels * _uvPixels * sizeof( complex<double> ) );

    // Generate w-planes. Unused so far.
    double *w_planes = (double *) malloc( sizeof(double) * _wPlanes );
    memset(w_planes, 0, sizeof(double) * _wPlanes);

    // Determine grid bounds
    double uv_min = -_uvPixels * _uvCellSize / 2;
    double uv_max = +_uvPixels * _uvCellSize / 2;

    // Make task call
    runTask("grid visibilities", task_grid_visibilities, {
        make_task_input("wavelength", _wavelength, "channel", _numChannels),
        make_task_input("baselines", &_sample->u,
                        "vis", _numSamples, "bl_coordinate", 3),
        make_task_input("visibilities", &_visibility[0],
                        "vis", _numSamples, "channel", _numChannels),
        make_task_input("w_planes", w_planes,
                        "w_plane", _wPlanes),
        make_task_input_output_bounds("grid", pGrid,
                                      "v", _uvPixels, uv_min, uv_max,
                                      "u", _uvPixels, uv_min, uv_max),
        make_task_input("kernels", _kernel,
                        "w_plane", _wPlanes,
                        "over_v", _oversample, "over_u", _oversample,
                        "v", _kernelSize, "u", _kernelSize),
    });

    delete [] w_planes;

    // return success/failure.
    return true;

} // gridVisibilities

//
//  performFFT()
//
//  CJS: 11/08/2015
//
//  Make a dirty image by inverse FFTing the gridded visibilites.
//

void performFFT( complex<double> * pGrid, int pSize, bool pInverse, bool pConvertToAbsolute, bool pFFTShift )
{

    printf( "performing fft.....\n" );

    // Make task call
    runTask("grid visibilities", task_fourier_transform, {
        make_task_input_output("grid", pGrid,
                               "v", pSize,
                               "u", pSize),
    });

} // performFFT

//
//  getDynamicRange()
//
//  CJS: 10/12/2015
//
//  Gets the dynamic range from the whole of the image.
//

void getDynamicRange()
{

    //const int REGION_SIZE = 50;

    // calculate the average pixel value.
    double total = 0;
    complex<double> * tmp = _grid;
    for ( int i = 0; i < _uvPixels * _uvPixels; i++ )
    {
        total = total + real( *tmp );
        tmp = tmp + 1;
    }
    double averageValue = total / (double)(_uvPixels * _uvPixels);

    // calculate the standard deviation.
    total = 0;
    tmp = _grid;
    for ( int i = 0; i < _uvPixels * _uvPixels; i++ )
    {
        total = total + pow( real( *tmp ) - averageValue, 2 );
        tmp = tmp + 1;
    }
    double standardDeviation = sqrt( total / (double)(_uvPixels * _uvPixels) );

    // calculate the peak value.
    double maxPixel = -1;
    int maxX = 0, maxY = 0;
    for ( int i = 0; i < _uvPixels; i++ )
        for ( int j = 0; j < _uvPixels; j++ )
            if ( real( _grid[ (j * _uvPixels) + i ] ) > maxPixel )
            {
                maxPixel = real( _grid[ (j * _uvPixels) + i ] );
                maxX = i;
                maxY = j;
            }
    printf( "position (%i, %i), value %f, rms %f, dynamic range %f\n", maxX, maxY, maxPixel, standardDeviation, maxPixel / standardDeviation );

} // getDynamicRange

//
//  hogbomClean()
//
//  CJS: 03/11/2015
//
//  Perform a Hogbom clean on our dirty image.
//

void hogbomClean()
{

    // create a clean image
    _cleanImage = (complex<double> *) malloc( _uvPixels * _uvPixels * sizeof( complex<double> ) );

    runTask("hogbom clean", task_hogbom_clean, {
        make_task_value("minor_cycles", int64_t(_minorCycles)),
        make_task_value("loop_gain", double(_loopGain)),
        make_task_value("psf_x", int64_t(_psfX)),
        make_task_value("psf_y", int64_t(_psfY)),

        make_task_input("dirty_beam", _psf, "m", _uvPixels, "l", _uvPixels),
        make_task_input("dirty_image", _grid, "m", _uvPixels, "l", _uvPixels),
        make_task_input("clean_beam", _cleanBeam, "m", _uvPixels, "l", _uvPixels),
        make_task_output("clean_image", _cleanImage, "m", _uvPixels, "l", _uvPixels),

    });

} // hogbomClean

//
//  loadData()
//
//  CJS: 07/08/2015
//
//  Load the data from the measurement set. We need to load a list of
//  sample (the uvw coordinates), a list of channels (we need the
//  frequency of each channel), and a list of visibility (we should
//  have one visibility for each sample/channel combination).
//

bool loadData( char * pInputUVWFilename, char * pInputVisFilename, char * pInputChannelFilename )
{

    bool ok = true;
    char input_a[50], input_b[50], input_c[50];
    char line[1024];

    // read lines from input channel file to count the number of channels.
    _numChannels = 0;
    FILE * channels = fopen( pInputChannelFilename, "rt" );
    while ( fgets( line, 1024, channels ) != NULL )
        _numChannels = _numChannels + 1;
    fclose( channels );

    // reserve memory for wavelengths.
    _wavelength = (double *) malloc( sizeof( double ) * _numChannels );

    // read lines from input channels file again, storing wavelengths.
    int channel = 0;
    channels = fopen( pInputChannelFilename, "rt" );
    while ( fgets( line, 1024, channels ) != NULL )
    {
        sscanf(line, "%s", input_a);
        _wavelength[ channel ] = (float)(CONST_C / atof( input_a ));
        channel = channel + 1;
    }
    fclose( channels );

    // read lines from input file to count the number of samples.
    _numSamples = 0;
    FILE * input = fopen( pInputUVWFilename, "rt" );
    while ( fgets(line, 1024, input) != NULL )
        _numSamples = _numSamples + 1;
    fclose(input);

    // reserve memory for visibilities.
    _sample = (VectorF *) malloc( sizeof( VectorF ) * _numSamples );

    // read lines from input file again, storing visibilies.
    int sample = 0;
    input = fopen( pInputUVWFilename, "rt" );
    while ( fgets(line, 1024, input) != NULL )
    {
        sscanf(line, "%s %s %s", input_a, input_b, input_c);
        _sample[sample].u = atof( input_a );
        _sample[sample].v = atof( input_b );
        _sample[sample].w = atof( input_c );
        sample = sample + 1;
    }
    fclose(input);

    // reserve memory for visibilities.
    _visibility = (complex<double> *) malloc( sizeof( complex<double> ) * _numSamples * _numChannels );

    // read lines from input file again, storing visibilies.
    int visibility = 0;
    input = fopen( pInputVisFilename, "rt" );
    while ( fgets(line, 1024, input) != NULL )
    {
        sscanf(line, "%s %s", input_a, input_b);
        _visibility[visibility] = complex<double>( atof( input_a ), atof( input_b ) );
        visibility = visibility + 1;
    }
    fclose(input);

    // return success/failure.
    return ok;

} // loadData

//
//  saveBitmap()
//
//  CJS: 10/08/2015
//
//  Save the current grid as a bitmap file.
//

bool saveBitmap( char * pOutputFilename, complex<double> * pGrid )
{

    unsigned char * image = NULL;

    const int HEADER_SIZE = 1078;

    // allocate and build the header.
    unsigned char * fileHeader = (unsigned char *) malloc( HEADER_SIZE );
    memset( fileHeader, 0, HEADER_SIZE );

    // file header.
    fileHeader[BIT_CONST] = 'B'; fileHeader[MAP_CONST] = 'M';                   // bfType
    int size = (_uvPixels * _uvPixels) + HEADER_SIZE; memcpy( &fileHeader[IMAGE_SIZE], &size, 4 );  // bfSize
    int offBits = HEADER_SIZE; memcpy( &fileHeader[FILE_HEADER_SIZE], &offBits, 4 );        // bfOffBits

    // image header.
    size = 40; memcpy( &fileHeader[BITMAP_INFO_HEADER], &size, 4 );                 // biSize
    memcpy( &fileHeader[IMAGE_WIDTH], &_uvPixels, 4 );                      // biWidth
    memcpy( &fileHeader[IMAGE_HEIGHT], &_uvPixels, 4 );                     // biHeight
    short planes = 1; memcpy( &fileHeader[COLOUR_PLANES], &planes, 2 );             // biPlanes
    short bitCount = 8; memcpy( &fileHeader[BIT_COUNT], &bitCount, 2 );             // biBitCount
    int coloursUsed = 256; memcpy( &fileHeader[COLOURS_USED], &coloursUsed, 4 );            // biClrUsed

    // colour table.
    for (unsigned int i = 0; i < 256; ++i)
    {
        unsigned int colour = (i << 16) + (i << 8) + i;
        memcpy( &fileHeader[54 + (i * 4)], &colour, 4 );
    }

    bool ok = true;

    // open file.
    FILE * outputFile = fopen( pOutputFilename, "w" );
    if (outputFile == NULL)
    {
        printf( "Could not open file \"%s\".\n", pOutputFilename );
        ok = false;
    }
    else
    {

        // write the file header.
        size_t num_written = fwrite( fileHeader, 1, 1078, outputFile );
        if (num_written != 1078)
        {
            printf( "Error: cannot write to file.\n" );
            ok = false;
        }

        // find the maximum and minimum pixel values.
        double min = abs( pGrid[0] );
        double max = abs( pGrid[0] );
        for ( int i = 0; i < _uvPixels * _uvPixels; i++ )
        {
            if (abs( pGrid[i] ) < min)
                min = abs( pGrid[i] );
            if (abs( pGrid[i] ) > max)
                max = abs( pGrid[i] );
        }

        printf("min - %f, max - %f\n", min, max );

        // add 1% allowance to max - we don't want saturation.
        max = ((max - min) * 1.01) + min;

        // construct the image.
        image = (unsigned char *) malloc( _uvPixels * _uvPixels * sizeof( unsigned char ) );
        for ( int i = 0; i < _uvPixels * _uvPixels; i++ )
            image[i] = (unsigned char)( (abs( pGrid[i] ) - min) * ((double)256 / max) );

        // write the data.
        if (ok == true)
        {

            size_t num_written = fwrite( image, 1, _uvPixels * _uvPixels, outputFile );
            if (num_written != size_t(_uvPixels * _uvPixels))
            {
                printf( "Error: cannot write to file.\n" );
                ok = false;
            }

        }

        // close file.
        fclose( outputFile );

    }

    // cleanup memory.
    free( (void *) fileHeader );
    if (image != NULL)
        free( image );

    // return success flag.
    return ok;

} // saveBitmap

//
//  generateDirtyBeam()
//
//  CJS: 04/11/2015
//
//  Generate the dirty beam by constructing the uv coverage, and then FFT'ing.
//

bool generateDirtyBeam( char * pDirtyBeamFilename )
{

    bool ok = true;
    struct timespec time1, time2;

    clock_gettime( CLOCK_REALTIME, &time1 );

    // create the psf
    _psf = (complex<double> *) malloc( _uvPixels * _uvPixels * sizeof( complex<double> ) );

    // create a new list of visibilities where all values are 1.
    complex<double> * tmpVisibility = (complex<double> *) malloc( sizeof( complex<double> ) * _numSamples * _numChannels );
    for ( int i = 0; i < _numSamples * _numChannels; i++ )
        tmpVisibility[ i ] = 1;

    // generate the uv coverage using the gridder. use a kernel of size 1, no oversampling, and turn off w-projection.
    printf( "\ngridding visibilities for psf.....\n" );
    ok = gridVisibilities( _psf, tmpVisibility, _oversample, _kernelSize, _support, _kernel, _wPlanes, _sample, _numSamples, _numChannels );

    clock_gettime( CLOCK_REALTIME, &time2 );
    float time = getTime( time1, time2 );
    double flops = 8.0 * _numSamples * _numChannels * _kernelSize * _kernelSize;
    printf ("kernel size = %d\n", _kernelSize);
    printf ("visibilities = %d\n", _numSamples);
    fprintf( stderr, "--- time (gridding for dirty beam): (%f ms, %f GOP/s) ---\n\n",
             time, float(flops) / time / 1e6);

    // FFT the uv coverage to get the psf.
    performFFT( _psf, _uvPixels, true, false, true );

    clock_gettime( CLOCK_REALTIME, &time1 );
    fprintf( stderr, "\n--- time (fft): (%f ms) ---\n", getTime( time2, time1 ) );

    // clean up memory.
    free( tmpVisibility );

    // get maximum pixel value.
    double max = -1;
    for ( int i = 0; i < _uvPixels * _uvPixels; i++ )
        if (abs( _psf[ i ] ) > max)
            max = abs( _psf[ i ] );

    // normalise the psf so that the maximum value is 1.
    for ( int i = 0; i < _uvPixels * _uvPixels; i++ )
    {
        _psf[ i ].real( real( _psf[ i ] ) / max );
        _psf[ i ].imag( imag( _psf[ i ] ) / max );
    }

    // save the psf.
    printf( "\npixel intensity range of dirty beam:\n" );
    ok = saveBitmap( pDirtyBeamFilename, _psf );

    clock_gettime( CLOCK_REALTIME, &time2 );
    fprintf( stderr, "\n--- time (save dirty beam): (%f ms) ---\n", getTime( time1, time2 ) );

    // return success flag.
    return ok;

} // generateDirtyBeam

//
//  generateCleanBeam()
//
//  CJS: 04/11/2015
//
//  Generate the clean beam by fitting elliptical Gaussian to the dirty beam.
//

bool generateCleanBeam( char * pCleanBeamFilename )
{

    printf( "\nconstructing clean beam.....\n" );

    struct timespec time1, time2;
    clock_gettime( CLOCK_REALTIME, &time1 );

    // create the clean beam
    _cleanBeam = (complex<double> *) malloc( _uvPixels * _uvPixels * sizeof( complex<double> ) );

    // make task call
    runTask("generate clean beam", task_generate_clean_beam, {
        make_task_value("max_size_for_psf_fitting", int64_t(MAX_SIZE_FOR_PSF_FITTING)),
        make_task_input("dirty_beam", _psf, "m", _uvPixels, "l", _uvPixels),

        make_task_output("clean_beam", _cleanBeam, "m", _uvPixels, "l", _uvPixels),
        make_task_output("psf_x", &_psfX),
        make_task_output("psf_y", &_psfY),
    });

    clock_gettime( CLOCK_REALTIME, &time2 );
    fprintf( stderr, "\n--- time (generate clean beam): (%f ms) ---\n", getTime( time1, time2 ) );

    // save the psf.
    printf( "\npixel intensity range of clean beam:\n" );
    bool ok = saveBitmap( pCleanBeamFilename, _cleanBeam );

    clock_gettime( CLOCK_REALTIME, &time1 );
    fprintf( stderr, "\n--- time (save clean beam): (%f ms) ---\n", getTime( time2, time1 ) );

    // return success flag.
    return ok;

} // generateCleanBeam

//
//  main()
//
//  CJS: 07/08/2015
//
//  Main processing.
//

int main( int pArgc, char ** pArgv )
{

    char DIRTY_BEAM_EXTENSION[] = "-dirty-beam.bmp";
    char CLEAN_BEAM_EXTENSION[] = "-clean-beam.bmp";
    char GRIDDED_EXTENSION[] = "-gridded.bmp";
    char DIRTY_IMAGE_EXTENSION[] = "-dirty-image.bmp";
    char CLEAN_IMAGE_EXTENSION[] = "-clean-image.bmp";
    char RESIDUAL_IMAGE_EXTENSION[] = "-residual-image.bmp";

    char outputDirtyBeamFilename[ 100 ];
    char outputCleanBeamFilename[ 100 ];
    char outputGriddedFilename[ 100 ];
    char outputDirtyImageFilename[ 100 ];
    char outputCleanImageFilename[ 100 ];
    char outputResidualImageFilename[ 100 ];

    struct timespec time1, time2;

    // read program arguments. we expect to see the program call (0), the input filename (1) and the output filename (2).
    if (pArgc != 5)
    {
        printf("Wrong number of arguments. Require the input uvw filename, the input visibility filename, the channel filename, and bitmap prefixes.\n");
        return 1;
    }

    char * inputUVWFilename = pArgv[1];
    char * inputVisFilename = pArgv[2];
    char * inputChannelFilename = pArgv[3];
    char * filenamePrefix = pArgv[4];

    // build filenames.
    strcpy( outputDirtyBeamFilename, filenamePrefix ); strcat( outputDirtyBeamFilename, DIRTY_BEAM_EXTENSION );
    strcpy( outputCleanBeamFilename, filenamePrefix ); strcat( outputCleanBeamFilename, CLEAN_BEAM_EXTENSION );
    strcpy( outputGriddedFilename, filenamePrefix ); strcat( outputGriddedFilename, GRIDDED_EXTENSION );
    strcpy( outputDirtyImageFilename, filenamePrefix ); strcat( outputDirtyImageFilename, DIRTY_IMAGE_EXTENSION );
    strcpy( outputCleanImageFilename, filenamePrefix ); strcat( outputCleanImageFilename, CLEAN_IMAGE_EXTENSION );
    strcpy( outputResidualImageFilename, filenamePrefix ); strcat( outputResidualImageFilename, RESIDUAL_IMAGE_EXTENSION );

    // get the parameters from file.
    getParameters();

    clock_gettime( CLOCK_REALTIME, &time1 );

    // load data. we need to load 'samples' (the uvw coordinates), 'channels' (frequency of each channel) and 'visibilities' (one for each sample/channel combination).
    bool ok = loadData( inputUVWFilename, inputVisFilename, inputChannelFilename );
    if (ok == true)
    {

        clock_gettime( CLOCK_REALTIME, &time2 );
        fprintf( stderr, "\n--- time (load data): (%f ms) ---\n\n", getTime( time1, time2 ) );

        // do phase correction.
        doPhaseCorrection();

        clock_gettime( CLOCK_REALTIME, &time1 );
        fprintf( stderr, "\n--- time (phase rotation): (%f ms) ---\n", getTime( time2, time1 ) );

        // calculate kernel size and other parameters.
        calculateKernelSize();

        clock_gettime( CLOCK_REALTIME, &time2 );
        fprintf( stderr, "\n--- time (calculate kernel size): (%f ms) ---\n\n", getTime( time1, time2 ) );

        // generate kernel.
        ok = generateKernel();

        clock_gettime( CLOCK_REALTIME, &time1 );
        fprintf( stderr, "--- time (generate kernel): (%f ms) ---\n", getTime( time2, time1 ) );

    }

    if (ok == true)
    {

        // generate the dirty beam.
        ok = generateDirtyBeam( outputDirtyBeamFilename );

        clock_gettime( CLOCK_REALTIME, &time2 );
        fprintf( stderr, "\n--- time (generate dirty beam): (%f ms) ---\n", getTime( time1, time2 ) );

    }

    if (ok == true)
    {

        // generate the clean beam.
        ok = generateCleanBeam( outputCleanBeamFilename );

        clock_gettime( CLOCK_REALTIME, &time1 );
        fprintf( stderr, "\n--- time (generate clean beam): (%f ms) ---\n", getTime( time2, time1 ) );

    }

    if (ok == true)
    {

        // create the grid.
        _grid = (complex<double> *) malloc( _uvPixels * _uvPixels * sizeof( complex<double> ) );

        // grid visibilities.
        printf( "\ngridding visibilities for dirty image.....\n" );
        ok = gridVisibilities( _grid, _visibility, _oversample, _kernelSize, _support, _kernel, _wPlanes, _sample, _numSamples, _numChannels );

        clock_gettime( CLOCK_REALTIME, &time2 );
        fprintf( stderr, "--- time (gridding for dirty image): (%f ms) ---\n\n", getTime( time1, time2 ) );

        // normalise the image following the FFT by dividing by the number of pixels.
        for ( int i = 0; i < _uvPixels * _uvPixels; i++ )
            _grid[ i ] = complex<double>(   real( _grid[ i ] ) / (float)(_uvPixels * _uvPixels),
                            imag( _grid[ i ] ) / (float)(_uvPixels * _uvPixels) );

    }

    if (ok == true)
    {

        // save image.
        printf( "pixel intensity range of uv data:\n" );
        ok = saveBitmap( outputGriddedFilename, _grid );

        clock_gettime( CLOCK_REALTIME, &time1 );
        fprintf( stderr, "\n--- time (save gridded image): (%f ms) ---\n\n", getTime( time2, time1 ) );

        // make dirty image.
        performFFT( _grid, _uvPixels, true, false, true );

        clock_gettime( CLOCK_REALTIME, &time2 );
        fprintf( stderr, "\n--- time (fft): (%f ms) ---\n", getTime( time1, time2 ) );

        // save dirty image.
        printf( "\npixel intensity range of dirty image:\n" );
        ok = saveBitmap( outputDirtyImageFilename, _grid );

        clock_gettime( CLOCK_REALTIME, &time1 );
        fprintf( stderr, "\n--- time (save dirty image): (%f ms) ---\n", getTime( time2, time1 ) );

        // get dynamic range.
        getDynamicRange();

        // do a Hogbom clean.
        hogbomClean();

        clock_gettime( CLOCK_REALTIME, &time2 );
        fprintf( stderr, "\n--- time (hogbom clean): (%f ms) ---\n", getTime( time1, time2 ) );

        // save the clean image.
        printf( "\npixel intensity range of clean image:\n" );
        ok = saveBitmap( outputCleanImageFilename, _cleanImage );

        clock_gettime( CLOCK_REALTIME, &time1 );
        fprintf( stderr, "\n--- time (save clean image): (%f ms) ---\n", getTime( time2, time1 ) );

        // save the residual image.
        printf( "\npixel intensity range of residual image:\n" );
        ok = saveBitmap( outputResidualImageFilename, _grid );

        clock_gettime( CLOCK_REALTIME, &time2 );
        fprintf( stderr, "\n--- time (save residual image): (%f ms) ---\n", getTime( time1, time2 ) );

    }

    // clear memory.
    if (_sample != NULL)
        free( _sample );
    if (_wavelength != NULL)
        free( _wavelength );
    if (_visibility != NULL)
        free( _visibility );
    if (_kernel != NULL)
        free( _kernel );
    if (_grid != NULL)
        free( _grid );
    if (_psf != NULL)
        free( _psf );
    if (_cleanBeam != NULL)
        free( _cleanBeam );
    if (_cleanImage != NULL)
        free( _cleanImage );

    return 0;

} // main
