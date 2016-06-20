
#ifndef C_ABI_H
#define C_ABI_H

#include <stdint.h>

#ifdef __cplusplus
#include <complex>
typedef std::complex<float> float_complex;
typedef std::complex<double> double_complex;
#else
#include <complex.h>
typedef double complex float_complex;
typedef double complex double_complex;
#endif


enum task_resource
{
    /// Number of nodes allocated to a task
    TASK_RESOURCE_NODE,

    /// Number of cores that are allocated to a task
    TASK_RESOURCE_CORE,
    /// Amount of main memory allocated to a task
    TASK_RESOURCE_MEMORY,

    /// Number of GPUs that are allocated to a task
    // (can be shared as required?)
    TASK_RESOURCE_GPU,
    /// Amount of GPU memory allocated to a task
    TASK_RESOURCE_GPU_MEMORY,

    TASK_RESOURCE_COUNT
};

struct task_resource_desc
{
    /// Amount of this resource allocated to the task
    int amount;
    /// Pointer to resource-specific data that further describes the
    // allocated data.
    void *handle;
};

/// A task instance.
struct task_instance
{
    /// Implemented task
    const char *name;
    /// User data pointer, for retaining data between calls
    void *user;
    /// Amount of resources allocated. (-1) means don't care.
    task_resource_desc resources[TASK_RESOURCE_COUNT];
};

enum task_param_kind
{
    /// Just a plain value. Expected to have no regions attached to it.
    TASK_PARAM_KIND_VALUE,
    /// Input buffer. Zero or more regions attached, might have data
    /// set on a START call.
    TASK_PARAM_KIND_INPUT,
    /// Output buffer. Zero or more regions attached
    TASK_PARAM_KIND_OUTPUT,
    /// Shared output buffer. This means that the buffer in question
    // already has data, which is expected to get updated by the task.
    TASK_PARAM_KIND_INPUT_OUTPUT,
};

enum task_param_type
{

    /// 32-bit signed little endian integer numbers
    TASK_PARAM_TYPE_INT_32,
    /// 64-bit signed little endian integer numbers
    TASK_PARAM_TYPE_INT_64,

    /// 32-bit IEEE 754 single precision floating point numbers
    TASK_PARAM_TYPE_FLOAT,
    /// 64-bit IEEE 754 double precision floating point numbers
    TASK_PARAM_TYPE_DOUBLE,

    /// Pairs of 32-bit IEEE 754 single precision floating point numbers
    TASK_PARAM_TYPE_COMPLEX_FLOAT,
    /// Pairs of 64-bit IEEE 754 double precision floating point numbers
    TASK_PARAM_TYPE_COMPLEX_DOUBLE,

};

/// Task region specification
//
// This declares how the array is layed out in memory as well as what
// the semantics of the array coordinates are. The idea here is that
// if we are encoding an underlying function "f", we have the
// following data mapping:
//
//   param.data[x * x_stride + y * y_stride] ==
//     encode(f(x_min + x * (x_max - x_min) / x_extent,
//              y_min + y * (y_max - y_min) / y_extent))
//
// This equation should hold for all integer x in [0..x_extent-1] and
// all integer y in [0..y_extent-1]. This is C-style array access, so
// the concrete pointer will be multiplied by the element-size of the
// data array, which depends on the underlying data type.
struct task_region
{
    /// Domain name
    const char *name;
    /// Extent of the array coordinates
    int64_t extent;
    /// Stride of the array coordinates
    int64_t stride;
    /// Minimum coordinate encoded
    double min;
    /// Maximum coordinate encoded
    double max;
};

struct task_param
{
    /// Parameter name
    const char *name;
    /// type of buffer
    task_param_kind kind;
    /// Parameter element data type
    task_param_type typ;
    /// Pointer to regions specifying the parameter layout (see
    /// documentation of task_region). Can be null for
    /// zero-dimensional constants.
    task_region *regions;
    /// Length of regions array. Must be 0 if regions is null.
    int region_count;
    /// Data.
    union {
        void *data;
        int32_t i32;
        int64_t i64;
        float f32;
        double f64;
    };
};

/// Task stage bitmask. When multiple bits are given, the strongest
// guarantees hold for the passed parameters.

/// Verify stage: Only check that enough resources are made available
// and that data layouts and kinds are acceptable.
//
// * Data is not expected to be valid.
// * Data layout can still change (downwards?)
// * User data will not be carried over
const int TASK_STAGE_VERIFY  = (1<<0);

/// Prepare for execution. This call will be made before any real data
// is passed to the kernel in order to allow for the kernel to prepare
// planning data. User data can be set and will be passed to
// subsequent calls.
//
// * Data is not expected to be valid.
// * Data layout can still change (downwards?)
const int TASK_STAGE_PREPARE = (1<<1);

/// Start the task. Data will be supplied with the exact data layout
// provided.
//
// * Data might still be incomplete (?)
// * Output buffers might already be set and can be updated in whatever
//   way the task sees fit.
const int TASK_STAGE_START   = (1<<2);

/// Finish processing. After this call has returned, the task must
// have updated all output buffers to completion.
//
// Note that additional START calls might happen after the last FINISH
// call, corresponding to the kernel getting used multiple times.
const int TASK_STAGE_FINISH  = (1<<3);

/// Free user data, perform cleanups
//
// Called when the task is to be freed.
const int TASK_STAGE_FREE    = (1<<4);


/// Trivial task interface: All stages at once. Calling a task with
// this bitmask means that the kernel should perform all verification
// and preparation, perform the work, write out the output and clean
// up before returning control to the execution engine.
const int TASK_STAGE_ALL =
    TASK_STAGE_VERIFY | TASK_STAGE_PREPARE |
    TASK_STAGE_START | TASK_STAGE_FINISH |
    TASK_STAGE_FREE;

/// The task function signature
typedef int (*task_function)(int stage,
                             task_instance inst,
                             int argc,
                             const task_param* argv
                             );

#endif // C_ABI_H
