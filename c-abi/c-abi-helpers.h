
#ifndef C_ABI_HELPERS_H
#define C_ABI_HELPERS_H

#include "c-abi.h"

#include <assert.h>
#include <string.h>
#include <stdio.h>

/// Helpers for writing tasks
inline const task_param *task_find_param(int argc, const task_param *argv, const char *name,
                                         task_param_kind kind, task_param_type typ, int regions)
{
    for (int i = 0; i < argc; i++) {
        if (!strcmp(name, argv[i].name)) {

            if (argv[i].kind != kind) {
                fprintf(stderr, "kind mismatch on parameter %s!\n", name);
                continue;
            }
            if (argv[i].typ != typ) {
                fprintf(stderr, "type mismatch on parameter %s!\n", name);
                continue;
            }

            if (argv[i].region_count < regions) {
                fprintf(stderr, "region count mismatch on parameter %s: %d < %d!\n", name,
                        argv[i].region_count, regions);
                continue;
            }

            return argv + i;
        }
    }
    fprintf(stderr, "parameter %s not found!\n", name);
    return 0;
}
inline const task_param *task_value(int argc, const task_param *argv,
                                    const char *name, task_param_type typ)
    { return task_find_param(argc, argv, name, TASK_PARAM_KIND_VALUE, typ, 0); }
inline const task_param *task_input(int argc, const task_param *argv,
                                    const char *name, task_param_type typ, int regions)
    { return task_find_param(argc, argv, name, TASK_PARAM_KIND_INPUT, typ, regions); }
inline const task_param *task_output(int argc, const task_param *argv,
                                     const char *name, task_param_type typ, int regions)
    { return task_find_param(argc, argv, name, TASK_PARAM_KIND_OUTPUT, typ, regions); }
inline const task_param *task_input_output(int argc, const task_param *argv,
                                           const char *name, task_param_type typ, int regions)
    { return task_find_param(argc, argv, name, TASK_PARAM_KIND_INPUT_OUTPUT, typ, regions); }

inline double &task_f64(const task_param *param) {
    assert (param->kind != TASK_PARAM_KIND_VALUE);
    double *array = reinterpret_cast<double*>(param->data);
    return array[0];
}

inline double &task_f64(const task_param *param, int x) {
    assert(param->region_count >= 1);
    double *array = reinterpret_cast<double*>(param->data);
    return array[x * param->regions[0].stride];
}

inline double &task_f64(const task_param *param, int x, int y) {
    assert(param->region_count >= 2);
    double *array = reinterpret_cast<double*>(param->data);
    return array[x * param->regions[0].stride +
                 y * param->regions[1].stride];
}

inline double_complex &task_c64(const task_param *param, int x, int y) {
    assert(param->region_count >= 2);
    double_complex *array = reinterpret_cast<double_complex*>(param->data);
    return array[x * param->regions[0].stride +
                 y * param->regions[1].stride];
}

inline double_complex &task_c64(const task_param *param, int x, int y, int z, int a, int b) {
    assert(param->region_count >= 5);
    double_complex *array = reinterpret_cast<double_complex*>(param->data);
    return array[x * param->regions[0].stride +
                 y * param->regions[1].stride +
                 z * param->regions[2].stride +
                 a * param->regions[3].stride +
                 b * param->regions[4].stride];
}

inline task_param make_task_value(const char *name, int64_t val) {
    task_param p = {};
    memset(&p, 0, sizeof(p));
    p.name = name;
    p.kind = TASK_PARAM_KIND_VALUE;
    p.typ = TASK_PARAM_TYPE_INT_64;
    p.i64 = val;
    return p;
}

inline task_param make_task_value(const char *name, double val) {
    task_param p = {};
    memset(&p, 0, sizeof(p));
    p.name = name;
    p.kind = TASK_PARAM_KIND_VALUE;
    p.typ = TASK_PARAM_TYPE_DOUBLE;
    p.f64 = val;
    return p;
}

inline task_param make_task_input(const char *name, double *val,
                                  const char *name_x, int64_t extent_x) {
    task_param p = {};
    memset(&p, 0, sizeof(p));
    p.name = name;
    p.kind = TASK_PARAM_KIND_INPUT;
    p.typ = TASK_PARAM_TYPE_DOUBLE;
    p.data = reinterpret_cast<void *>(val);
    p.region_count = 1;
    p.regions = new task_region[p.region_count];
    memset(p.regions, 0, p.region_count * sizeof(*p.regions));
    p.regions[0].name = name_x;
    p.regions[0].extent = extent_x;
    p.regions[0].stride = 1;
    return p;
}

inline task_param make_task_input(const char *name, double *val,
                                  const char *name_x, int64_t extent_x,
                                  const char *name_y, int64_t extent_y) {
    task_param p = {};
    memset(&p, 0, sizeof(p));
    p.name = name;
    p.kind = TASK_PARAM_KIND_INPUT;
    p.typ = TASK_PARAM_TYPE_DOUBLE;
    p.data = reinterpret_cast<void *>(val);
    p.region_count = 2;
    p.regions = new task_region[p.region_count];
    memset(p.regions, 0, p.region_count * sizeof(*p.regions));
    p.regions[0].name = name_x;
    p.regions[0].extent = extent_x;
    p.regions[0].stride = extent_y;
    p.regions[1].name = name_y;
    p.regions[1].extent = extent_y;
    p.regions[1].stride = 1;
    return p;
}

inline task_param make_task_input(const char *name, double_complex *val,
                                  const char *name_x, int64_t extent_x,
                                  const char *name_y, int64_t extent_y) {
    task_param p = {};
    memset(&p, 0, sizeof(p));
    p.name = name;
    p.kind = TASK_PARAM_KIND_INPUT;
    p.typ = TASK_PARAM_TYPE_COMPLEX_DOUBLE;
    p.data = reinterpret_cast<void *>(val);
    p.region_count = 2;
    p.regions = new task_region[p.region_count];
    memset(p.regions, 0, p.region_count * sizeof(*p.regions));
    p.regions[0].name = name_x;
    p.regions[0].extent = extent_x;
    p.regions[0].stride = extent_y;
    p.regions[1].name = name_y;
    p.regions[1].extent = extent_y;
    p.regions[1].stride = 1;
    return p;
}

inline task_param make_task_input(const char *name, double_complex *val,
                                  const char *name_x, int64_t extent_x,
                                  const char *name_y, int64_t extent_y,
                                  const char *name_z, int64_t extent_z,
                                  const char *name_a, int64_t extent_a,
                                  const char *name_b, int64_t extent_b) {
    task_param p = {};
    memset(&p, 0, sizeof(p));
    p.name = name;
    p.kind = TASK_PARAM_KIND_INPUT;
    p.typ = TASK_PARAM_TYPE_COMPLEX_DOUBLE;
    p.data = reinterpret_cast<void *>(val);
    p.region_count = 5;
    p.regions = new task_region[p.region_count];
    memset(p.regions, 0, p.region_count * sizeof(*p.regions));
    p.regions[0].name = name_x;
    p.regions[0].extent = extent_x;
    p.regions[0].stride = extent_y * extent_z * extent_a * extent_b;
    p.regions[1].name = name_y;
    p.regions[1].extent = extent_y;
    p.regions[1].stride = extent_z * extent_a * extent_b;
    p.regions[2].name = name_z;
    p.regions[2].extent = extent_z;
    p.regions[2].stride = extent_a * extent_b;
    p.regions[3].name = name_a;
    p.regions[3].extent = extent_a;
    p.regions[3].stride = extent_b;
    p.regions[4].name = name_b;
    p.regions[4].extent = extent_b;
    p.regions[4].stride = 1;
    return p;
}

inline task_param make_task_output(const char *name, int64_t *val) {
    task_param p = {};
    memset(&p, 0, sizeof(p));
    p.name = name;
    p.kind = TASK_PARAM_KIND_OUTPUT;
    p.typ = TASK_PARAM_TYPE_INT_64;
    p.data = reinterpret_cast<void *>(val);
    p.region_count = 0;
    p.regions = 0;
    return p;
}

inline task_param make_task_output(const char *name, double_complex *val,
                                   const char *name_x, int64_t extent_x,
                                   const char *name_y, int64_t extent_y) {
    task_param p = {};
    memset(&p, 0, sizeof(p));
    p.name = name;
    p.kind = TASK_PARAM_KIND_OUTPUT;
    p.typ = TASK_PARAM_TYPE_COMPLEX_DOUBLE;
    p.data = reinterpret_cast<void *>(val);
    p.region_count = 2;
    p.regions = new task_region[p.region_count];
    memset(p.regions, 0, p.region_count * sizeof(*p.regions));
    p.regions[0].name = name_x;
    p.regions[0].extent = extent_x;
    p.regions[0].stride = extent_y;
    p.regions[1].name = name_y;
    p.regions[1].extent = extent_y;
    p.regions[1].stride = 1;
    return p;
}

inline task_param make_task_output(const char *name, double_complex *val,
                                   const char *name_x, int64_t extent_x,
                                   const char *name_y, int64_t extent_y,
                                   const char *name_z, int64_t extent_z,
                                   const char *name_a, int64_t extent_a,
                                   const char *name_b, int64_t extent_b) {
    task_param p = {};
    memset(&p, 0, sizeof(p));
    p.name = name;
    p.kind = TASK_PARAM_KIND_OUTPUT;
    p.typ = TASK_PARAM_TYPE_COMPLEX_DOUBLE;
    p.data = reinterpret_cast<void *>(val);
    p.region_count = 5;
    p.regions = new task_region[p.region_count];
    memset(p.regions, 0, p.region_count * sizeof(*p.regions));
    p.regions[0].name = name_x;
    p.regions[0].extent = extent_x;
    p.regions[0].stride = extent_y * extent_z * extent_a * extent_b;
    p.regions[1].name = name_y;
    p.regions[1].extent = extent_y;
    p.regions[1].stride = extent_z * extent_a * extent_b;
    p.regions[2].name = name_z;
    p.regions[2].extent = extent_z;
    p.regions[2].stride = extent_a * extent_b;
    p.regions[3].name = name_a;
    p.regions[3].extent = extent_a;
    p.regions[3].stride = extent_b;
    p.regions[4].name = name_b;
    p.regions[4].extent = extent_b;
    p.regions[4].stride = 1;
    return p;
}

inline task_param make_task_input_output(const char *name, double *val,
                                         const char *name_x, int64_t extent_x,
                                         const char *name_y, int64_t extent_y) {
    task_param p = {};
    memset(&p, 0, sizeof(p));
    p.name = name;
    p.kind = TASK_PARAM_KIND_INPUT_OUTPUT;
    p.typ = TASK_PARAM_TYPE_DOUBLE;
    p.data = reinterpret_cast<void *>(val);
    p.region_count = 2;
    p.regions = new task_region[p.region_count];
    memset(p.regions, 0, p.region_count * sizeof(*p.regions));
    p.regions[0].name = name_x;
    p.regions[0].extent = extent_x;
    p.regions[0].stride = extent_y;
    p.regions[1].name = name_y;
    p.regions[1].extent = extent_y;
    p.regions[1].stride = 1;
    return p;
}

inline task_param make_task_input_output(const char *name, double_complex *val,
                                         const char *name_x, int64_t extent_x,
                                         const char *name_y, int64_t extent_y) {
    task_param p = {};
    memset(&p, 0, sizeof(p));
    p.name = name;
    p.kind = TASK_PARAM_KIND_INPUT_OUTPUT;
    p.typ = TASK_PARAM_TYPE_COMPLEX_DOUBLE;
    p.data = reinterpret_cast<void *>(val);
    p.region_count = 2;
    p.regions = new task_region[p.region_count];
    memset(p.regions, 0, p.region_count * sizeof(*p.regions));
    p.regions[0].name = name_x;
    p.regions[0].extent = extent_x;
    p.regions[0].stride = extent_y;
    p.regions[1].name = name_y;
    p.regions[1].extent = extent_y;
    p.regions[1].stride = 1;
    return p;
}

inline task_param make_task_input_output_bounds(
        const char *name, double_complex *val,
        const char *name_x, int64_t extent_x, double min_x, double max_x,
        const char *name_y, int64_t extent_y, double min_y, double max_y) {

    task_param p = {};
    memset(&p, 0, sizeof(p));
    p.name = name;
    p.kind = TASK_PARAM_KIND_INPUT_OUTPUT;
    p.typ = TASK_PARAM_TYPE_COMPLEX_DOUBLE;
    p.data = reinterpret_cast<void *>(val);
    p.region_count = 2;
    p.regions = new task_region[p.region_count];
    memset(p.regions, 0, p.region_count * sizeof(*p.regions));
    p.regions[0].name = name_x;
    p.regions[0].extent = extent_x;
    p.regions[0].stride = extent_y;
    p.regions[0].min = min_x;
    p.regions[0].max = max_x;
    p.regions[1].name = name_y;
    p.regions[1].extent = extent_y;
    p.regions[1].stride = 1;
    p.regions[1].min = min_y;
    p.regions[1].max = max_y;
    return p;
}

inline void task_param_free(task_param *param)
{
    if (param->kind != TASK_PARAM_KIND_VALUE) {
        delete param->regions;
        param->regions = 0;
    }
}

/// Returns the byte size of an array element of the given parameter
inline size_t task_param_elem_size(const task_param *param)
{
    switch(param->typ) {
    case TASK_PARAM_TYPE_INT_32: return sizeof(int32_t);
    case TASK_PARAM_TYPE_INT_64: return sizeof(int64_t);
    case TASK_PARAM_TYPE_FLOAT: return sizeof(float);
    case TASK_PARAM_TYPE_DOUBLE: return sizeof(double);
    case TASK_PARAM_TYPE_COMPLEX_FLOAT: return sizeof(float_complex);
    case TASK_PARAM_TYPE_COMPLEX_DOUBLE: return sizeof(double_complex);
    }
    // Should not happen
    assert(false);
    return 1;
}

/// Returns the size in bytes that the given parameter array occupies.
inline size_t task_param_size(const task_param *param)
{
    size_t elem_size = task_param_elem_size(param);
    int64_t max_length = 0;
    // Determine maximum length
    for (int64_t i = 0; i < param->region_count; i++) {
        int64_t len = param->regions[i].extent * param->regions[i].stride;
        if (len > max_length) {
            max_length = len;
        }
    }
    // Multiply out and return.
    return elem_size * size_t(max_length);
}

#endif // C_ABI_HELPERS_H
