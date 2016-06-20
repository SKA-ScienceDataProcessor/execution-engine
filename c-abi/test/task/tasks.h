
#ifndef TASKS_H
#define TASKS_H

#include "c-abi.h"

// Domains of the oversampled convolution kernel stack
enum K_DOMAINS {
    K_PLANE,
    K_OVER_V,
    K_OVER_U,
    K_V,
    K_U,

    K_COUNT
};

#ifdef __cplusplus
extern "C" {
#endif

int task_phase_correction(int stage,
                          task_instance inst,
                          int argc,
                          const task_param* argv);

int task_generate_w_kernel(int stage,
                           task_instance inst,
                           int argc,
                           const task_param* argv);

int task_grid_visibilities(int stage,
                           task_instance inst,
                           int argc,
                           const task_param* argv);

int task_fourier_transform(int stage,
                           task_instance inst,
                           int argc,
                           const task_param* argv);

int task_generate_clean_beam(int stage,
                             task_instance inst,
                             int argc,
                             const task_param* argv);

int task_hogbom_clean(int stage,
                      task_instance inst,
                      int argc,
                      const task_param* argv);
#ifdef __cplusplus
}
#endif

#endif
