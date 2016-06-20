
TASK_SOURCES = \
	phase_correction.cpp phasecorrection.cpp generate_kernel.cpp grid_visibilities.cpp \
	fourier_transform.cpp hogbom_clean.cpp generate_clean_beam.cpp

MAIN_SOURCES = \
	cppGridderNoKernel.cpp

TASK_OBJS = $(TASK_SOURCES:%.cpp=task/%.o)
MAIN_OBJS = $(MAIN_SOURCES:%.cpp=%.o)

CFLAGS +=  -Wall -ggdb -O0
CXXFLAGS += -Wall -I.. -Itask -ggdb -O0 -std=c++11
LDFLAGS += -lfftw3 -lfftw3_omp -fopenmp -ggdb -O0

gridder : $(OSKAR_OBJS) $(MAIN_OBJS) $(TASK_OBJS)
	$(CXX) $(LDFLAGS) $^ -o $@

.PHONY: clean
clean :
	rm -f $(OSKAR_OBJS) $(MAIN_OBJS) $(TASK_OBJS)
