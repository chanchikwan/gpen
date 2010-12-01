CUDA = /usr/local/cuda
NVCC = $(CUDA)/bin/nvcc
LINK = $(addprefix -Xlinker , -rpath $(CUDA)/lib)
FLGS = $(addprefix --compiler-options , -Wall)

double:
	$(NVCC) *.cu -O3 $(LINK) $(FLGS) -o run.x -DOUBLE -arch sm_13

float:
	$(NVCC) *.cu -O3 $(LINK) $(FLGS) -o run.x

clean:
	-rm run.x
	-rm *~
