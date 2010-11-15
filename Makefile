CUDA = /usr/local/cuda
NVCC = $(CUDA)/bin/nvcc
LINK = $(addprefix -Xlinker , -rpath $(CUDA)/lib)
FLGS = $(addprefix --compiler-options , -Wall)

compile:
	$(NVCC) *.cu -O3 $(LINK) $(FLGS) -o run.x

clean:
	-rm run.x
	-rm *~
