# Directory for each of the libraries required
DIR_MAT := /lustre/apps/MATLAB/R2014b

DIR_ARPACK := /hpcfs/home/ciencias/fisica/postgrado/fj.gomez34/ARPACK
DIR_TNT := /hpcfs/home/ciencias/fisica/postgrado/fj.gomez34/.lib/tnt_v1.2.1
DIR_INT := /hpcfs/home/ciencias/fisica/postgrado/fj.gomez34/.pkgs/intel/parallel_studio_xe_2017.7.065/compilers_and_libraries_2017/linux

DIR_MKL = ${DIR_INT}/mkl/lib/intel64
DIR_MTH = ${DIR_INT}/compiler/lib/intel64

# Define the linked MATLAB libraries and MATLAB flags
MATLIBS := -L$(DIR_MAT)/bin/glnxa64 -L$(TMW_ROOT)/sys/os/glnxa64 -lmat -lmx
RPATH := -Wl,-rpath-link,$(DIR_MAT)/bin/glnxa64,-rpath-link,$(DIR_MAT)/sys/os/glnxa64
MATFLAGS := -I$(DIR_MAT)/extern/include 

# Define the linked TNT libraries and TNT flags
TNTLIBS := $(DIR_TNT)/lib/icc_apk_mat/libtnt.a $(DIR_TNT)/lib/icc_apk_mat/libtntMps.a
TNTFLAGS := -I$(DIR_TNT)/include  

# Define the NAG flags
MKLFLAGS := -DLONG_CHARLEN -DRETURN_COMPLEX_PARAM

# Define all linked libraries
ALLLIBS := $(MATLIBS) -Wl,--start-group $(DIR_MKL)/libmkl_intel_lp64.a $(DIR_MKL)/libmkl_intel_thread.a
ALLLIBS := $(ALLLIBS) ${DIR_MKL}/libmkl_core.a ${DIR_MTH}/libiomp5.a ${DIR_ARPACK}/libarpack.a $(TNTLIBS) -Wl,--end-group -lpthread -nofor-main
ALLLIBS := $(ALLLIBS) $(RPATH)

# Define all library flags
LIBFLAGS := $(MATFLAGS) $(TNTFLAGS) $(MKLFLAGS)

# Specify the compilers to use and any additional compiler flags
CC := icc
FF := ifort

CFLAGS := $(LIBFLAGS)
FFLAGS := $(ALLLIBS)

# Extra flags
#CFLAGS += -DTNT_OMP -openmp -DTNTARPACKAVAIL -DTNT_SAFESVD
#FFLAGS += -DTNT_OMP -openmp -DTNTARPACKAVAIL -DTNT_SAFESVD

CFLAGS += -DTNT_OMP -qopenmp -DTNTARPACKAVAIL -DTNT_SAFESVD
FFLAGS += -DTNT_OMP -qopenmp -DTNTARPACKAVAIL -DTNT_SAFESVD

MAKE := make
MKFLAGS =
