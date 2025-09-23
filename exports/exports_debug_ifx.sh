FPM_FC="ifx"
export FPM_FC

FPM_FFLAGS_DEBUG="-DDEBUG -O0 -g -check bounds -check uninit -check pointer -traceback"
FPM_FFLAGS="$FPM_FFLAGS_DEBUG -DUSING_IFX"
export FPM_FFLAGS

OMP_NUM_THREADS=8
export OMP_NUM_THREADS
