CC = mpicxx
HEADER = evp.h evpFlags.h V3math.h 
CFLAGS = -Wall -O3
SEEKFLAG = -DMPICH_IGNORE_CXX_SEEK
#INCLUDES = -I/public/home/pzhao/Soft_mu01/include 
#LFLAGS = -L/public/home/pzhao/Soft_mu01/lib
MPI_CFLAGS = -I/path/mpich/3.2.1/include
MPI_LIBS = -L/path/mpich/3.2.1/lib
FFTW_CFLAGS = -I/path/fftw/2.5.1/include
FFTW_LIBS = -L/path/fftw/2.5.1/lib
LIBS = -lfftw_mpi -lfftw -lgsl -lgslcblas -lm
SRCS = evp.c io.c init.c kinematics.c evolution.c constitutive.c   testf.c
OBJS = $(SRCS:.c=.o)
MAIN = evp
.PHONY: depend clean
#all: $(MAIN)

$(MAIN): $(OBJS) Makefile
	$(CC) $(SEEKFLAG) $(MPI_CFLAGS) $(FFTW_CFLAGS) $(CFLAGS)  -o $(MAIN) $(OBJS) $(MPI_LIBS) $(FFTW_LIBS)  $(LIBS)

%.o: %.c $(HEADER) Makefile
	$(CC) $(SEEKFLAG) $(MPI_CFLAGS) $(FFTW_CFLAGS) -c $< -o $@

clean:
	$(RM) *.o *~ $(MAIN)

depend: $(SRCS)
	makedepend $^
