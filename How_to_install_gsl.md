##How to install gsl:       Credit [Gss-installation](https://coral.ise.lehigh.edu/jild13/2016/07/11/hello/)



step1: wget ftp://ftp.gnu.org/gnu/gsl/gsl-2.6.tar.gz                   *Help: this will download gsl-2.6  or you can download other newer available version of gsl*
step2: tar -zxvf ~/Downloads/gsl-2.6.tar.gz                            *Help: Extract the gz file*
step3  cd gsl-2.6/                                                     *Help: change the directory to gsl-2.6*
step4  ./configure --prefix=/mnt/home/student/camit/.local/lib/python3.7/site-packages/gsl_installationp
                                                                       *Help:  ./configure command checks our system for the required software needed to build the program,prefix to direct where you cwant to intstall gsl*
step5   make                                                           *Help: make is typically used to build executable programs and libraries from source code*
step6   make check                                                     *Help: make check is a command to a makefile*
step7   make install                                                   *Help: "make install", the make program takes the binaries from the previous step and copies them into some appropriate locations so that they can be accessed.*
step8   vim example.c                                                  *Help: Create a simple c file or Copy the content of example.c(already present in this repository)*
step9   gcc -Wall -I/mnt/home/student/camit/.local/lib/python3.7/site-packages/gsl_installation/include -c ~/example.c
                                                                       *Help: compiling example.c*
step10  gcc -L/mnt/home/student/camit/.local/lib/python3.7/site-packages/gsl_installation/lib  example.o -lgsl -lgslcblas -lm
                                                                       *Help: compiling example.o* 
step11  LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/mnt/home/student/camit/.local/lib/python3.7/site-packages/gsl_installation/lib
                                                                       *Help: LD_LIBRARY_PATH is a colon-separated set of directories where libraries should be searched*
step12  export LD_LIBRARY_PATH                                         *Help: export LD_LIBRARY_PATH*
###Note###: Commands corresponding to Step 11 & 12 sets the path for libraries to be search for at first.It is exteremely recommended to put these commads(step 11&12) in your ###.bashrc### file if you dont wish to manually type them again and again.  
 
 
step3  ./a.out                                                         *Help: Will work without error, if everything went  fine*
