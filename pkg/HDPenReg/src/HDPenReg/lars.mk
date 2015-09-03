## lars makefile for R package

#include $(R_HOME)/etc${R_ARCH}/Makeconf
all:lib

#-----------------------------------------------------------------------
# Variables
# 
PKGLIB = ./lib/libHDPenReg.a

#-----------------------------------------------------------------------
# Sources files
#
SRCS =  ./lars/Lars.cpp \
  ./lars/Path.cpp \
	./lars/PathState.cpp \
	./lars/functions.cpp \
  ./lars/Fusion.cpp \
	./lars/Cvlars.cpp \



#-------------------------------------------------------------------------
# generate the variable OBJS containing the names of the object files
#
OBJS= $(SRCS:./lars/%.cpp=%.o)

#-------------------------------------------------------------------------
# rule for compiling the cpp files
#
%.o: ./lars/%.cpp
	$(CXX) $(HD_CPPFLAGS1) $(HD_CXXFLAGS)  $< -c -o $@

#-----------------------------------------------------------------------
# The rule lib create the library 
#
lib: $(PKGLIB)

$(PKGLIB): $(OBJS)
	$(AR) -r $@ $?
  
mostlyclean: clean

clean:
	@-rm -rf .libs _libs $(PKGLIB)
	@-rm -f *.o