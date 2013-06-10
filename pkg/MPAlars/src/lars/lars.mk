## lars makefile for R package

include $(R_HOME)/etc${R_ARCH}/Makeconf

all:lib

#-----------------------------------------------------------------------
# Variables
# 
LIB = ../libMPA.a

#SRC_DIR = ..
STK_INC_DIR = -I../

#-----------------------------------------------------------------------
# Sources files
#
SRCS =  Lars.cpp \
	Path.cpp \
	PathState.cpp \
	functions.cpp \
	Cvlars.cpp \



#-------------------------------------------------------------------------
# generate the variable OBJS containing the names of the object files
#
OBJS= $(SRCS:%.cpp=%.o)

#-------------------------------------------------------------------------
# rule for compiling the cpp files
#
%.o: %.cpp
	$(CXX) $(CXXFLAGS) ${CPICFLAGS} ${STK_INC_DIR} $< -c -o $@

#-----------------------------------------------------------------------
# The rule lib create the library 
#
lib: $(LIB)

$(LIB): $(OBJS)
	$(AR) -r $@ $?
  
mostlyclean: clean

clean:
	@-rm -rf .libs _libs $(LIB)
	@-rm -f *.o
