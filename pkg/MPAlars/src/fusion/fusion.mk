## fusion makefile for R package

include $(R_HOME)/etc${R_ARCH}/Makeconf

all:lib

#-----------------------------------------------------------------------
# Variables
# 
LIB = ../libMPA.a
STK_INC_DIR = -I../

#-----------------------------------------------------------------------
# Sources files
#
SRCS =  Fusion.cpp \


#-------------------------------------------------------------------------
# generate the variable OBJS containing the names of the object files
#
OBJS= $(SRCS:%.cpp=%.o)

#-------------------------------------------------------------------------
# rule for compiling the cpp files
#
%.o: %.cpp
	$(CXX) $(CXXFLAGS)    ${CPICFLAGS} ${STK_INC_DIR} $< -c -o $@

#-----------------------------------------------------------------------
# The rule lib create the library 
#
lib: $(LIB)

$(LIB): $(OBJS)
	$(AR) $(ARFLAGS) $@ $?
  
mostlyclean: clean

clean:
	@-rm -rf .libs _libs $(LIB)
	@-rm -f *.o
