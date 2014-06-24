## lars makefile for R package

ifeq ($(WIN),64)
  include $(R_HOME)/etc${R_ARCH}/Makeconf
endif
all:lib

#-----------------------------------------------------------------------
# Variables
# 
LIB = ./lib/libHDPenReg.a

STK_INC_DIR = -I./  -I../

#-----------------------------------------------------------------------
# Sources files
#
SRCS =  ./lassoModels/LassoPenalty.cpp \
	./lassoModels/LassoSolver.cpp \
	./lassoModels/FusedLassoPenalty.cpp \
	./lassoModels/FusedLassoSolver.cpp \
	./lassoModels/LogisticLassoPenalty.cpp \
	./lassoModels/LogisticLassoSolver.cpp \
	./lassoModels/LogisticFusedLassoPenalty.cpp \
	./lassoModels/LogisticFusedLassoSolver.cpp \
	./lassoModels/CV.cpp \
#	./lassoModels/EnetPenalty.cpp \


#-------------------------------------------------------------------------
# generate the variable OBJS containing the names of the object files
#
OBJS= $(SRCS:./lassoModels/%.cpp=%.o)

#-------------------------------------------------------------------------
# rule for compiling the cpp files
#
%.o: ./lassoModels/%.cpp
	$(CXX) $(CXXFLAGS)  ${CPICFLAGS} $(STK_INC_DIR) $< -c -o $@

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