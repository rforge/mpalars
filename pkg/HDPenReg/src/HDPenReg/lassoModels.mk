## lars makefile for R package

include $(R_HOME)/etc${R_ARCH}/Makeconf

all:lib

#-----------------------------------------------------------------------
# Variables
# 
#LIB = ./libHD.a
LIB = ./lib/libHDPenReg.a

#SRC_DIR = ..
STK_INC_DIR = -I./  -I../

#-----------------------------------------------------------------------
# Sources files
#
SRCS =  ./lassoModels/LassoPenalty.cpp \
	./lassoModels/EnetPenalty.cpp \


#-------------------------------------------------------------------------
# generate the variable OBJS containing the names of the object files
#
OBJS= $(SRCS:./lassoModels/%.cpp=%.o)

#-------------------------------------------------------------------------
# rule for compiling the cpp files
#
%.o: ./lassoModels/%.cpp
	$(CXX) $(CXXFLAGS) ${CPICFLAGS} $(STK_INC_DIR) $< -c -o $@

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
