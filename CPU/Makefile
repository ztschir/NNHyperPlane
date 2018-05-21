.SUFFIXES: .cpp .o

CC = g++

srcdir = .
CFLAGS = -O3
LFLAGS = 

all: NNHyperPlaneSplit

clean:
	@rm -f *.o NNHyperPlaneSplit

##################################################################
# Rules
##################################################################

NNHyperPlaneSplit: NNHyperPlaneSplit.o
	@echo  "Linking" $@ "... "
	${CC} -o NNHyperPlaneSplit NNHyperPlaneSplit.o ${LFLAGS}

##################################################################
# Generic rules
##################################################################

.cpp.o:
	${CC} -o $@ -c ${CFLAGS} $<

