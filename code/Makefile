##########################  Macro Definitions  ############################
# Let's try to auto-detect what platform we're on.  
# If this fails, set PLATFORM manually in the else block.
AUTOPLATFORM = Failed
ifeq ($(MSYSTEM),MINGW32)
	AUTOPLATFORM = Succeeded
	PLATFORM = __MSYS__
endif
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
	AUTOPLATFORM = Succeeded
	PLATFORM = __linux__
endif
ifeq ($(UNAME_S),Darwin)
	AUTOPLATFORM = Succeeded
	PLATFORM = __APPLE__
endif
ifeq ($(AUTOPLATFORM),Failed)
	# Autodetect failed.  Set platform manually.
	#PLATFORM = __APPLE__
	#PLATFORM = __linux__
	PLATFORM = __MSYS__
endif

NAME = SOCRATES-Retrieval

# Basic directories
PROJDIR = ../
OUT = ../
CODE = $(PROJDIR)code/
OBJ = $(CODE)object/
INC = $(CODE)include/
SRC = $(CODE)source/
42KIT = $(SRC)42/

ifeq ($(PLATFORM),__APPLE__)
	# Mac Macros
	CINC = -I /usr/include
	LIBS = -framework System -lnetcdf
	LFLAGS = -bind_at_load
	# ARCHFLAG = -arch i386
	# ARCHFLAG = -arch x86_64
	ARCHFLAG = -arch arm64
	EXENAME = $(NAME)
	CC = mpicc
endif

ifeq ($(PLATFORM),__linux__)
	# Linux Macros
	CINC =
	LIBS = -ldl -lm -lnetcdf
	LFLAGS =
	ARCHFLAG =
	EXENAME = $(NAME)
	CC = mpicc
endif

ifeq ($(PLATFORM),__MSYS__)
	CINC =
	LIBS =  -lws2_32 -lnetcdf
	LFLAGS =
	ARCHFLAG =
	EXENAME = $(NAME).exe
	CC = mpicc
endif

SOOPOBJ = $(OBJ)socratesret.o $(OBJ)init.o $(OBJ)geometry.o $(OBJ)ground.o $(OBJ)product.o $(OBJ)antenna.o $(OBJ)util.o $(OBJ)42kit.o

RETOBJ = $(OBJ)retrieval.o

CFLAGS = -Wall -Wshadow -Wno-deprecated -g $(CINC) -I $(INC) -I $(SRC) -O0 $(ARCHFLAG) 

##########################  Rules to link smatret  #############################

SOCRATES-retrieval : $(SOOPOBJ) $(RETOBJ)
	$(CC) $(LFLAGS) -o $(OUT)$(EXENAME) $(SOOPOBJ) $(RETOBJ) $(LIBS)

####################  Rules to compile objects  ###########################

$(OBJ)socratesret.o      : $(SRC)socratesret.c $(INC)socratesret.h
	$(CC) $(CFLAGS) -c $(SRC)socratesret.c -o $(OBJ)socratesret.o

$(OBJ)init.o             : $(SRC)init.c $(INC)socratesret.h
	$(CC) $(CFLAGS) -c $(SRC)init.c -o $(OBJ)init.o

$(OBJ)geometry.o         : $(SRC)geometry.c $(INC)socratesret.h
	$(CC) $(CFLAGS) -c $(SRC)geometry.c -o $(OBJ)geometry.o

$(OBJ)ground.o           : $(SRC)ground.c $(INC)socratesret.h
	$(CC) $(CFLAGS) -c $(SRC)ground.c -o $(OBJ)ground.o

$(OBJ)product.o          : $(SRC)product.c $(INC)socratesret.h
	$(CC) $(CFLAGS) -c $(SRC)product.c -o $(OBJ)product.o

$(OBJ)antenna.o          : $(SRC)antenna.c $(INC)socratesret.h
	$(CC) $(CFLAGS) -c $(SRC)antenna.c -o $(OBJ)antenna.o

$(OBJ)util.o             : $(SRC)util.c $(INC)util.h $(INC)socratesret.h
	$(CC) $(CFLAGS) -c $(SRC)util.c -o $(OBJ)util.o

$(OBJ)retrieval.o        : $(SRC)retrieval.c $(INC)socratesret.h
	$(CC) $(CFLAGS) -c $(SRC)retrieval.c -o $(OBJ)retrieval.o
 
$(OBJ)42kit.o            : $(42KIT)42kit.c $(INC)socratesret.h
	$(CC) $(CFLAGS) -c $(42KIT)42kit.c -o $(OBJ)42kit.o
########################  Miscellaneous Rules  ############################
clean :
ifeq ($(PLATFORM),_WIN32)
	del .\obj\*.o .\$(OUT)$(EXENAME)
else ifeq ($(PLATFORM),_WIN64)
	del .\obj\*.o .\$(OUT)$(EXENAME)
else
	rm -f $(OBJ)*.o .\$(OUT)$(EXENAME)
endif


