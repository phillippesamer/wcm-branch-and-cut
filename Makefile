.RECIPEPREFIX = >

UNAME_S := $(shell uname -n)
ifeq ($(UNAME_S),djervheta)
    GRB_PATH    = /opt/gurobi1002/linux64
    GRB_INCLUDE = -I$(GRB_PATH)/include/
    GRB_LINK    = -L$(GRB_PATH)/lib/ -lgurobi_g++5.2 -lgurobi100
else
    ifeq ($(UNAME_S),ii3102747.dyn.uib.no)
        GRB_PATH    = /scratch/gurobi1002/linux64
        LEMON_PATH  = /scratch/lemon-1.3.1/lemon
        GRB_INCLUDE = -I$(GRB_PATH)/include/ -I$(LEMON_PATH)
        GRB_LINK    = -L$(GRB_PATH)/lib/ -lgurobi_g++5.2 -lgurobi100
    else
        GRB_PATH    = /opt/gurobi1002/linux64
        GRB_INCLUDE = -I$(GRB_PATH)/include/ 
        GRB_LINK    = -L$(GRB_PATH)/lib/ -lgurobi_g++5.2 -lgurobi100
    endif
endif

CC          = g++ -Wall -Wextra -O3 -m64

FILES_CC    = graph.cpp io.cpp wcm_model.cpp wcm_cutgenerator.cpp main.cpp

BINARY      = wcm

all: clean compile

clean:
> find . -name '*.o' -exec rm -f '{}' ';'
> rm -f $(BINARY);

compile:
> $(CC)  -o $(BINARY)  $(FILES_CC)  $(GRB_INCLUDE)  $(GRB_LINK)  -lm -lemon
