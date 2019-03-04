DIR_SRC = ./src
DIR_OBJ = ./obj
BINDIR = ./bin

SRC = $(wildcard ${DIR_SRC}/*.cpp)  
OBJ = $(patsubst %.cpp,${DIR_OBJ}/%.o,$(notdir ${SRC})) 

TARGET = fuspipe

BIN_TARGET = ${TARGET}

CC = g++
CFLAGS = -std=c++11 -g

${BIN_TARGET}:${OBJ}
	$(CC) $(OBJ) -L. -lpthread -o $@
    
${DIR_OBJ}/%.o:${DIR_SRC}/%.cpp make_obj_dir
	$(CC) $(CFLAGS) -O3 -c $< -o $@
.PHONY:clean
clean:
	rm -rf $(DIR_OBJ)
	rm $(TARGET)
	rm -rf $(BINDIR)

make_obj_dir:
	@if test ! -d $(DIR_OBJ) ; \
	then \
		mkdir $(DIR_OBJ) ; \
	fi

install:
	@if test ! -d $(BINDIR); \
	    then \
	    mkdir $(BINDIR); \
	fi
	install $(TARGET) $(BINDIR)/$(TARGET)
	@echo "Installed."
