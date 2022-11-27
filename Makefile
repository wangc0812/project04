SRC_DIR = ./src
SOURCE  = $(wildcard $(SRC_DIR)/*.c) # Search all the .c files in the current directory, and return to SOURCE

OBJS    = $(patsubst %.c, %.o, $(SOURCE)) # Replace all .c files with .o files 
TARGET  = main
INCLUDE = -I./inc # -I means search files in the specificed folder
# all .h files are in inc
# all .cpp files are in src

CC = gcc
CFLAGES = -o3 -c -Wall -g
CFLAG = -o3 -g #  -g gdb debuger
OMPFLAGS = -fopenmp
SIMD = -mavx2

# options pass to the compiler
# -c generates the object file
# -Wall displays complier warning
#  $@: object files
#  $^: all the prerequisitites files
#  $<: the first prerequisite file


$(TARGET) : $(OBJS)
	$(CC) $(SIMD) $(OMPFLAGS) $(CFLAG) -o $@ $(OBJS) 

%.o : %.c # all the .o objects depend on the .c files
	$(CC) $(SIMD) $(OMPFLAGS) $(CFLAGES) $< -o $@ $(INCLUDE)


.PHONY : clean # .PHONY will prevent making from confusing the phony target with a file name
clean : 
	rm -f $(SRC_DIR)/*.o $(TARGET)