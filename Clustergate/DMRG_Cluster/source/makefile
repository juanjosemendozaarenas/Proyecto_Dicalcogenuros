include makefile_mkl.inc

.PHONY: all clean Ground_State

all: Ground_State

Ground_State: ../bin/Ground_State

../bin/%: ../obj/%.o  
	mkdir --parents ../bin
	${FF} -o $@ $< $(FFLAGS) 
	
../obj/%.o: %.c $(DEPS)
	mkdir --parents ../obj
	$(CC) -c -o $@ $< $(CFLAGS) 

clean:
	rm ../obj/* ../bin/*
