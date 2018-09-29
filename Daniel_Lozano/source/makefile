include makefile.inc

OBJMAIN := $(patsubst %,$(OBJDIR)/%,$(OBJMAIN))

.PHONY: all clean $(SUBDIRS) $(MKFILE)

$(MKFILE): all

all: ${OBJMAIN} $(SUBDIRS) $(EXDIR)/Ground_State_Ladder

$(EXDIR)/Ground_State_Ladder: $(OBJ)
	${FF} -o $@ $(OBJ) $(NAGFLAGS) $(CLIBS) $(FFLAGS)

$(SUBDIRS):
	mkdir --parents $(OBJDIR)/$@
	cd $@; $(MAKE) $(MKFLAGS) $(MKFILE)

$(OBJDIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

clean:
	-rm $(OBJ) $(EXDIR)/Ground_State_Ladder
