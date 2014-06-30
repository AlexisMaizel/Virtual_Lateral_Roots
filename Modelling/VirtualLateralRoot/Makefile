OBJECTS=model.o bezier.o surface.o ModelExporter.o ModelUtils.o

default: model.vve

standalone: model

include $(shell vveinterpreter --makefile)
# Uncomment the next line to compile with debug options
#CXXFLAGS=$(CXXFLAGS_DEBUG)
# Add extra compilation options here
CXXFLAGS+=
# Add extra libraries here
LIBS+=
# Add extra linking options here
#    for compiling the model as a library (.vve)
LD_SO_FLAGS+=
#    for compiling the model as a stand-alone program
LD_EXE_FLAGS+=

model.o: model.moc

bezier.o: bezier.cpp bezier.h

surface.o: surface.cpp surface.h

model.vve: $(OBJECTS)

model${EXESUFFIX}: $(OBJECTS)

run:
	vveinterpreter model.vve

clean:
	$(RM) *.o model model.vve

