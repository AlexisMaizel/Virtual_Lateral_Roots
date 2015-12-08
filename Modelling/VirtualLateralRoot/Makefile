OBJECTS=model.o bezier.o surface.o ModelExporter.o SurfacePoints.o RealSurface.o ModelUtils.o DivisionSetting.o GraphicsClass.o InitSurface.o SurfaceBaseClass.o

default: model.vve

standalone: model

include $(shell vveinterpreter --makefile)
# Uncomment the next line to compile with debug options
#CXXFLAGS=$(CXXFLAGS_DEBUG)
# Add extra compilation options here
CXXFLAGS+= -w -std=c++11
# Add extra libraries here
LIBS+= -lftgl
# Add extra linking options here
#    for compiling the model as a library (.vve)
LD_SO_FLAGS+=
#    for compiling the model as a stand-alone program
LD_EXE_FLAGS+=

model.o: model.moc model.h ModelHeader.h PrincipalComponentAnalysis.h

DivisionSetting.o: DivisionSetting.cpp DivisionSetting.h

GraphicsClass.o: GraphicsClass.cpp GraphicsClass.h

SurfaceBaseClass.o: SurfaceBaseClass.cpp SurfaceBaseClass.h

ModelUtils.o: ModelUtils.cpp ModelUtils.h

bezier.o: bezier.cpp bezier.h

InitSurface.o: InitSurface.cpp InitSurface.h

SurfacePoints.o: SurfacePoints.cpp SurfacePoints.h

RealSurface.o: RealSurface.cpp RealSurface.h

surface.o: surface.cpp surface.h

model.vve: $(OBJECTS)

model${EXESUFFIX}: $(OBJECTS)

run:
	vveinterpreter model.vve

clean:
	$(RM) *.o model model.vve


