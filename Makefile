TARGET = nbody
SRCDIR = ./src
OBJDIR = ./obj
INCDIR = ./inc

COMPILER = clang++
mpi:COMPILER = mpic++
CFLAGS = -std=c++11 
mpi:CFLAGS = -std=c++11 -DMPI
INCLUDE = -I $(INCDIR)
LDFLAGS = 
SRCS = $(SRCDIR)/nbody.cpp $(SRCDIR)/Particle.cpp
OBJS = $(addprefix $(OBJDIR)/,$(notdir $(SRCS:.cpp=.o)))

$(TARGET): $(OBJS)
	$(COMPILER) $(CFLAGS) -o ./bin/nbody $^ $(LDFLAGS)

nbody0 : $(SRCDIR)/nbody0.cpp $(SRCDIR)/Particle.cpp $(SRCDIR)/ParticleList.cpp
	$(COMPILER) $(CFLAGS) -o ./bin/nbody0 $^ $(INCLUDE) $(LDFLAGS)

mpi:$(OBJS)
	$(COMPILER) $(CFLAGS) -o ./bin/mpinbody $^ $(LDFLAGS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(COMPILER) $(CFLAGS) $(INCLUDE) -o $@ -c $<

.PHONY : clean cleanall
clean :
	rm -f $(OBJS) $(TARGET)
