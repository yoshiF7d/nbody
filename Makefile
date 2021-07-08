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
SRCS = $(wildcard $(SRCDIR)/*.cpp)
OBJS = $(addprefix $(OBJDIR)/,$(notdir $(SRCS:.cpp=.o)))

$(TARGET): $(OBJS)
	$(COMPILER) $(CFLAGS) -o ./bin/nbody $^ $(LDFLAGS)

mpi:$(OBJS)
	$(COMPILER) $(CFLAGS) -o ./bin/mpinbody $^ $(LDFLAGS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(COMPILER) $(CFLAGS) $(INCLUDE) -o $@ -c $<

.PHONY : clean cleanall
clean :
	rm -f $(OBJS) $(TARGET)
