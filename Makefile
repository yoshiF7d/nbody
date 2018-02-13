TARGET = ./bin/nbody
SRCDIR = ./src
OBJDIR = ./obj
INCDIR = ./inc

COMPILER = clang++
CFLAGS = 
INCLUDE = -I $(INCDIR)
LDFLAGS = 
SRCS = $(wildcard $(SRCDIR)/*.cpp)
OBJS = $(addprefix $(OBJDIR)/,$(notdir $(SRCS:.cpp=.o)))

$(TARGET): $(OBJS)
	echo $(SRCS)
	$(COMPILER) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(COMPILER) $(CFLAGS) $(INCLUDE) -o $@ -c $<

.PHONY : clean cleanall
clean :
	rm -f $(OBJS) $(TARGET)
