# ray tracing makefile to build every .c, .cpp or .cxx file in the directory

# list of all .o files that could be generated from .c, .cpp or .cxx 
# uses all source files in this directory
# typically, you would list these out by hand to say exactly what to use
OBJS = $(patsubst %.c,   %.o, $(wildcard *.c)) \
       $(patsubst %.cpp, %.o, $(wildcard *.cpp)) \
       $(patsubst %.cxx, %.o, $(wildcard *.cxx))

# set C and C++ compiler flags to give all warnings
# warnings are telling you something is likely wrong. 
# if you get any, FIX YOUR CODE to get rid of the warning
CFLAGS = -Wall
CXXFLAGS = -Wall

# optimization flags (switch to -g for debugging)
OPT = -O4

# first rule says what to build if you say "make" instead of "make balls2.png"
balls2.png:

# make a .png file from a .ppm file with the same name
# $@ is the current target (the png), $< is the (first) thing it depends on
%.png: %.ppm
	convert $< $@

# how to make a ppm from program and nff file
# $@ will be the ppm, $* is the base name without the .nff or .ppm
# make as trace.ppm, then rename
%.ppm: %.nff trace
	./trace < $*.nff
	mv trace.ppm $@

# how to build trace executable from OBJS files
# link with c++ compiler
# $@ is the current target (trace)
# set OBJ variable to -g or -O for debugging or optimization
# set LDFLAGS to any library directorys (-L...) to search
# set LDLIBS to any libraries to use (-l...)
# should be OK to leave all of those blank
trace: $(OBJS)
	$(CXX) $(OPT) -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)

# how to build a .o from a .c if the .c is newer
# $@ is the current target, $< is the .o it depends on
# you can set OPT as above to -g or -O if you want
# set CFLAGS to any additional C compiler options
%.o: %.c
	$(CC) $(OPT) -c -o $@ $(CFLAGS) $<

# how to build a .o from a .cxx
# set CXXFLAGS to any additional C++ compiler options
%.o: %.cxx
	$(CXX) $(OPT) -c -o $@ $(CXXFLAGS) $<

# how to build a .o from a .cpp
%.o: %.cpp
	$(CXX) $(OPT) -c -o $@ $(CXXFLAGS) $<

# target to remove all of the object files (say "make clean")
clean:
	rm -f *.o

# target to do the stuff in clean + remove trace and any ppm or png images
clobber: clean
	rm -f trace *.ppm *.png

# I recommend adding source and header dependencies by running
#   g++ -MM *.c, g++ -MM *.cpp, or g++ -MM *.cxx
# and pasting the output here
