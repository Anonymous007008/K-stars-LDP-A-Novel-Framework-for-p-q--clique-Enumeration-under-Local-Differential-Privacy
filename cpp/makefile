CXXFLAGS	= -march=native -O3
CXXFLAGS	+= -std=c++11
CXXFLAGS	+= -isystem ./include
CXXFLAGS	+= -MMD
CC	= g++


k-stars_LDP: k-stars_LDP.o

clean:
	rm -f *~ $(TARGET) *.d *.o

-include *.d
