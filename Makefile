CXX		= g++
CXXFLAGS	= -O3 -g -Wall -std=c++0x -I ../eigen -I ../boost_1_55_0
LDFLAGS		= -static
LIBS		= -lm -lgsl

BIN		= synth
OBJS	= main.o Tract.o SpeechSynthesizer.o

.SUFFIXES: .cxx .o

all:		$(BIN)

$(BIN):		$(OBJS)
			$(CXX) $(OBJS) $(LDFLAGS) $(LIBS) -o $(BIN)

main.cxx:	Tract.h SpeechSynthesizer.h
Tract.cxx:	Tract.h Const.h
SpeechSynthesizer.cxx:	SpeechSynthesizer.h Const.h

.cxx.o:
			$(CXX) $(CXXFLAGS) -c $<

clean:
			rm *.o $(BIN)
