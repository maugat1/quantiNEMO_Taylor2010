CC=g++ -O3 

SOURCES=$(shell ls *.cpp)
OBJECTS=$(shell for file in $(SOURCES);\
                do echo -n $$file | sed -e "s/\(.*\)\.cpp/\1\.o/"; echo -n " ";\
                done)

PRGNAME=quantiNemo
INCS	=	-I. #-D_DEBUG

all:  bin

bin : objects
	echo $(OBJECTS)
	$(CC) *.o -o $(PRGNAME) $(LIBS)

objects : $(OBJECTS)

%.o : %.cpp
#	removeLF.sh
	$(CC) -c -w $< -o $@ $(INCS)

clean:
	rm -f *.o $(PRGNAME) 

depend:
#	removeLF.sh
	$(CC) -M *.cpp > $@ $(COMDIR)






#include depend
