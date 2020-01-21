COPTS = $(shell root-config --cflags)

LIBS = $(shell root-config --libs)

ConvertToBinary: ConvertToBinary.o genStruct.o genField.o
	g++ -o ConvertToBinary ConvertToBinary.o genStruct.o genField.o $(LIBS)
ConvertToBinary.o: ConvertToBinary.cc genStruct.hh
	g++ $(COPTS) -c -o $@ $<
genStruct.o: genStruct.cc genStruct.hh genField.hh
	g++ $(COPTS) -c -o $@ $<
genField.o: genField.cc genField.hh
	g++ $(COPTS) -c -o $@ $<

.PHONY: clean
clean:
	rm -f ConvertToBinary ConvertToBinary.o genStruct.o genField.o

