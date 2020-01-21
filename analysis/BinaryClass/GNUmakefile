testStruct: testStruct.o genStruct.o genField.o
	g++ -o testStruct testStruct.o genStruct.o genField.o
testStruct.o: testStruct.cc genStruct.hh
	g++ -c -o $@ $<
genStruct.o: genStruct.cc genStruct.hh genField.hh
	g++ -c -o $@ $<
genField.o: genField.cc genField.hh
	g++ -c -o $@ $<

.PHONY: clean
clean:
	rm -f testStruct testStruct.o genStruct.o genField.o test_Struct.dat

