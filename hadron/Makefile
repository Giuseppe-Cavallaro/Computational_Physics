MAIN = main
EXEC = $(MAIN)

$(MAIN): $(MAIN).cpp
	g++ $(MAIN).cpp -o $(EXEC) -I`./.lhapdf-config --incdir` -L`./.lhapdf-config --libdir` -static -lLHAPDF

run: $(MAIN)
	./$(EXEC)
