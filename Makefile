CC=g++
CFLAGS=-O3 -lgsl -lgslcblas -pthread -lboost_math_c99
TARGET=main
DATA=wavefunction
SCRIPT=pl.py
OUTPUT=output_plot_python.png
OUTPUT1=output1.png
OUTPUT2=output2.png
OUTPUT3=output3.png

all: $(TARGET)

$(TARGET): main.cpp
	$(CC) $^ $(CFLAGS) -o $(TARGET)

run: $(TARGET)
	./$(TARGET) $(n) $(l) $(m)

plot: run
	python $(SCRIPT)

view: plot
	kitten icat $(OUTPUT)

clean:
	rm -f $(TARGET) 

.PHONY: all run plot view clean

