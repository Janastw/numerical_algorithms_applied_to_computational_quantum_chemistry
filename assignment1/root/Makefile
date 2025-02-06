all: compile run clean

compile:
	g++ -std=c++20 src/main.cpp src/atom.cpp src/cluster.cpp src/lennard_jones.cpp src/derivative_approximation.cpp -o main -Iinclude -I/opt/homebrew/include -L/opt/homebrew/lib -larmadillo -w

run:
	./main || true

clean:
	rm main