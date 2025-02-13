all: compile run clean

compile:
	g++ -std=c++20 root/src/main.cpp root/src/atom.cpp root/src/cluster.cpp root/src/lennard_jones.cpp root/src/derivative_approximation.cpp -o main -Iroot/include -I/opt/homebrew/include -L/opt/homebrew/lib -larmadillo -w

run:
	./main || true

clean:
	rm main