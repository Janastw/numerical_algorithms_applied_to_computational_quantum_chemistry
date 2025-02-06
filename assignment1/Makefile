all: compile run clean

compile:
	g++ -std=c++20 h1.cpp -o h1 -I/opt/homebrew/include -L/opt/homebrew/lib -larmadillo -w

run:
	./h1

clean:
	rm h1