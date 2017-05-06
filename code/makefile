main: main.cpp Grid.o FluidSimulation.o Particle.o
	g++ -o main Grid.o FluidSimulation.o Particle.o main.cpp -lGL -lglut -lm

FluidSimulation: FluidSimulation.cpp FluidSimulation.hpp Grid.o
	g++ -c FluidSimulation.cpp FluidSimulation.hpp Grid.o -lm

Grid: Grid.cpp Grid.hpp
	g++ -c Grid.cpp Grid.hpp

Particle: Particle.cpp Particle.hpp
		g++ -c Particle.cpp Particle.hpp

clean:
	rm *.o
	rm main
	rm *.gch
