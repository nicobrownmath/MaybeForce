CC = g++
CFLAGS = -std=c++17 -Wall -Wextra -pedantic -O3 -march=native -mtune=native -fno-stack-protector -fomit-frame-pointer #-flto
LDFLAGS =

all: MaybeForce
MaybeForce: MaybeForce.cpp LifeAPI.h
	$(CC) $(CFLAGS) $(INSTRUMENTFLAGS) -o MaybeForce MaybeForce.cpp $(LDFLAGS)