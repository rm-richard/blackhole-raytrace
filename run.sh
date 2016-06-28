#!/bin/bash
_FILENAME=main.cpp
g++ $_FILENAME -w -g -W -Wall -o output -lGL -lglut -lm

if [[ $? == 0 ]]; then
  ./output
fi
