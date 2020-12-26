#!/bin/bash

# Compile the C code into a library.
gcc -o base_aux.so -shared -fPIC -O2 base_aux.c