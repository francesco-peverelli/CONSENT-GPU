#!/bin/bash

cd BMAN-GPU
./install.sh;
cd ..;
cd minimap2;
make;
cd ..;
make;
