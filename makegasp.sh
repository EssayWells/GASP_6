#!/bin/bash
COMPILER="g++ -O2"
echo "Compiling binaries from cpp files"
$COMPILER -c basic_utils.cpp Vectors.cpp Rotors.cpp Geometry.cpp Structure.cpp GASPcontrol.cpp GASPmain.cpp
$COMPILER -c duper.cpp  super.cpp coordinates2xtl.cpp xtl2cif.cpp polycomp.cpp shifter.cpp rotordotcom.cpp
$COMPILER -c xtlduper.cpp xtlsuper.cpp cif2xtl.cpp xtlshifter.cpp netdr2.cpp polyana.cpp transformCell.cpp
echo "Compiling gasp executable"
$COMPILER -o gasp basic_utils.o Vectors.o Rotors.o Geometry.o Structure.o GASPcontrol.o GASPmain.o
echo "Compiling duper executable"
$COMPILER -o duper basic_utils.o Vectors.o Rotors.o Geometry.o Structure.o duper.o 
$COMPILER -o xtlduper basic_utils.o Vectors.o Rotors.o Geometry.o Structure.o xtlduper.o 
echo "Compiling super executable"
$COMPILER -o super basic_utils.o Vectors.o Rotors.o Geometry.o Structure.o super.o 
$COMPILER -o xtlsuper basic_utils.o Vectors.o Rotors.o Geometry.o Structure.o xtlsuper.o 
echo "Compiling coordinates2xtl executable"
$COMPILER -o coordinates2xtl basic_utils.o Vectors.o Rotors.o Geometry.o Structure.o coordinates2xtl.o 
echo "Compiling xtl2cif and cif2xtl executable"
$COMPILER -o xtl2cif basic_utils.o Vectors.o Rotors.o Geometry.o Structure.o xtl2cif.o
$COMPILER -o cif2xtl basic_utils.o Vectors.o Rotors.o Geometry.o Structure.o cif2xtl.o
echo "Compiling polyhedron comparison polycomp and polyhedron analysis polyana executables."
$COMPILER -o polycomp basic_utils.o Vectors.o Rotors.o polycomp.o 
$COMPILER -o polyana basic_utils.o Vectors.o Rotors.o polyana.o 
echo "Compiling shifter utility"
$COMPILER -o shifter basic_utils.o Vectors.o Rotors.o Structure.o shifter.o
$COMPILER -o xtlshifter basic_utils.o Vectors.o Rotors.o Structure.o xtlshifter.o
echo "Compiling rotordotcom comparison utility"
$COMPILER -o rotordotcom basic_utils.o Vectors.o rotordotcom.o
echo "Compiling netdr2 comparison utility"
$COMPILER -o netdr2 basic_utils.o Vectors.o Rotors.o Geometry.o Structure.o netdr2.o
echo "Compiling transformCell redraw utility"
$COMPILER -o transformCell basic_utils.o Vectors.o Rotors.o Geometry.o Structure.o transformCell.o


