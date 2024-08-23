
mkdir -p build

## Compile each module of RNMC
pwd
echo "building GMC"
cd GMC
make clean
make GMC
cd ..
pwd

echo "building NPMC"
cd NPMC
make clean
make NPMC
cd ..

echo "building LGMC"
cd LGMC
make clean
make LGMC
cd ..