
mkdir -p build

## Compile each module of RNMC
echo "building GMC"
cd ../GMC
make clean
make GMC

# echo "building RNMC"
# cd ../NPMC
# make clean
# make NPMC

# echo "building LGMC"
# cd ../LGMC
# make clean
# make LGMC