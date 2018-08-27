cd ~
rm SLiM.zip

#download SLiM
wget http://benhaller.com/slim/SLiM.zip
unzip SLiM.zip

#make build directory
mkdir slim_build
cd slim_build

#load newer version of gcc
module load gcc/4.9.3

#tell cmake which one to use and compile
CC=gcc CXX=g++ cmake ../SLiM
make slim