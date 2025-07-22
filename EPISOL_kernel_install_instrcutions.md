Navigate to your directory of choice and create a folder to do installation
```
cd $HOME/
mkdir episol
cd episol
```
Download and extract Fourier transform libraries
```
wget https://github.com/EPISOLrelease/EPISOL/raw/refs/heads/main/src/fftw/fftw-3.3.8.tar.gz
tar -xzf fftw-3.3.8.tar.gz
cd fftw-3.3.8/
#!./configure --prefix=/home/fftw-3.3.8
```
Configure and make libraries
```
./configure --prefix=$HOME/episol/fftw-3.3.8
make
make install
```
return to your build directory and download and extract the episol kernel
```
cd ../
wget https://github.com/EPISOLrelease/EPISOL/raw/refs/heads/main/src/kernel/release.tar.gz
tar -xzf release.tar.gz
cd release/
```
Configure using the path of your installed FFTW libraries and install :)
```
./configure --with-fftw=$HOME/episol/fftw-3.3.8
make
make install
```
