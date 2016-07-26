# RESOLUTE_AC
Create an attenuation correction map from UTE TE sequences

## Introduction
RESOLUTE is a tool for creating MR attenuation maps from UTE TE images.
The TE images are segmented for brain, CSF air, soft tissue and continuous bone values.
The result is a DCM dataset.

## Usage
The program, after compilation, can simply be called by
<pre><code>
RESOLUTE < path_to_TE1 > < path_to_TE2 > < path_to_umap > < output_folder >
</code></pre>
The program takes about 3 minutes to run.

## Installation (Ubuntu and CentOS)

### Ubuntu 64 Bit
Install below dependencies first. Then:
<pre><code>
git clone --recursive git://github.com/claesnl/RESOLUTE_AC.git RESOLUTE_AC
cd RESOLUTE_AC/linux
mkdir build && cd build
cmake ..
make && sudo make install
</code></pre>

### CentOS
If you wish to use the configure file to install the dependencies (third line), run the following code. During the `./configure` step, you will be prompted for sudo password and ccmake GUI input several times. The total run-time to install is about 20 minutes.

Otherwise, install dependencies first (see below) and skip the third line.
<pre><code>
git clone --recursive git://github.com/claesnl/RESOLUTE_AC.git RESOLUTE_AC
cd RESOLUTE_AC/centos
./configure
mkdir build && cd build
cmake .. -DDCMTK_DIR=/usr/local/include/dcmtk/
make && sudo make install
</code></pre>

## Install from existing RESOLUTE version
<pre><code>
cd RESOLUTE_AC
git clean -f -d
git pull origin master
</code></pre>
and repeat the build step from above. 

## Dependencies

### Ubuntu:

#### Add MNI to source list to install autoreg
<pre><code>
echo 'deb http://packages.bic.mni.mcgill.ca/ubuntu-maverick/ ./' | sudo tee -a /etc/apt/sources.list
sudo apt-get update
</code></pre>

#### Install minc tools, DCMTK, ANTS and other compile dependencies
<pre><code>
sudo apt-get install libnetcdf-dev libhdf5-dev hdf5-tools \
   libdcmtk2-dev libdcmtk2 dcmtk minc-tools ants \
   mni-autoreg mni-models-icbm152-lin bicpl libmni-perllib-perl \
   g++ cmake libx11-dev bicpl libmni-perllib-perl
   </code></pre>

### CentOS:
All of the below dependencies and programs can be installed with `./configure` under the centos folder. Manual install instructions shown below.

#### Dependencies:
<pre><code>
sudo yum install \
 cmake flex bison \
 git wget unzip \
 freeglut freeglut-devel \
 libXi-devel libXi \
 libXmu libXmu-devel \
 libXrandr libXrandr.devel \
 libXpm-devel libXft-devel \
 libXres-devel libXcomposite-devel \
 redhat-lsb
</code></pre>

#### ANTs:
<pre><code>
git clone git://github.com/stnava/ANTs.git
mkdir antsbin && cd antsbin
cmake ../ANTs/
make -j 4
echo 'export PATH=$PATH:$HOME/antsbin/bin' >> ~/.bash_profile
</code></pre>

#### DCMTK:
<pre><code>
git clone https://github.com/commontk/DCMTK.git
cd DCMTK
./configure
make all
sudo make install-all
</code></pre>

#### minc-toolkit:
<pre><code>
wget http://packages.bic.mni.mcgill.ca/minc-toolkit/RPM/minc-toolkit-1.0.08-20160205-CentOS_6.7-x86_64.rpm
sudo rpm -Uvh minc-toolkit-1.0.08-20160205-CentOS_6.7-x86_64.rpm

echo 'export PATH=$PATH:/opt/minc/bin' >> ~/.bash_profile
echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/minc/lib' >> ~/.bash_profile
echo 'export PERL5LIB=$PERL5LIB:/opt/minc/perl' >> ~/.bash_profile
source ~/.bash_profile
</code></pre>

#### mni-models:
<pre><code>
mkdir mni-models && cd mni-models
wget http://packages.bic.mni.mcgill.ca/mni-models/icbm152/mni_icbm152_lin_minc1.zip
unzip mni_icbm152_lin_minc1.zip
rm README COPYING
mincblur -fwhm 2 -gradient icbm_avg_152_t1_tal_lin.mnc icbm_avg_152_t1_tal_lin_2
mincblur -fwhm 4 -gradient icbm_avg_152_t1_tal_lin.mnc icbm_avg_152_t1_tal_lin_4
mincblur -fwhm 8 -gradient icbm_avg_152_t1_tal_lin.mnc icbm_avg_152_t1_tal_lin_8
mincblur -fwhm 16 -gradient icbm_avg_152_t1_tal_lin.mnc icbm_avg_152_t1_tal_lin_16
cp icbm_avg_152_t1_tal_lin_mask.mnc icbm_avg_152_t1_tal_lin_2_mask.mnc
cp icbm_avg_152_t1_tal_lin_mask.mnc icbm_avg_152_t1_tal_lin_4_mask.mnc
cp icbm_avg_152_t1_tal_lin_mask.mnc icbm_avg_152_t1_tal_lin_8_mask.mnc
cp icbm_avg_152_t1_tal_lin_mask.mnc icbm_avg_152_t1_tal_lin_16_mask.mnc
cd ..
sudo mv mni-models /opt/minc/share/
</code></pre>

#### Add links and install ctime.pl:
<pre><code>
sudo mkdir /opt/minc/etc/mni_autoreg
sudo ln -s /opt/minc/etc/mritotal.default.cfg /opt/minc/etc/mni_autoreg/mritotal.default.cfg
sudo ln -s /opt/minc/etc/mritotal.icbm.cfg /opt/minc/etc/mni_autoreg/mritotal.icbm.cfg
get http://perl5.git.perl.org/APC/perl-5.10.x/lib/ctime.pl
sudo mv ctime.pl /usr/share/perl5/
</code></pre>

## Trouble shooting
If libhdf5.so.6 is missing, but libhdf5.so.7 (or higher) is installed, do the following, or similar depending on installation location:
<pre><code>
sudo ln -s /usr/local/lib/libhdf5.so.7 /usr/local/bic/lib/libhdf5.so.6
</code></pre>
