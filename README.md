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

## Dependencies
### Add MNI to source list to install autoreg
<pre><code>
echo 'deb http://packages.bic.mni.mcgill.ca/ubuntu-maverick/ ./' | sudo tee -a /etc/apt/sources.list
sudo apt-get update
</code></pre>

### Install minc tools, DCMTK, ANTS and other compile dependencies</b></p>
<pre><code>
sudo apt-get install libnetcdf-dev libhdf5-dev hdf5-tools \
   libdcmtk2-dev libdcmtk2 dcmtk minc-tools ants \
   mni-autoreg mni-models-icbm152-lin bicpl libmni-perllib-perl \
   g++ cmake libx11-dev bicpl libmni-perllib-perl
   </code></pre>

## Installation (currently only for linux)
<pre><code>
git clone --recursive git://github.com/claesnl/RESOLUTE_AC.git RESOLUTE_AC
cd RESOLUTE_AC/linux
mkdir build && cd build
cmake ..
make && sudo make install
</code></pre>

## Trouble shooting
If libhdf5.so.6 is missing, but libhdf5.so.7 (or higher) is installed, do the following, or similar depending on installation location:
<pre><code>
sudo ln -s /usr/local/lib/libhdf5.so.7 /usr/local/bic/lib/libhdf5.so.6
</code></pre>
