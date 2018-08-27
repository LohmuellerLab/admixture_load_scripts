: <<'COMMENT'
This script sets up a conda environment for msprime (python3)
The lines that are commented out are for reference if you want to do optional things
COMMENT

#install tools if oyu don't have them
#sudo apt-get install -y bzip2 git wget

cd $HOME

#silent install of miniconda
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash $HOME/Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3

#source conda, adds conda binary folder to your paths
. $HOME/miniconda3/etc/profile.d/conda.sh

#add channels
conda config --add channels conda-forge
conda config --add channels bioconda

#create environment for msprime and activate it. everything you install in the 
#env is local and doesn't require sudo. all the installed packages,programs,libs, etc.
#will only be available in when the specific env is activated
#https://anaconda.org/anaconda/repo
#https://anaconda.org/conda-forge/repo
#https://bioconda.github.io/recipes.html

#you specify the python but note that it will up/downgrade versions as you install
#so be careful with python packages so you don't break anything, e.g. install python2
#packages into a separate python2 env
conda create -n msprime python=3
conda activate msprime

#install msprime stuff
conda install -y -c conda-forge msprime

#install python package from github
pip install git+https://github.com/tskit-dev/pyslim.git

#export conda environment to yaml file that can be shared between users
#conda env export > msprime.yml

#close conda environment
#conda deactivate 

#add source lines to bashrc so conda will be automatically available when you login
echo "" >> ~/.bashrc
echo ". ${HOME}/miniconda3/etc/profile.d/conda.sh" >> ~/.bashrc
. ~/.bashrc

#the statement ". ${HOME}/miniconda3/etc/profile.d/conda.sh" should be written to 
#job submission scripts. then  activating a local conda env should be no prob