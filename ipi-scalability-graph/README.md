Usage
======= 

The following are instructions to use the workflow from the hirshlab group.

Install Snakemake
------------------

First make sure that a module of python3 is loaded. E.g.:

```bash
module load python/python-anaconda_3.7
```

Then isntall snakemake (should be done only once per user):

```bash
pip install --user snakemake
```

Sanity check
-------------- 

Check that the workflow works on local node for small run:	
Edit `config.yaml` to use only 16 bosons:

```yaml
#boson_numbers: [16, 32, 64, 128, 256, 512, 1024]
boson_numbers: [16]
```

setup the environment to use the modified i-pi version:
```bash
source "/hirshblab/data/yotamfe1/mi-pi/env.sh"
```
This step will be eliminated when those env settings will be translated to a modulefile.

Run the workflow using:  
```bash
snakemake --cores 1
```
This should create the file `boson_scalabilty.png` after about 2 minutes


Next, perform the sanity check by sending the same workflow to the queue.  Note that no changes in the Snakemake file are required.  

```bash
git clean -x -f -d . # to clean all outputs from the previous check

snakemake --cluster 'qsub -q hirshb -l ncpus={threads}' --cluster-cancel 'qdel' -j 7&
```

To avoid the lenghty command line, the profile functionality may be used.

First, install the supplied profile directory:
```bash
mkdir ~/.config/snakemake
cp -rfp ../.config/hirshblab ~/.config/snakemake
```
Then the above snakemake command line turns into:
```bash
snakemake --profile hirshblab
```

Setup Snakemake syntax highlighting for VIM
--------------------------------------------

```bash
mkdir -p ~/.vim/syntax
wget https://github.com/snakemake/snakemake/raw/main/misc/vim/syntax/snakemake.vim
mv snakemake.vim ~/.vim/syntax

mkdir ~/.vim/ftdetect
wget https://github.com/snakemake/snakemake/raw/main/misc/vim/ftdetect/snakemake.vim
mv snakemake.vim ~/.vim/ftdetect
```

