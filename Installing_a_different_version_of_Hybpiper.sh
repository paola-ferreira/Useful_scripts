# How to Install an Older Version of HybPiper
# I spent quite some time figuring out how to install an older version of HybPiper (version 2.0.2). Here’s a straightforward way to manage it:

# 1. Download and install Miniconda (if not installed yet)
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
# Follow the prompts to complete installation, then restart your shell or open a new terminal

# 2. Install mamba in the base environment (faster conda alternative)
conda install -c conda-forge mamba

# 3. Create a new environment named 'hybpiper' and install HybPiper
mamba create -n hybpiper hybpiper

# 4. Activate the new environment
eval "$(mamba shell hook --shell bash)"
mamba activate hybpiper

# 5. Install specific HybPiper version 2.0.2 directly from GitHub
pip install git+https://github.com/mossmatters/HybPiper.git@2.0.2

# 6. Verify installation and check version
hybpiper -v
