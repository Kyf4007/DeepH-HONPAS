# DeepH-HONPAS
Equivariant neural network deep learning Hamiltonian code supported by HONPAS

# Introduction
DeepH-HONPAS is a computational package designed for electronic structure calculations. It integrates DeepH (https://github.com/mzjb/DeepH-pack?tab=readme-ov-file) with Hefei Order-N Packages for Ab initio Simulations (HONPAS, https://github.com/xmqin/HONPAS?tab=GPL-3.0-1-ov-file) to provide efficient and accurate DFT results. This code support DeepH and DeepH-E3.

# Installation

To prepare the environment on a Windows/Linux/MacOS machine, directly use:
```
conda create -n DeepH-HONPAS python=3.11
```
For denpendencies, both 'pip' and 'conda' are possible. The required packages can be found in https://github.com/mzjb/DeepH-pack?tab=readme-ov-file. The following packages are mandatory:
- NumPy
- SciPy
- e3nn
- pymatgen
- h5py
- pathos
- psutil
- PyTorch
- PyTorch Geometric

The version labels are not strict, but here we provide an example on a personal computer. The following is the hardware and system configuration of the personal computer used during the development and testing of this project. This information is provided for reference only; the code should be compatible with a wide range of systems that meet the software requirements.
CPU: Intel® Core™ i9-10850K @ 3.60GHz
RAM: 64 GB
GPU: NVIDIA GeForce RTX 4090 (24 GB VRAM)
Operating System: 64-bit Windows 10
CUDA Version: 12.8

Here are the installation commands for DeepH-HONPAS environment:
```
pip install numpy==2.2.6 -i https://pypi.org/simple
pip install scipy==1.15.3 -i https://pypi.org/simple
pip install e3nn -i https://pypi.org/simple
pip install pymatgen -i https://pypi.org/simple
pip install h5py -i https://pypi.org/simple
pip install pathos -i https://pypi.org/simple
pip install psutil -i https://pypi.org/simple
```
For the installation of PyTorch, please refer to the official website and select the suitable command for your system. For example, to install the latest version on a compatible GPU with CUDA installed:
```
pip3 install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu126
```
The version of PyTorch Geometric needs to be compatible with PyTorch version. To install PyTorch Geometric, additional libraries are necessary: 'pyg-lib', 'torch-cluster', 'torch-sparse', 'torch-scatter', 'torch-spline-conv'. 
Additionally, if electronic band structure calculation is required ('task = [5]' in the inference step), a Julia interpreter is required. Please refer to [https://github.com/mzjb/DeepH-pack] to install Julia. Note that Julia code works only in band calculation, and the predicted Hamiltonian already comes in handy after 'task = [4]', thus the Julia interpreter can be installed elsewhere.

When the environment activated, install DeepH-HONPAS with:
```
git clone
cd DeepH-HONPAS
pip install .
```

HONPAS is required for dataset generation. Please refer to https://github.com/xmqin/HONPAS?tab=GPL-3.0-1-ov-file for HONPAS installation and usage.
If you only want to try the machine learning code, the dataset for twisted bilayer graphene (TBG) and twisted bilayer molybdenum disulfide (TBMoS2) are also available at Zenodo [https://doi.org/10.5281/zenodo.14979348].

# Usage

For preprocess, train, and inference, please edit the .ini files in the deeph/scripts/ directory and run the following commands. Here are the commands for using DeepH with HONPAS.
To prepare the dataset, we recommend using the preprocess.ini in the deeph/scripts/ directory. Set the paths for 'raw_dir' and 'processed_dir' and run:
```
deeph-preprocess --config preprocess.ini
```
To train the model, set the paths for 'graph_dir', 'save_dir' and 'raw_dir' for the graph files, training log files, and processed dataset files. Begin training by:
```
deeph-train --config train.ini
```
To inference, set the paths for 'work_dir', 'OLP_dir', and 'trained_model_dir' for prediction results, overlap matrix files, and trained model files (training log in the train step). To get the predicted band structures, set 'task = [1, 2, 3, 4, 5]' and run:
```
deeph-inferenece --config inference.ini
```

For DeepH-E3, the main difference is the neural network structure used in the training step. Here are the commands for using DeepH-E3 with HONPAS.
In DeepH-E3, the dataset process and graph generation steps are merged into a single step. After editing the base.ini in the DeepH-E3/scripts/, run:
```
python deephe3-preprocess.py base.ini
```
After setting the proper paths in .ini files, the training and inference tasks can be done by:
```
python deephe3-train.py train.ini
python deephe3-eval.py eval.ini
```
Note that after excecuting 'deephe3-eval', the predicted Hamiltonians are available in the working directory, but 'deeph-inference' is still needed to calculate the band structure. Thus set the 'task = [5]' in inference.ini, set the 'work_dir' and run 'deeph-inferenece --config inference.ini' to sparse matrix calculation for band structures.

