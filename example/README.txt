
This is the example for the Hamiltonian prediction of twisted bilayer graphene (twist angle 7.34°).
This example is tested under DeepH-HONPAS code.
To run this example, please install the environment and DeepH-HONPAS code, refering to [https://github.com/Kyf4007/DeepH-HONPAS].

1. Preprocess

a. Please download the dataset named "Bilayer_graphene_PBE.zip" from [https://zenodo.org/records/14979348] 
and extract this .zip file in the directory: DeepH-HONPAS\example\Twisted_bilayer_graphene
This will generate a directory named "official_DZ" which contains DFT results of 300 bilayer graphene structures.

b. Change directory to: DeepH-HONPAS\example\Twisted_bilayer_graphene\scripts.
Activate the installed environment if not, and modify the following paths in preprocess_example.ini:
raw_dir = \absolute\path\to\official_DZ
processed_dir = \absolute\path\to\processed_dir

c. Type the command:
deeph-preprocess --config preprocess_example.ini

2. Train

a. In the directory: DeepH-HONPAS\example\Twisted_bilayer_graphene\scripts, modify the train_example.ini according
to your machine:
graph_dir = G:\\DeepH-HONPAS-test\\DeepH-HONPAS-main\\example\\Twisted_bilayer_graphene\\graph_dir
save_dir = G:\\DeepH-HONPAS-test\\DeepH-HONPAS-main\\example\\Twisted_bilayer_graphene\\train_dir
raw_dir = G:\\DeepH-HONPAS-test\\DeepH-HONPAS-main\\example\\Twisted_bilayer_graphene\\processed_dir

b. Type the command:
deeph-train --config train_example.ini

c. The training process can be stopped with ctrl+c, and the best model is in the train_dir.

3. Inference

a. Download the overlap matrix file of 7.34° TBG, "TBG-7.34-DZ-PBE.zip", from [https://zenodo.org/records/14979348] (This file has already been downloaded)
Extract this file in the OLP_dir, and get TBG-7.34-DZ directory. This directory contains the overlap matrix of large structure,
but not the converged Hamiltonian.
Edit the path in the inference_example.ini, make sure that trained_model_dir point to the directory named by date:
work_dir = G:\DeepH-HONPAS-test\DeepH-HONPAS-main\example\Twisted_bilayer_graphene\work_dir\TBG-7.34
OLP_dir = G:\DeepH-HONPAS-test\DeepH-HONPAS-main\example\Twisted_bilayer_graphene\OLP_dir\TBG-7.34-DZ
trained_model_dir = ["G:\\DeepH-HONPAS-test\\DeepH-HONPAS-main\\example\\Twisted_bilayer_graphene\\train_dir\\2025-05-21_22-21-01"]
sparse_calc_config = G:\DeepH-HONPAS-test\DeepH-HONPAS-main\example\Twisted_bilayer_graphene\sparse_calc\info.JSON

b. In the script directory, type:
deeph-inference --config inference_example.ini

c. Please make sure julia is functional.
In the work_dir, for example in "G:\DeepH-HONPAS-test\DeepH-HONPAS-main\example\Twisted_bilayer_graphene\work_dir\TBG-7.34",
a band.mat file will contain the electronic bands, which can be read in MATLAB.
Also, the standard HONPAS band output, honpas.band, will also available in work_dir.
