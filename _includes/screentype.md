## Introduction
<a name="Introdcution"></a>

High-throughput screening (HTS) is characterized by experimentation on a large or even ultra-large (> 100,000 experiments) scale, which is achieved through a high degree of automation with regards to pipetting and microtiter plate handling as well as through miniaturization of experimental procedures. 
Data generated in HTS are subject to variation caused by different sources and therefore requires appropriate normalization. In general, the following steps are taken:

1. **Input Data:** The columns of the input data need to be assigned correctly, for HiTSeekR to proceed with the analysis.
2. **Quality Control:** Sources of variation are identified and the quality of replicates and controls are assessed.
3. **Normalization:** Bias in the data is countered by a suitable normalization method. Control based and plate based normalization methods are differentiated.
4. **Hit Discovery:** A hit detection method, either based on controls or on plate statistics, is used to identify active samples in the screen.
5. **Down-stream analysis:** The hit list generated in the step above is used as input for further analyses based on systems biology methods. The goal here is to learn what biological processes are affected rather than relying on individual observations that can also be due to outliers or off-target effects of the samples.

HiTSeekR guides the user through each of these steps and enables well informed decisions based on informative plots and direct method comparisons. 
One decisive advantage of the HTS technology is its versatility, i.e. almost any kind of cell based assay can in theory be adapted. 
Since the early days of the technology, where it was mainly used for small molecule screens (drug screens) in pharmaceutical companies, this flexibility has opened up new use cases, such as
RNA interference (gene silencing) or miRNA screens (inhibitors or mimicks). Each of these screen types demands a different strategy for down-stream analysis. All three types of screens are supported in HiTSeekR and will be elaborated on in this tutorial.
However, since the first four steps are identical, they will only be shown exemplary on an RNAi screen aimed at the discovery of vorinostat resistance genes obtained from PubChem (https://pubchem.ncbi.nlm.nih.gov/bioassay/743454).    

## RNA interference / CRISPR CAS9 screens (Gene silencing)
<a name="siRNA"></a>

### Step 1: Input Data

### Step 2: Quality Control

### Step 3: Normalization

### Step 4: Hit Discovery

### Step 5: Downstream Analysis

## microRNA inhibitor / mimicks screens
<a name="miRNA"></a>

### Steps 1-4

See RNA interference.

### Step 5: Downstream Analysis

## Small compound (drug) screen
<a name="compound"></a>

### Steps 1-4 

See RNA interference

### Step 5: Downstream Analysis