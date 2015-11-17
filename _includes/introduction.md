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
RNA interference (gene silencing) or miRNA screens (inhibitors or mimicks). Each of these screen types demands a different strategy for down-stream analysis. 
All three types of screens are supported in HiTSeekR and due to differences in the downstream analysis, the first step in a HiTSeekR analysis is always to choose a screen type as illustrated here:


![The choice of the screen type is always the first step in HiTSeekR]({{ "/images/tutorial/select_screen.png" | prepend:site.baseurl }} "The choice of the screen type is always the first step in HiTSeekR")

We will elaborat on all three screen types in this tutorial. However, since the first four steps are identical, they will only be shown exemplary for the first type of screen focused on gene silencing.    

## RNA interference / CRISPR CAS9 screens (Genes)
<a name="siRNA"></a>

### Step 1: Input Data

After selecting _genes_ on the initial page, the input data tab is shown:

![For the input data you can choose between demo data sets, uploading your own data, or downloading a dataset via its PubChem Assay ID]({{ "/images/tutorial/input_options.png" | prepend:site.baseurl }} "For the input data you can choose between demo data sets, uploading your own data, or downloading a dataset via its PubChem Assay ID")

Here, the user is offered three choices:

1. Selecting one of the available demo datasets for this screen type

2. Uploading a new dataset in a comma, semicolon or tab separated format

3. Selecting and downloading a dataset from PubChem via its assay id (AID)

For now, we will select _RNAi screen for vorinostat resistance genes_ as input. This screen was originally obtained from PubChem (https://pubchem.ncbi.nlm.nih.gov/bioassay/743454). This action will change the view:

![Options after selecting a screen]({{ "/images/tutorial/input_options2.png" | prepend:site.baseurl }} "Options after selecting a gene")

At the bottom of the page, a table is shown that represents a raw representation of the input data. This allows the user to check if the data was 
read in correctly and to identify the content of each column. 
The user now has to decide whether the measured values need to be log2 transformed (stabilizes variance between lower and higher signal values), 
if the B-score should be calculated (This is a computationally expensive score suitable to counter position related effects), and, most importantly,
the user can also assign columns of the input data to properties in HiTSeekR. This step can be omitted for demo data since the column assignment is 
pre-configured. For custom uploaded and PubChem downloaded screens, however, this step is necessary before the analysis can be continued:

![Assigning columns in the raw data]({{ "/images/tutorial/input_options3.png" | prepend:site.baseurl }} "Assigning columns in the raw data")

Please take
your time to do the column assignment right, because experience shows that most of the problems with using HiTSeekR are caused by incorrect column assignments. A few notes:
 
* For the well position, three different formats are accepted. Row / column indices, alphanumeric well names, e.g. A01, and absoulte well indices, e.g. (1-384)
* For each screen type, we support the most common identifiers, which HiTSeekR will try to map to its default accession type. We advise users to check if this mapping is successful, since outdated identifiers
cannot always be mapped like ambigious miRNA identifiers that do not yet follow the 3p and 5p notation. 
* Several measurement, replicate and experiment columns can be selected. The background of this is that in many data sets this information is not fully normalized. This header, for example is not normalized:

| Sample | Experiment1 | Experiment2 | Experiment3 | Accession |
| :----: | :-------:   | :---------: | :---------: | :-------: |
|...		|...			  |...			 |...		      | ...			|

In contrast, this header is normalized:

| Sample | Experiment | Accession |
| :----: | :--------: | :-------: |
|...     |Experiment1 |...        |
|...     |...         |...        |
|...     |Experiment2 |...        |
|...     |...         |...        |
|...     |Experiment3 |...        |


If multiple columns of experiments, readouts, or measurements are selected, HiTSeekR will attempt to construct a normalized table. 
Note that this multiple column feature is currently restricted to one of those three types.
After choosing the appropriate columns the user is encouraged to click on _Process raw data_ to continue. After successful processing, new options will appear in the tab bar on top of the page. 
In addition, continue buttons like the one highlighted here are supposed to guide the user through the process:

![After the processing of the raw data, new options become available]({{ "/images/tutorial/continue_to_qc.png" | prepend:site.baseurl }} "After the processing of the raw data, new options become available")

### Step 2: Quality Control

In the next step of the analysis, various plots are produced to illuminate various aspects of the screening quality. The following plot, for instance, gives an impression of the signal distribution over
the entire screen. Each box plots here depicts the signal spread of one plate, allowing it to easily spot differences of plate means and signal variance. Please note the checkbox _show help text_, which, when clicked, will provide 
additional information on each plot to help the user with the interpretation:

![Quality control]({{ "/images/tutorial/continue_to_qc.png" | prepend:site.baseurl }} "Quality control")

The following plots are available:

* Plate Signal Variation: Whole screen box plots of each plate and replicate, divided into different readouts and experiments.
* Row and Column Effect: Bar plot of avarge signal of each row and column over the entire screen. Helps to spot positional bias.
* Control Signal Spread: Depicts in box plots the signal variation of the controls (if included).
* Control Separability: Depicts strictly standardized mean scores (SSMD, Zhang et al.) scores for assessing the separability of positive and negative controls (if given).
* Replicate Correlation: Pairwise comparison of screening replicates including linear regression and R2 value to assess replicate correlation. 

The main purpose of the quality control step is to learn about bias in the raw screening data. This allows the user to make informed decisions on the choice of normalization method in the next step.

### Step 3: Normalization

In this step, various normalization methods can be compared by investigating their effect on the data. First, the user can select one of the following methods:

* centered by mean: Divide each plate by its mean signal value 
* centered by median: Divide each plate by its median signal value
* z-score: Subtract for each plate the mean value and divide by the standard deviation
* robust z-score: Subtract for each plate the median value and divide by the median absolute deviation
* percentage of negative control: Use the negative controls as reference point (1).
* normalized percentage inhibition: Use positive and negative control as reference points (positive 100, negative 0).

The plots will then be updated to reflect the chose normalization. As an example, we here show the plate signal variation plot again, which, with the selected z-score, guarantees that the plate means are now aligned:

![Normalization]({{ "/images/tutorial/normalization.png" | prepend:site.baseurl }} "Normalization")

In addition, to a downloadable table with all of the normalized values, the following plots are available to study the effect of a normalization method:

* Plate Signal Variation: Same as above
* Whole Screen Scatterplot: Similar to the plate signal variation plot, but shows each sample as individual point. For smaller screens an interactive javascript version is rendered.
* Signal Distribution: This plot shows the signal density plot for each replicate. Here, normalization effects can be assessed very effectively by checking how well replicates agree.
* QQ plot: The qq plot shows if the data is normally distributed. In this case most of the data points will form a line from which only active samples will deviate at the beginning and end.
* Plate Viewer: Here, individual plates can be investigated closely by an interactive scatter plot and heatmap. 
* Replicate Correlation: Same as above

When the user feels confident with the different normalization methods and their effect, we continue with the identification of active samples, also called hits. 

### Step 4: Hit Discovery

On this page, the user might feel a bit overwhelmed with the different options at first, so let's go through them one by one, starting with the select boxes:

* Select experiment: The user can select one or seveal experiments to check for hits
* Select readout: The user can select one or several readouts to be included
* Normalization: Here, the normalization should be selected that appeared most appropriate after quality control and the normalization step before. 
* Hit discovery method: There are several ways to identify hits:
  * Fixed cutoffs: This is the most primitive hit detection method in which the user defines an upper and lower threshold to be considered as hard cutoff.
  * Standard deviation: Probably the most common hit detection method in which the general variability of the screening data is used as a measure of activity. A common threshold here is a margin of 3 standard deviations.
  * Median absolute deviation: Similar to standard deviation method, but more robust due to the use of the median instead of the mean.
  * Inter-quartile range: Another variant in which the difference of the 3rd and 1st quartile are used instead to assess variability.
  * SSMD: The strictly standardized mean difference proposed by Zhang et al. is a robust alternative to multiple t-tests between the sample and the controls, which is not affected by the number of replciates / samples. A SSMD value of +-3 is considered good, a SSMD value of +-6 is considered excellent. Naturally, this method only works if replicates are available.
  * Bayesian statistics: The probably most sophisticated hit detection method takes plate as well as screen-wide statistics into account as priors to assess sample activity as a posterior probability. It is the only method that allows setting a threshold based on a p-value and works best on raw signal values. It depends on the availablity of negative control samples on each plate.
Similar to the normalization step, we encourage users to try out different strategies and thresholds. Comparing the results will reveal what works best for a particular dataset.
* Effect: This is interesting if the user is interested in only one effect direction, e.g. only promotors or suppressors.
* Categories: HiTSeekR depicts each hit with one of the following three categories:
  * Promotor (blue): The sample shows increased signal, e.g. proliferation effect
  * Suppressor (red): The samples shows decreased signal, e.g. killing effect
  * Included (yellow): The sample has been manually included in the hit list by the user
  
Below the options panel, a tabular hit list is shown. This table, like all tables in HiTSeekR, can be downloaded (button below the table), sorted after arbitrary columns and searched either on the entire table or on a specific column. 

The hits plot tab provides a graphical representation of the screening hits in a scatter plot:

![Hits plot]({{ "/images/tutorial/hits_plot.png" | prepend:site.baseurl }} "Hits plot")

Finally, the Heatmap is a graphical representation of the entire screen. Hits are depicted by small black arrows. 
This plot allows the user to quickly spot if hits accumulate on a specific plate or a general plate location (e.g. mostly in the first row or column etc.), which might indicate problems with the selected normalization or hit identification strategy:

![Whole screen heatmap with hits indicated by arrows.]({{ "/images/tutorial/heatmap.png" | prepend:site.baseurl }} "Whole screen heatmap with hits indicated by arrows.")

When the user has successfully selected hits we are ready for downstream analysis. Therefore, additional tabs are unlocked as soon as the hit discovery tab is selected.

### Step 5: Downstream Analysis

#### Gene sets

#### De novo network enrichment

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