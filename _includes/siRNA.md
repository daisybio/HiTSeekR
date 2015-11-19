## RNA interference / CRISPR CAS9 screens (Gene silencing)
<a name="siRNA"></a>

Genes and their products are the basic unit of the cellular machinery. Systematically inhibiting gene function is thus an ideal tool for functional and drug studies based on cellular responses. 
This type of screen is supported in HiTSeekR and begins by selecting the appropriate screen type:

![The choice of the screen type is always the first step in HiTSeekR]({{ "/images/tutorial/select_screen_genes.png" | prepend:site.baseurl }} "The choice of the screen type is always the first step in HiTSeekR") 

### Step 1: Input Data

After selecting _genes_ on the initial page, the input data tab is shown:

![For the input data you can choose between demo data sets, uploading your own data, or downloading a dataset via its PubChem Assay ID]({{ "/images/tutorial/input_options.png" | prepend:site.baseurl }} "For the input data you can choose between demo data sets, uploading your own data, or downloading a dataset via its PubChem Assay ID")

Here, the user is offered three choices:

1. Selecting one of the available demo datasets for this screen type

2. Uploading a new dataset in a comma, semicolon or tab separated format

3. Selecting and downloading a dataset from PubChem via its assay id (AID)

For now, we will select _RNAi screen identifies Caspase4 as factor for TNFa signaling_ as input (obtained from http://mcb.asm.org/content/32/17/3372.abstract). This action will change the view:

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

In this step, the user is required to select a normalization and hit identification strategy to identify those samples that are considered active, i.e. that change the phenotype in the experiment. 

![Hit Discovery]({{ "/images/tutorial/hit_discovery.png" | prepend:site.baseurl }} "Hit Discovery")

Here, users might feel a bit overwhelmed with the different options at first, so let's go through them one by one, starting with the select boxes.

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

A number of check boxes on the left side provide additional options:

* Show sample position: Add the well location to the output table
* Show all computed scores: Add normalized values other than the one selected to the table
* Sample filter options: Provides users with the possiblity to use positive and negative, i.e. inclusion and exclusion filters using regular expression on the sample names.
* Differential screen: Enables users to select one readout or experiment as baseline screen that is considered as reference. Hits identified in the reference screen will be excluded from the remaining screens. The background of this is that often such reference readouts are included to separate specific from unspecific effects. This can be used, for instance, to exclude gene knockouts that are generally lethal in a particular cell line. 
  
Below the options panel, a tabular hit list is shown. This table, like all tables in HiTSeekR, can be downloaded (button below the table), sorted after arbitrary columns and searched either on the entire table or on a specific column. 

The hits plot tab provides a graphical representation of the screening hits in a scatter plot:

![Hits plot]({{ "/images/tutorial/hits_plot.png" | prepend:site.baseurl }} "Hits plot")

Finally, the Heatmap is a graphical representation of the entire screen. Hits are depicted by small black arrows. 
This plot allows the user to quickly spot if hits accumulate on a specific plate or a general plate location (e.g. mostly in the first row or column etc.), which might indicate problems with the selected normalization or hit identification strategy:

![Whole screen heatmap with hits indicated by arrows.]({{ "/images/tutorial/heatmap.png" | prepend:site.baseurl }} "Whole screen heatmap with hits indicated by arrows.")

When the user has successfully selected hits we are ready for downstream analysis. Therefore, additional tabs are unlocked as soon as the hit discovery tab is selected.

### Step 5: Downstream Analysis

The hit list in a gene silencing experiment such as RNA interference or CRISPR / CAS9 screens consists of a number of genes that are associated with a studied phenotype. While each of these hits might be
an interesting gene for follow-up experiments, limited resources usually do not permit to do this. Moreover, effects observed with individual genes might in some cases be attributed to off-target effects etc.
It is therefore of interest to consider the observed changes in phenotype on the level of systems biology. This provides functional annotation of, for instance, pathways and biological processes involved in these changes. 
In addition, it allows for new complex hypotheses to be developed and to be tested in future experiments. In HiTSeekR we distinguish two general types of systems biology analysis:
 
#### Gene sets

Gene sets are manually curated collections of genes that are annotated with biological properties. This includes gene ontology terms and pathway sources such as Reactome or KEGG. Two methods are integrated 
in HiTSeekR to test if a given hit list is significantly associated with given gene sets, namely gene set overrepresentation analysis, which relies on hypergeometric testing, and gene set enrichment analysis, which
is based on the computation of a running sum statistic and evaluates significance via generating permutations of the ranked genes. The minimum gene set size that can be adjusted here will help to reduce the
number of tests by ignoring smaller gene sets.

![Gene set overrepresentation and enrichment analyses]({{ "/images/tutorial/gene_set_analysis.png" | prepend:site.baseurl }} "Gene set overrepresentation and enrichment analyses")

Because a complete list of tested samples together with their achieved scores
are required for gene set enrichment analysis, this is only supported for gene silencing experiments, where information on all genes is available. Moreover, due to the computational demands of this method,
the number of permutations in the server is currently limited to 100. The user can select the results of the computed gene sets. Moreover, the results are divided into:

* Hypergeometric test results (gene set overrepresentation analysis)
* Gene set enrichment analysis 
* Gene sets that achieved a significant (corrected) p-value in both analyses 

#### De novo network enrichment

In contrast to gene set analysis, which depends on the quality of existing functional annotations, de novo network enrichment operates directly on large interaction networks that are not subject to this bias. 
The idea is to extract as large as possible sub-networks that are enriched with active genes. These can be the hits of a gene silencing screen or the miRNA or drug target genes, respectively, for the other two types of screens.
KeyPathwayMiner is a tool that enables efficient extraction of such subnetworks from arbitrary networks. In HiTSeekR, the web version of KeyPathwayMiner was integrated via its webservice API. This
has the advantage that computations are outsourced and that the user can continue working with HiTSeekR until the results are reported back.

To get started with a KeyPathwayMiner analysis you should first select a network and the number of exception genes you want to allow (between 0 and 3). This corresponds to the number of non-hit genes that 
KeyPathwayMiner is allowed to add to connect smaller sub-networks into larger ones. This parameter thus enables the user to take direct influence on the network enrichment process through an intuitive parameter.

When computation is complete, KeyPathwayMiner will show a slider that allows you to go through the top 20 subnetworks found. You can also select to draw a union graph composed of all 20 solutions in which
reoccuring genes are colored in a gradient from grey (appeared once) to bright green (was part in several or all solutions). All solutions can be exported in a cytoscape compatible SIF format and as
a tab delimited table. 

Note that pressing the "Start KeyPathwayMiner Analysis" button again will trigger a new analysis and erase existing results without a warning. 

![De novo network enrichment with gene silencing screens]({{ "/images/tutorial/de_novo_genes.png" | prepend:site.baseurl }} "De novo network enrichment with gene silencing screens")

