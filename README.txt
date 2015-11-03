Data for Cursons et al. (2015). Spatially-transformed fluorescence 
image data for ERK-MAPK and selected proteins within human epidermis. 
GigaScience. Submitted Aug 2015.
	doi: to-be-known

For further information please contact:
The University of Melbourne Systems Biology Laboratory:
	Edmund Crampin - edmund.crampin@unimelb.edu.au
	Joe Cursons - joseph.cursons@unimelb.edu.au


Last Modified 01/09/15
-----------------------------------------------------------------

::: Data Overview:

For a more comprehensive description, please refer to the GigaScience
Data Note published with these data (cited at the top of this README).

The folders/directory structure for the immunofluorescence 
single-target data (detailed below) are as follows:
	/image/<folder>			**IF Target**
		CALM			Calmodulin
		ERK			ERK-1/2, MAPK3/1
		ERK_ph			phospho-ERK-1/2 (pT183/pY185)
		FOSC			c-Fos
		FRA2			Fra-2
		ITGB1			Integrin beta 1, CD29
		ITGB4			Integrin beta 4
		JUNB			Jun-B
		JUNC			c-Jun
		K10			Keratin-10
		K14			Keratin-14
		MEK			MEK-1/2, MAPKK1/2
		MEK_ph			phospho-MEK-1/2 (pS218/pS222)
		RAF			Raf-1, c-Raf
		RAF_ph			phospho-Raf-1 (pS338)
		SFN			Stratifin, 14-3-3 sigma

There are also a small subset of co-labelled data:
- /image_multi/RAF_ph+MEK_ph+ERK_ph:
	- containing Alexa 350 conjugated phospho-Raf; Alexa 555
		conjugated phospho-MEK and Alexa 488 conjugated
		phospho-ERK simultaneously labelled, together with
		appropriate cross-control and secondary control 
		images.
- /image_multi/MEK_ph+K14:
	- containing a small subset of Alexa 488 conjugated
		phospho-MEK and Alexa 555 conjugated K14; note that
		these data do not have control images.

There is also a collection of processed image data; containing
sample locations for image data quantification; tissue layer masks;
intermediate and fully-processed output (as both MATLAB .mat data
files and .csv files).
- /processed/<target> (as above)
* note that the processed files for the multi-labelled data are
	contained within the /image_multi folder

-----------------------------------------------------------------

::: Relative Folder Structures - Raw Image Data:

Image data are sorted by the immunofluorescence target (listed below),
and then by patient (Pat_x; 1-3), establishing the following folder
structures:
e.g. for Patient One, Keratin 10:
	/data/image/K10/Pat_1/image_data
e.g. for Patient Three, phospho-MEK (pS218/pS222):
	/data/image/MEK_ph/Pat_3/image_data

Within each folder, there is a .txt file with the microscope settings, 
such as:
	- pixel resolution
	- imaging gain and offset (these were manually adjusted between
		image sets due to the time associated with collecting all
		the data, but consistent with the corresponding 'secondary 
		control images' (described below)
	- fluorescence spectra (wavelengths) captured


The image data come with corresponding secondary control images (/sec_ctrl/), 
collected using the same microscope settings and made available for users that
wish to correct for background fluorescence. In some cases, a 'high gain'
image is also taken to indicate the position of the epidermis within the
image.
e.g. for Patient One, Keratin 10:
	/data/image/K10/Pat_1/sec_ctrl
e.g. for Patient Three, phospho-MEK (pS218/pS222):
	/data/image/MEK_ph/Pat_3/sec_ctrl
*NB*: Unfortunately the secondary control data for patient 3 c-Fos and ITGB1 
	have been lost.

The number in the image data prefix (and the name of the corresponding .txt 
file) is irrelevant, it simply reflects the order in which image data were
collected on a given day.
Images named in the format:
<#>_Image<###>_ch<##>.tif - reflect individual images taken during gain and
				offset calibration
<#>_Series<###>_z<###>_ch<##>.tif - reflect image stacks (through the z-
					direction), the data used for
					image processing and segmentation.

The fluorescence channel is recorded in ch<##>, starting at ch00. This is 
ch00 (green; collected around 520 nm) for all data except phospho-ERK Pat2.
These image data are still labelled ch00 however, as they were collected 
during preliminary experiments that trialled Alexa-555.

Some co-labelling experiments were also performed, however there were no
consistent data collection efforts as it was not practicable to address
large numbers of antibody combinations. These include:
	/data/image_multi/RAF_ph+MEK_ph+ERK_ph
	/data/image_multi/MEK_ph+K14
* Note that the co-labelled K14 and phospho-MEK data show a 'weird'
	tissue morphology, in that there are a large number of
	dermal protrusions, creating 'dermal bead' structures
	(rings of basal cells separated ['more suprabasal'] from 
	the identifiable tissue basal layer). This is relatively
	interesting from the perspective of considering
	phospho-MEK-1/2 regulation within basal keratinocytes.

-----------------------------------------------------------------

::: Relative Folder Structures - Processed Data:

The processed data, including segmentation masks and sample locations, are
stored under /data/processed/<IF target>/Pat<#>, in a similar fashion to the
raw image data.

Within these folders, there are tissue layer boundary segmentation masks, 
stored as binarised TIFF files:
	im_basal_memb_z<#>.tiff
	im_bound1_z<#>.tiff
	im_bound2_z<#>.tiff
	im_outer_memb_z<#>.tiff

Together with these, is a combined tissue segmentation mask, 
im_layers_z<#>.tiff, which contains data values stored within the 
corresponding epidermal tissue layers:
	0 = outside epidermis (including dermis), or stratum corneum
	1 = basal layer
	2 = spinous and granular layers
	3 = transitional layer

There is also a /C/, /N/, and/or /M/ subfolder, corresponding to the
cytoplasm, nucleus or plasma membrane. These folders contain the sample
location data and the corresponding sampled fluoresecence data, stored as:
.mat files - MATLAB data files that can be loaded as intermediate output
		when executing the MATLAB scripts to improve run-time;
	sample_loc.mat - sample locations
	sample_ana.mat - sampled fluorescence data
.TIF file - binarised image file containing sample locations
	im_cytoSamples_z<#>.tiff
	im_nucSamples_z<#>.tiff
.CSV files - sampled fluorescence intensity data at the specified sample
		locations (using the specified sampling kernel);
	intensity_xy.csv - signal intensity data at specified (x,y) co-
				ordinates within the image data,
		** single column
	intensity_dnorm.csv - signal intensity data at the specified 
				normalised-distance position,
		** multiple (variable) columns

For a subset of data, individual cell masks for the cytoplasm and 
nucleus, as used to generate Fig. AF4.1. These are generally labelled
in a consistent manner as the image from which they were produced,
for example:
	6_Series004_z004_ch00.tif (a phospho-MEK Patient Two data file)
	--> 6_Series004_z004_nuc.tif
	    6_Series004_z004_cyto.tif
	under \data\processed\MEK_ph\Pat_2\
Also:
	11_Series033_z019_ch0<1-3>.tif (co-labelled pRaf/pMEK/pERK data)
	--> 11_Series033_z019_cyto.tif
	    11_Series033_z019_nuc.tif

In many cases, these masks were created using GIMP with multiple layers
to control opacity, and the path tool to perform segmentation. In these
cases, the GIMP .xcf files are provided for users to easily extend or
modify the segmentation.

For a small subset of the data, there are also '3D whole cell' 
segmentation data, used to derive Fig. S10 in [Cursons et al. (2015),
BMC Systems Biology], where we investigated changes in nuclear-to-
cytoplasmic volume ratio over the depth of the epidermis.
These data are labelled as 'CytoSegVol.mat' and 'NucSegVol.mat' files,
found within the /processed/<target>/Pat_<#>/ folders.

-----------------------------------------------------------------

::: Relative Folder Structures - Image Processing Scripts:

A collection of MATLAB scripts and functions are also provided for
image processing, and transformation of the image data onto a 
normalised distance co-ordinate, as described in:
	Cursons et al. (2015). Regulation of ERK-MAPK signaling 
		in human epidermis. BMC Systems Biology.
		doi: 10.1186/s12918-015-0187-6

For a more detailed description of these image processing scripts,
please refer to the PDF documentation hosted upon the GitHub
project page:
	http://github.com/uomsystemsbiology/epidermal_data

These scripts should be extracted to the base sub-folder: 
	/code/
Functions which are used by these scripts are in the sub-folder:
	/code/functions/

These scripts (and dependent functions) require a number of MATLAB
toolboxes, including the Image Processing Toolbox, the Curve Fitting
Toolbox and the Neural Network Toolbox.

The scripts which will form a useful starting point for users
looking to modify these analyses, include:
	sample_to_loessDiscCentroids.m - a script which reads in
		the sample location data, extracts pixel intensities
		from the raw image data, maps to the normalised
		distance co-ordinate, performs loess smoothing, and
		outputs .csv/.tiff at various stages.
	recreate_pMEK_heterogeneity_fig.m - a script which will 
		reproduce the plots shown in Fig. AF4.1 (note that
		this output was tidied in GIMP to produce the
		final figure).

These scripts will be expanded over the following weeks (during
September) to ensure a variety of analyses are available to users.

A Virtual Reference Environment is also available to users who do
not have access to MATLAB, but can satisfy the MathWorks Academic
Usage Agreement:
	http://github.com/uomsystemsbiology/epidermal_data_vagrant

For further details on Virtual reference environments, please
refer to:
	Hurley, DG, Budden, DB, & Crampin, EJ. (2014). Virtual 
		Reference Environments: a simple way to make 
		research  reproducible. Briefings in bioinformatics.
		doi:10.1093/bib/bbu043
	http://uomsystemsbiology.github.io/research/reproducible-research/

-----------------------------------------------------------------