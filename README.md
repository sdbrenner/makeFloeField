# Floe-field generation

## Introduction

These MATLAB scripts/functions can be used to generate synthetic areal fields of sea ice floes with specified sea ice concentration (SIC) and floe size distribution (FSD) within a target domain/boundary.
These fields can be used, for example, as input files to sea ice discrete element models (DEMs) such as FloeDyn ([Rabatel et al., 2015](doi.org/10.1002/2015JC010909)) or for other analyses of surface heterogeneity in the marginal ice zone.

<img src='/Figures/floeField.png' width='400'>


## Files and dependencies

The file `makeFloeField.m` is a master script that specifies the target parameters and calls other functions.
Supporting functions are stored...

Floe shapes are selected from an inventory of polygonal floe shapes saved in a `.mat` file. 
The floe shape inventory `Files/floesInventory.mat` is provided as a limited example and uses a set of shapes found by tracing images of sea ice floes. 
This file can be replaced by an equivalently structured file containing any number of floe shapes (which could be, for example, a more extensive satellite-derived set of floe outlines, or any other arbitrary closed polygons).

These scripts/functions rely on the following MATLAB toolboxes being installed:
* Image Processing Toolbox



## Generation procedure

The floe field generation is divided into two main steps, each of which are independently handled by separate functions.
They are:

1. Generate a vector of floe sizes corresponding to a given FSD (`getFloeSizes`)
2. Randomly distribute the floes within some predefined boundary (`placeFloes`)

The boundary in which the floes are placed can be any arbitrarily-shaped polygon (defined only by a set of vertices).
This allows it to be distinct from, for example, the model domain used in a DEM.
The boundary can be closed (floes are fully contained within the boundary area), or can be doubly-periodic (floes can overlap the boundary and are tiled across opposite boundary sides).
Note that the periodic-boundary option will likely lead to unexpected results for a non-rectangular boundary shape.
Singly-periodic boundaries are not currently suppported, but could be added in the future. 


<img src='/Figures/floeField_circular.png' width='250'>  <img src='/Figures/floeField_periodic.png' width='250'>

### Floe size vector generation (`getFloeSizes`)

#### Summary 

A list of floe sizes is generated following a prescribed FSD, such that the total surface area of the floes satisfies some target value.
(A target SIC can be acheived by setting the target ice area to SIC\*A_d, where A_d is the domain area over which SIC is calculated; e.g., the area enclosed by the boundary). 

Support is included for four different FSD types:

* Power-law (Pareto)
* Lognormal
* Uniform
* Single floe size (Dirac Delta)

#### Procedure

1. Caculate the number of floes for a given continuous FSD, where:
	* N(r) is the number distribtuion of floes *per unit area*, normalized such that ∫πr^2 N(r)dr = SIC
	* (A_d)∫N(r)dr = # of floes where A_d is the total area
2. Numerically revise the estimate (to account for discretization of the FSD)
3. Generate a list of floe sizes by numerically inverting the CCDF
4. Rescale floes to achieve target ice area (to account for discretization of the FSD)

### Floe placement (`placeFloes`)

#### Summary 

Floes are randomly selected from a predefined inventory of floe shapes, then scaled, randomly rotated and placed within the boundary, starting from the largest floe in the floe-size vector and working down the list. 
The random placement location follows a probability distribution that increases with increasing distance away from other previously placed floes (which involves binarizing the floe field on some specified grid and updating the probability matrix after each placement).

After each attempted placement, the floe position is rejected if it overlaps with another floe, or (depending on periodicity) the edge of the boundary. 
To avoid infinite loops, a limit is placed on the number of "rejections" allowed for any floe--if the limit is reached then placement fails. 

Convergence is not guaranteed, and is generally difficult for high concentrations and/or small size-spans. In particular, getting beyond 80-85% SIC can be challenging unless allowing for a very wide size span/small minimum floe size.
Chances for convergence can be increased slightly by increasing the number of "rejection" attempts allowed per floe and (to some extent) by increasing the grid resolution of the probability distribution matrix, both at the expense of computation time.

#### Procedure

Starting from the first entry in the floe size list (the largest floe):

1. Generate a set of candidate floe positions based on the probability matrix <img src='/Figures/probability.png' width='400'>
2. Loop through each candidate floe position
3. Select a random floe from the floe inventory, scale floe to desired floe size and place in candidate position
4. Check that floe doesn't overlap with other floes ("hit test")
	* Check if any points are outside the boundary
	* If using periodic boundary conditions, add "ghost floes" for floes crossing the boundary (tiled across opposite boundary sides)
	* Identify "potential" overlapping floes based on bounding boxes 
	* Loop through potential hits and check for overlapping points
	* <img src='/Figures/hittest_pass.png' width='400'> 
5. If the placement fails, try the next candidate position (return to step 3), otherwise save floe placement and proceed to the next step. If *all* candidate positions fail, floe placement fails and the algorithm exits unsucessfully.
6. Update probability matrix
	* Update the binary mask 
	<img src='/Figures/iceMask.png' width='400'> 
	* get distances to existing floes
	* decrease probability in buffer-zone around each floe
7. Proceed to next floe in the list (step 1)


## Usage and output

To use these functions:
1. Open `makeFloeField.m` in MATLAB.
2. Modify the parameters at the beginning of the script to suit your target domain, desired SIC, and FSD type. Comments in the script provide guidance on how to set these parameters.

The functions output:

* A vector list of floe sizes (1xN)
* A cell array of the x-y coordinates of the floe outlines (1xN)
* A double array of the x-y coordinates of the floe centroids (2xN)
* A logical array of the gridded ice cover binary mask 

These outputs can be converted to an HDF5 file compatible with FloeDyn by using the provided function `makeFloeDynInput`.

The `makeFloeField.m` script visualizes the floe field by converting the cell array of floe outlines to an array of polyshapes, and plotting those.
The script can be easily modified to loop over different input parameters (e.g., different SIC) and generate sets of files


## Performance

No systematic performance/benchmarking tests have completed; however, I have found that placements of up to a few thousand floes in moderate to high SIC can be completed on a laptop within a couple of minutes.


