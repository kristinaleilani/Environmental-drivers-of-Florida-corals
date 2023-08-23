# Environmental-drivers-of-Florida-corals
Scripts, genotype, and environment data to accompany the manuscript "Major environmental drivers of genetic adaptation in Florida corals"


Sequences are deposited on the Sequence Read Archive under Bioproject PRJNA812916. Pipelines for processing raw 2bRAD sequencing data are located here: https://github.com/z0on/2bRAD_denovo. All processing was done on Lonestar 6 at the Texas Advanced Computing Center.
In-depth tutorial on RDA forest can be found here: https://github.com/z0on/RDA-forest




The following R scripts reproduce figures from the manuscript:

**Agaricia_Popgen.R** explores population structure including admixture, PCoA of IBS genetic distance, isolation by distance, and depth partitioning of lineages in Agaricia agaricites.

**Porites_Popgen.R** explores population structure (including all aspects listed above) in Porites astreoides.

**Relatedness.R** explores relatedness structure in both coral species.

**SERC_interpolation.R** performs kriging interpolation to produce rasters for all water quality observations collected by the SERC water quality monitoring network.

**RDA_forest.R** performs gradient forest on genetic distances corrected for population structure (i.e. admixture and spatial autocorrelation) in both species and all lineages.

**Compare_lineages.R** compares adaptation patterns between lineages by measuring distance between turnover grids.





For questions/clarification, please to reach out to Kristina Black (kblack@utexas.edu).
