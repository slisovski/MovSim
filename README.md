# Animal tracks on the globe <img src="man/figures/globeTracks.png" align="right" />

## Background

The animation shows the epic migration of Bar-tailed godwits (Limosa lapponica) from New Zealand to Alaska via the Yellow Sea during northward migration, and the amazing non-stop flight from Alaska to New Zealand (or in some cases via Australia or the southern Pacific Islands) during southward migration. Thanks to Jesse Conklin and Phil Battley for sharing the tracks and for valuable input and discussion sourrounding the migration of this species. The animation does not only provide the individual mirgation track, but puts them in relation to wind (from ECWMF ERA5 dataset) and snow cover (NOAA IMS Northern Hemisphere Daily Snow Cover dataset). On the globe, you see the wind speed over the ocean (median wind speed across 10m, 750p, 800p and 925p) and simulated particle tracks drifting with the wind. The individials circles change size and color according to the wind support during movements (no movement indicated with grey circles). The green inland color shows the topography (based on NOAA ETOPO2v2). On the right panel, you see the distribution of wind support across all individuals separated by the major fligth bouts (legs).

Tutorial is coming soon.

## Installation

``` r
# To install the latest version from Github:
# install.packages("devtools")
devtools::install_github("tylermorganwall/rayshader")
```




