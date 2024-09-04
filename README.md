# Home range data for "The impacts of anthropogenic linear features on the space-use patterns of two sympatric ungulates"

[https://doi.org/10.5061/dryad.c2fqz61k6](https://doi.org/10.5061/dryad.c2fqz61k6)

## Description of the data and file structure

These data contain 1) the shapefiles of the k-LoCoH home ranges/occurrence distributions, 2) the home ranges' used/available data , and 3) the metadata for each pronghorn and mule deer home range used in this analysis. Because the GPS location data for these animals are sensitive and protected by the Utah Division of Wildlife Resources (UDWR), we have only included the derived home ranges. We included the used/available attributed data as we could not include the full dataset of environmental and linear feature covariates that we used to find these used/available attributes, however all of those data are publicly available.

### Files and variables

#### File: raw\_data.zip

The structure of these data are as follows:

##### `HRs` (home ranges)

These are the shapefiles of the k-LoCoH home ranges. They are projected to the coordinate reference system  WGS 84 / UTM zone 12N (EPSG:32612)

* `hr_id`: the unique ID identifying each home range. Each ID is made up of four parts (XXXXX_SPYYSXXXX_YYYY_SN) 
  * XXXXX: unique numerical ID padded with 0s
  * SPYYSXXXX: unique animal ID;
    * SP: species code (MD for mule deer and PR for pronghorn)
    * YY: last two digits of the year the animal was collared
    * S: sex of the animal
    * XXXX: a unique numerical ID padded with 0s
  * YYYY: year of the home range
  * SN: season of the home range (either "02" for February/winter or "07" for July/summer)
* `level`: the isopleth level of that home range (i.e., how intensely that region of the home range was used; ranges from 15-95%)
* `area`: the area of that isopleth in m2
* `what`: indicates that these are estimates of home ranges (each row is labled "estimate")

##### `used_avail_long`

These are the average (weighted by intensity of use) of the used environmental and linear features within the home range and the average availabilities within a 1km buffer around the centroid of the home range.

* `hr_id`: the unique ID identifying each home range (see above)
* `attribute`: the environmental or linear feature attribute
  * `rd.pvd`: paved roads
  * `rd.unpvd`: unpaved roads
  * `fence`: fence
  * `cvr.shrub`: shrub cover
  * `cvr.tree`: tree cover
  * `elev`: elevation
  * `rough`: Terrain Roughness Index
  * `pdsi`: Palmer Drought Severity Index (PDSI)
  * `forage`: an index created by multiplying NDVI by herbaceous cover
  * `snd`: snow depth
* `used`; `avail`: the average used and available of the given attribute 

##### `hr_info`

Metadata for each home range

* `hr_id`: the unique ID identifying each home range (see above)
* `animal_id`: the unique ID identifying each animal (see above) that the home range belongs to
* `species`: the species (either "Mule Deer" or "Pronghorn") that the home range belongs to
* `sex`: the sex of the individual animal (can be "F" for female or "M" for male)
* `year`: the year the animal used this home range
* `month`: the month the animal used this home range (can be either "2" for February or "7" for July)
* `season`: the season corresponding to the home ranges' month (either Winter or Summer)
* `n_pts`: how many GPS points were used to delineate the home range
* `n_days`: how many days the GPS collar recorded location data for this animal on its home range (minimum 20)

## Code/software

To see how these data were used in the analysis for this paper, see the github repository for this paper: [https://github.com/ronanhart/homerangeALFs](https://github.com/ronanhart/homerangeALFs)

## Access information

Data was derived from GPS location data that were provided by the Utah Division of Wildlife Resources (UDWR)
