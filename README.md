NASA API Vignette
================

-   [Package requirements](#package-requirements)
-   [Custom functions](#custom-functions)
    -   [annualExoDiscoveries()](#annualexodiscoveries)
    -   [calculateHZ()](#calculatehz)
    -   [hzFluxCalculator()](#hzfluxcalculator)
    -   [habitableExoFinder()](#habitableexofinder)
-   [Exploratory Data Analysis](#exploratory-data-analysis)
    -   [Annual discoveries](#annual-discoveries)
    -   [Discovery methods](#discovery-methods)
    -   [Metallicity correlations](#metallicity-correlations)
    -   [Mass-radius diagram](#mass-radius-diagram)
    -   [Exoplanet habitability](#exoplanet-habitability)
-   [References](#references)

## Package requirements

To re-create this vignette in R, users are required to install:

-   **[tidyverse](https://www.tidyverse.org/)** - encompasses packages
    such as `dplyr` for subsetting data and `ggplot2` for creating
    layered graphics.  
-   **[httr](https://httr.r-lib.org/)** - supplies the `GET()` function
    to programmatically retrieve content from NASA Exoplanet Archive’s
    [TAP
    service](https://exoplanetarchive.ipac.caltech.edu/docs/TAP/usingTAP.html#examples).
-   **[jsonlite](https://cran.r-project.org/web/packages/jsonlite/vignettes/json-aaquickstart.html)** -
    supplies the `fromJSON()` function to simplify JSON content into an
    atomic vector.
-   **[latex2exp](https://cran.r-project.org/web/packages/latex2exp/vignettes/using-latex2exp.html)** -
    “an R package that parses and converts LaTeX math formulas to R’s
    plotmath expressions.” This is useful for labeling masses and radii
    on plots using [solar system
    symbols](https://solarsystem.nasa.gov/resources/680/solar-system-symbols/).
-   **[ggrepel](https://ggrepel.slowkow.com/)** - eliminate overlapping
    text labels in the `ggplot2` exoplanet mass-radius diagram.

## Custom functions

### annualExoDiscoveries()

This custom function programmatically retrieves data from the [NASA
Exoplanet Archive’s TAP
service](https://exoplanetarchive.ipac.caltech.edu/docs/TAP/usingTAP.html#examples).
Users can specify one of two tables - Planetary Systems (*ps*) or
Planetary Systems Composite Parameters (*pscomppars*) - as well as a
range for the year(s) in which planets were discovered. The default
values for this function are:

-   `tableName = "pscomppars"` - according to the [table
    definitions](https://exoplanetarchive.ipac.caltech.edu/docs/API_PS_columns.html)
    for *ps* and *pscomppars*, “PSCompPars is a more filled-in table,
    with only one row per planet, enabling a more statistical view of
    the known exoplanet population and their host environments. This
    table provides a more complete, though not necessarily
    self-consistent, set of parameters.”
-   `startYear = 1989` - the earliest listing in the PSCompPars table,
    attributed to the planet [HD 114762
    b](https://exoplanetarchive.ipac.caltech.edu/overview/HD%20114762%20b#planet_data_HD-114762-b).  
-   `endYear = as.integer(format(Sys.Date(), "%Y"))` - the current
    calendar year as understood by the user’s computer.
-   `controversial = 0` - exclude planets for which the confirmation
    status “has been questioned in the published literature.”
-   `cb_flag = 0` - exclude (0) planets which orbits a binary system. An
    option of `1` lists only planets which orbit two or more stars.
-   `format = json` - return queries in JSON format. You may also
    request data as a comma-separated-value file (CSV).

``` r
# Retrieve names, discovery years, discovery methods and various other
# planetary and stellar parameters.
# The default search begins in 1989 (the earliest date in the pscomppars table)
# and ends in the current calendar year: format(Sys.Date(), "%Y").
annualExoDiscoveries <- function(tableName = "pscomppars", 
                                 startYear = 1989, 
                                 endYear = as.integer(format(Sys.Date(), "%Y")), 
                                 controversialFlag = 0,
                                 cb_flag = 0,
                                 format = "json"){
  # Create URL string
  urlString <- paste0("https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query=select+pl_name,disc_year,discoverymethod,pl_orbper,pl_rade,pl_bmasse,pl_radj,pl_bmassj,pl_eqt,pl_dens,st_spectype,st_teff,st_lum,pl_controv_flag,pl_orbeccen,pl_orbsmax,st_mass,st_metratio,st_met+from+", 
                      tableName, "+where+disc_year+between+", 
                      startYear, "+and+", endYear, 
                      "+and+pl_controv_flag+=+", controversialFlag, 
                      "+and+cb_flag+=+", cb_flag,
                      "&format=", format)
  # Provide string to httr GET function
  apiCall <- GET(urlString)
  
  if(format == "json"){
    # Convert JSON content to data frame, rename columns
    apiContent <- apiCall$content %>% rawToChar() %>% fromJSON() %>%
      mutate(luminosityRatio = 10^(st_lum)) 
  } else {
    # Specify format as CSV, convert to data frame
    apiContent <- as.data.frame(read_csv(urlString))
  }
  # Return formatted data frame
  return(apiContent)
}
```

### calculateHZ()

This function calculates planetary habitable zones and their associated
stellar flux boundaries for a host star with effective temperature
`tempEff` and a stellar luminosity `luminosityRatio`. The calculations
are based on formulae defined by Kopparapu et al., whereby the effective
solar flux *S*<sub>*e**f**f*</sub> equates to
*S*<sub>*e**f**f*</sub> = *S*<sub>*e**f**f*⊙</sub> + (*a* ⋅ *T*<sub>⋆</sub>) + (*b* ⋅ *T*<sub>⋆</sub><sup>2</sup>) + (*c* ⋅ *T*<sub>⋆</sub><sup>3</sup>) + (*d* ⋅ *T*<sub>⋆</sub><sup>4</sup>)
and the corresponding habitability zone distances, *d*, equate to
*d* = (*L*/*L* ⊙ )/(*S*<sub>*e**f**f*</sub>)<sup>0.5</sup> AU (Kopparapu
et al., 2014).

The required parameters for this function are:

-   `tempEff` - the effective temperature of a host star such that
    <sub>*e**f**f*</sub> − 5780 = *T*<sub>⋆</sub>.
-   `luminosityRatio` - stellar luminosity, *L*/*L*⊙, required to
    calculate habitability zone distances, *d*. These values may be
    calculated with the `annualExoDiscoveries()` function, which finds
    the inverse logarithm of the stellar luminosity in the *PSCompPars*
    table (`st_lum`).

The output of this function is a list with four numeric parameters -
*optimisticInnerDist*, *optimisticOuterDist*, *optimisticInnerFlux*, and
*optimisticOuterFlux* - which may be used in the function
`hzFluxCalculator()`.

``` r
# Calculate habitable stellar flux boundaries for exoplanetary habitable zones. 
# Distances are returned in Astronomical Units (AU).
# Formula for s_eff and its coefficients is provided by Kopparapu et al.
# https://iopscience.iop.org/article/10.1088/2041-8205/787/2/L29
# Re-factored to R from John Armstrong's Python code at
# https://depts.washington.edu/naivpl/sites/default/files/hzcalc.py.txt
calculateHZ <- function(tempEff, luminosityRatio){
  
  # Initiate empty vectors
  s_eff <- vector()
  distanceFromStar <- vector()
  
  starTemp <- vector()
  recentVenus <- vector()
  runawayGreenhouse <- vector()
  maxGreenhouse <- vector()
  earlyMars <- vector()
  fivemeRunaway <- vector()
  tenthmeRunaway <- vector()
  
  # Populate variables with coefficients from research paper by Kopparapu et al. 
  s_eff_sun  = c(1.776, 1.107, 0.356, 0.320, 1.188, 0.99)
  a <- c(2.136e-4, 1.332e-4, 6.171e-5, 5.547e-5, 1.433e-4, 1.209e-4)
  b <- c(2.533e-8, 1.580e-8, 1.698e-9, 1.526e-9, 1.707e-8, 1.404e-8)
  c <- c(-1.332e-11, -8.308e-12, -3.198e-12, -2.874e-12, -8.968e-12, -7.418e-12)
  d <- c(-3.097e-15, -1.931e-15, -5.575e-16, -5.011e-16, -2.084e-15, -1.713e-15)
  
  t_star <- tempEff-5780
  
  
  for (i in 1:length(a)){
    # Calculate effective solar flux (s_eff) using formula
    # from research paper by Kopparapu et al.
    s_eff[i] <- s_eff_sun[i] + 
      a[i]*t_star + b[i]*t_star^2 + c[i]*t_star^3 + d[i]*t_star^4
    
    # Calculate corresponding inner/outer habitability zone distances
    distanceFromStar[i] <- (luminosityRatio/s_eff[i])^0.5
    
    optimisticInnerDist <- distanceFromStar[1]
    optimisticOuterDist <- distanceFromStar[4]
    
    # Calculate effective solar flux incident on the planet
    optimisticInnerFlux <- s_eff[1]
    optimisticOuterFlux <- s_eff[4]
    
    return(list(optimisticInnerDist = optimisticInnerDist, 
                optimisticOuterDist = optimisticOuterDist, 
                optimisticInnerFlux = optimisticInnerFlux, 
                optimisticOuterFlux = optimisticOuterFlux))
  }
}
```

### hzFluxCalculator()

This custom function calculates the minima and the maxima for a planet’s
habitability zone (in units of *AU*) and stellar flux (in units of
*dex*). It requires the name of a data set and operates with the
following default parameters:

-   `earthMassCol = "pl_bmasse"` - a vector with planetary masses in
    units of Earth mass (*M*⊕).
-   `starSpecTypeCol = "st_spectype"` - a vector listing the spectral
    type of host stars.
-   `effectiveTempCol = "st_teff"` - a vector listing the effective
    temperatures of host stars.
-   `luminosityRatioCol = "luminosityRatio"` - a vector with the stellar
    luminosity ratios, *L*/*L*⊙.

These default column names are based on the variables in the [Planetary
Systems Composite Parameters
(PSCompPars)](https://exoplanetarchive.ipac.caltech.edu/docs/API_PS_columns.html)
table. The *spectralClass* column uses the `substr()` function to
distill the complete classification of a star into its main class (O, B,
A, etc.).

Habitability zone distances and incident flux are calculated only for
planets with a mass of 10*M*⊕ or less. This is the hypothetical
threshold for planets with an appreciable composition of “volatiles”
such as *H*<sub>2</sub>*O* and *N**H*<sub>3</sub> (Kuchner, 2003).

``` r
# Customer function to calculate values for 
# inner and outer habitable zone, flux
hzFluxCalculator <- function(data, earthMassCol = "pl_bmasse", 
                             starSpecTypeCol = "st_spectype", 
                             effectiveTempCol = "st_teff",
                             luminosityRatioCol = "luminosityRatio"){
  
  data %>% mutate(innerHZ = NA, outerHZ = NA, innerFlux = NA, 
                  outerFlux = NA, spectralClass = NA)
  
  earthMassCol <- data[ , earthMassCol]
  starSpecTypeCol <- data[ , starSpecTypeCol]
  effectiveTempCol <- data[ , effectiveTempCol]
  luminosityRatioCol <- data[ , luminosityRatioCol]
  
  for(i in 1:length(earthMassCol)){
    if(!is.na(starSpecTypeCol[i])){
      data$spectralClass[i] <- substr(starSpecTypeCol[i], 1, 1)
    } else {
      data$spectralClass[i] <- NA
    }
    
    
    if(!is.na(earthMassCol[i]) & earthMassCol[i] <= 10 & 
       earthMassCol[i] >= 0.1){
      
      
      hzVars <- calculateHZ(effectiveTempCol[i], 
                            luminosityRatioCol[i])
      
      data$innerHZ[i] <- hzVars[[1]]
      data$outerHZ[i] <- hzVars[[2]]
      
      data$innerFlux[i] <- hzVars[[3]]
      data$outerFlux[i] <- hzVars[[4]]
      
    } else {
      data$innerHZ[i] <- NA
      data$outerHZ[i] <- NA
      data$innerFlux[i] <- NA
      data$outerFlux[i] <- NA
    }
  }
  
  return(data)
}
```

### habitableExoFinder()

This function produces a data frame of potentially habitable exoplanets
from a set of general habitability criteria, many of which are inspired
by the University of Puerto Rico’s [Planetary Habitability
Laboratory](http://phl.upr.edu/projects/habitable-exoplanets-catalog/methods).

These parameters can be tuned using arguments from published research.
For simplicity, the default values are:

-   `minEarthMass = 0.1`
-   `maxEarthMass = 5`
-   `minEarthRadius = 0.5`
-   `maxEarthRadius = 1.5`
-   `maxInnerFlux = 1.5`
-   `maxOuterFlux = 0.20`
-   `minTemp = 273` - units of Kelvin
-   `maxTemp = 340` - units of Kelvin

``` r
# Create function to identify potentially habitable exoplanets. 
# Default function parameters provided by Planetary Habitability Laboratory, 
# http://phl.upr.edu/projects/habitable-exoplanets-catalog
habitableExoFinder <- function(data, minEarthMass = 0.1, maxEarthMass = 5, 
                               minEarthRadius = 0.5, maxEarthRadius = 1.5,
                               maxInnerFlux = 1.5, maxOuterFlux = 0.20,
                               minTemp = 273, maxTemp = 340){
  
  # Subset data using provided parameters
  habitablePlanets <- data %>% select(pl_name, pl_eqt, spectralClass, 
                                      pl_bmasse, pl_rade, pl_orbeccen, 
                                      pl_orbsmax, innerHZ, outerHZ, 
                                      innerFlux, outerFlux) %>% 
    filter(spectralClass %in% c("F", "G", "K", "M") & 
             (pl_orbsmax >= innerHZ) & (pl_orbsmax <= outerHZ) & 
             (pl_bmasse >= minEarthMass) &
             (pl_bmasse <= maxEarthMass) & 
             (pl_rade >= minEarthRadius) & (pl_rade <= maxEarthRadius) &
             (innerFlux <= maxInnerFlux) & (outerFlux >= maxOuterFlux) &
             (pl_eqt <= maxTemp | (is.na(pl_eqt))) & 
             (pl_eqt >= minTemp | (is.na(pl_eqt == NA))))
  
  return(habitablePlanets)
}
```

## Exploratory Data Analysis

### Annual discoveries

The `annualExoDiscoveries()` function retrieves the latest exoplanet
data from NASA’s Exoplanet Archive.

``` r
# Retrieve latest exoplanet data
exoplanetData <- annualExoDiscoveries()
exoplanetData
```

    ## # A tibble: 4,462 × 20
    ##    pl_name    disc_year discoverymethod
    ##    <chr>          <int> <chr>          
    ##  1 OGLE-2016…      2020 Microlensing   
    ##  2 GJ 480 b        2020 Radial Velocity
    ##  3 Kepler-27…      2013 Transit        
    ##  4 Kepler-82…      2016 Transit        
    ##  5 K2-283 b        2018 Transit        
    ##  6 Kepler-47…      2016 Transit        
    ##  7 HAT-P-15 b      2010 Transit        
    ##  8 HD 149143…      2005 Radial Velocity
    ##  9 HD 210702…      2007 Radial Velocity
    ## 10 HIP 12961…      2010 Radial Velocity
    ## # … with 4,452 more rows, and 17 more
    ## #   variables: pl_orbper <dbl>,
    ## #   pl_rade <dbl>, pl_bmasse <dbl>,
    ## #   pl_radj <dbl>, pl_bmassj <dbl>,
    ## #   pl_eqt <dbl>, pl_dens <dbl>,
    ## #   st_spectype <chr>, st_teff <dbl>,
    ## #   st_lum <dbl>, …

``` r
# Print a subset of observations
exoplanetData
```

    ## # A tibble: 4,462 × 20
    ##    pl_name    disc_year discoverymethod
    ##    <chr>          <int> <chr>          
    ##  1 OGLE-2016…      2020 Microlensing   
    ##  2 GJ 480 b        2020 Radial Velocity
    ##  3 Kepler-27…      2013 Transit        
    ##  4 Kepler-82…      2016 Transit        
    ##  5 K2-283 b        2018 Transit        
    ##  6 Kepler-47…      2016 Transit        
    ##  7 HAT-P-15 b      2010 Transit        
    ##  8 HD 149143…      2005 Radial Velocity
    ##  9 HD 210702…      2007 Radial Velocity
    ## 10 HIP 12961…      2010 Radial Velocity
    ## # … with 4,452 more rows, and 17 more
    ## #   variables: pl_orbper <dbl>,
    ## #   pl_rade <dbl>, pl_bmasse <dbl>,
    ## #   pl_radj <dbl>, pl_bmassj <dbl>,
    ## #   pl_eqt <dbl>, pl_dens <dbl>,
    ## #   st_spectype <chr>, st_teff <dbl>,
    ## #   st_lum <dbl>, …

As of Tue Oct 5 21:15:51 2021, the NASA Exoplanet Archive’s [Planetary
Systems Composite
Parameters](https://exoplanetarchive.ipac.caltech.edu/docs/API_PS_columns.html)
(PSCompPars) table lists 4501 confirmed exoplanet observations. The
stacked bar plot below enumerates the annual number of exoplanet
findings since 1989 and highlights 2014 and 2016 as the most prolific
years for discovery.

``` r
# Retrieve latest exoplanet data
annualDiscoveries <- exoplanetData 

# Create horizontal bar plot with annual number of discoveries
annualDiscoveryBar <- ggplot(annualDiscoveries, aes(x = as.character(disc_year)))
annualDiscoveryBar + geom_bar(aes(fill = discoverymethod), 
                             alpha = 0.8, position = "stack") +
  labs(x = "Discovery year", y = "Count",
       title = "Exoplanet discoveries over time", 
       subtitle = "Grouped by discovery method") +
  theme(axis.text.x = element_text(angle = 45)) + 
  geom_text(stat="count", aes(label=..count..), 
            vjust=0.5, hjust = -0.01, size = 2.5) +
  theme(legend.position = c(0.65, 0.33)) +
  # Split legend into two columns
  guides(fill=guide_legend(ncol=2)) +
  # Change color palette for bars based on
  # http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
  scale_fill_brewer(palette = "Spectral", name = "Discovery method") +
  coord_flip() 
```

![](README_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

The contingency table below summarizes the cumulative number of
observations for each discovery method.

``` r
# Contingency table 
# Total number of exoplanets found with each discovery method
discoveriesByMethod <- table(annualDiscoveries$discoverymethod)

# Display as data frame for aesthetics
knitr::kable(discoveriesByMethod,
             col.names = c("Discovery method", "Frequency"))
```

| Discovery method              | Frequency |
|:------------------------------|----------:|
| Astrometry                    |         1 |
| Disk Kinematics               |         1 |
| Imaging                       |        45 |
| Microlensing                  |       113 |
| Orbital Brightness Modulation |         4 |
| Pulsar Timing                 |         6 |
| Pulsation Timing Variations   |         2 |
| Radial Velocity               |       870 |
| Transit                       |      3398 |
| Transit Timing Variations     |        22 |

Of the known 4462 exoplanets, 75.5% were observed while transiting their
host star and temporarily reducing its brightness. Another 19.3% were
observed indirectly via the radial velocity method, whereby the planet
and its star orbit around a common center of gravity and prompt
noticeable Doppler shifts in the stellar spectrum.

``` r
circumbinaryPlanets <- annualExoDiscoveries(cb_flag = 1) 
circumbinaryPlanets
```

    ## # A tibble: 39 × 20
    ##    pl_name    disc_year discoverymethod
    ##    <chr>          <int> <chr>          
    ##  1 NSVS 1425…      2019 Eclipse Timing…
    ##  2 RR Cae b        2012 Eclipse Timing…
    ##  3 2MASS J19…      2015 Eclipse Timing…
    ##  4 SR 12 AB c      2010 Imaging        
    ##  5 MXB 1658-…      2017 Eclipse Timing…
    ##  6 Kepler-16…      2016 Transit        
    ##  7 VHS J1256…      2015 Imaging        
    ##  8 2MASS J01…      2013 Imaging        
    ##  9 Kepler-45…      2015 Transit        
    ## 10 Kepler-35…      2011 Transit        
    ## # … with 29 more rows, and 17 more
    ## #   variables: pl_orbper <dbl>,
    ## #   pl_rade <dbl>, pl_bmasse <dbl>,
    ## #   pl_radj <dbl>, pl_bmassj <dbl>,
    ## #   pl_eqt <dbl>, pl_dens <dbl>,
    ## #   st_spectype <chr>, st_teff <dbl>,
    ## #   st_lum <dbl>, …

The table above lists
`r`length(binaryStarPlanets$pl\_name)`planets which are thought to orbit two or more stars - the so-called "circumbinary planets" [@kepler-16]. They are presented in the table above for completeness but, to reduce the complexity of two-body problems, they are otherwise excluded by default from the`annualExoDiscoveries()\`
function.

### Discovery methods

Each observation method excels in specific scenarios. The transit and
radial velocity detection methods favor planets which orbit their star
at an average distance of 0.12-1.6 AU.

``` r
# Subset data to include only detection methods with a relatively large number
# of exoplanet discoveries
extendedDiscoveryProp <- exoplanetData
extendedDiscoveryProp <- extendedDiscoveryProp %>% 
  filter(discoverymethod %in% c("Transit", "Radial Velocity", "Microlensing", "Imaging") &
           !is.na(pl_orbsmax))

# Display mean and median semi-major axes (in units of AU) for 
# planets observed by these methods
discoverySummaries <- extendedDiscoveryProp %>% group_by(discoverymethod) %>%
  summarise(meanSMA = mean(pl_orbsmax), medianSMA = median(pl_orbsmax))
knitr::kable(discoverySummaries, 
             col.names = c("Discovery method", "Mean SMA (AU)",
                           "Median SMA (AU)"))
```

| Discovery method | Mean SMA (AU) | Median SMA (AU) |
|:-----------------|--------------:|----------------:|
| Imaging          |   607.8113636 |        181.0000 |
| Microlensing     |     2.6661518 |          2.2500 |
| Radial Velocity  |     1.6044780 |          1.0200 |
| Transit          |     0.1261318 |          0.0788 |

Direct imaging, on the other hand, requires planets to be relatively far
from a star in order for the stellar brightness not to overwhelm the
planetary dimness. In this data set, the median distance for a planet
that was directly observed is 162 AU. The boxplot below shows the
distribution of semi-major axes for four of the most productive
exoplanet observation methods.

``` r
orbsmaxBoxPlot <- ggplot(extendedDiscoveryProp, 
                         aes(x = discoverymethod, y = pl_orbsmax,
                             fill = discoverymethod))
orbsmaxBoxPlot + geom_boxplot() +
  # Remove legend after coloration
  theme(legend.position = "none") +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
   labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(x = "Discovery method", y = "Orbit semi-major axis (log)",
       title = "Distribution of semi-major axes for various discovery methods") +
  # Add log ticks to left side
  annotation_logticks(sides="l")
```

![](README_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

Direct imaging also favors young stars, which tend to be “self-luminous
due to ongoing contraction and…accretion” (service), 2016). The
combination of large semi-major axes and a luminous nature is generally
attributed to giants that match or exceed the mass of Jupiter. This is
corroborated by the scatter plot below, whereby the bulk of
directly-imaged planets have a mass of approximately 10*M*<sub>*J*</sub>
and reside at a distance of 10-10,000 astronomical units from their
star.

``` r
# New vector with temporary data
orbsmaxMassData <- exoplanetData 

# Scatter plot of masses/radii for discovered exoplanets
# Use LaTeX to denote the standard astronomical symbol for the Earth
orbsmaxMassScatter <- ggplot(extendedDiscoveryProp, aes(x = pl_orbsmax, y = pl_bmassj))
orbsmaxMassScatter + geom_point(aes(color = pl_orbeccen, shape = discoverymethod), 
                                alpha = 0.6, position = "jitter") +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
   labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
   labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  # Select color palette for eccentricity from 
  # https://rdrr.io/r/grDevices/palettes.html
  scale_color_gradientn(colours = terrain.colors(5)) +
  labs(x = "Semi-major axis [AU]", y = TeX(r'(Planet mass $(M_{J})$)'),
       title = "Planetary mass versus semi-major axis", col = "Orbit eccentricity",
       shape = "Discovery method") 
```

    ## Warning: Removed 17 rows containing missing
    ## values (geom_point).

![](README_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

### Metallicity correlations

We can observe at least one more obvious trend from the data - that of
the “giant planet–metallicity correlation” - whereby giant planets “tend
to appear around metal-rich stars” (Adibekyan, 2019). COonversely,
“there is little or no dependence on metallicity for low-mass planets
such as super-Earths” (Hasegawa & Pudritz, 2014). If we use 10*M*⊕ as
the threshold for super-Earths and create a histogram for the
metallicity distribution of stars in the exoplanet database, we find
that lanets with masses of 10*M*⊕ or less are centered around stars with
\[Fe/H\] = 0. On the other hand, giant plants (those which exceed
10*M*⊕) exhibit a slightly right-skewed distribution with a substantial
concentration near \[Fe/H\] = 0.13.

``` r
metallicityData <- extendedDiscoveryProp %>% filter(st_metratio == "[Fe/H]" &
                                                      !is.na(pl_bmasse))
metallicityData %>% mutate(giantPlFlag = NA)
```

    ## # A tibble: 3,459 × 21
    ##    pl_name      disc_year discoverymethod
    ##    <chr>            <int> <chr>          
    ##  1 Kepler-276 c      2013 Transit        
    ##  2 Kepler-829 b      2016 Transit        
    ##  3 K2-283 b          2018 Transit        
    ##  4 Kepler-477 b      2016 Transit        
    ##  5 HAT-P-15 b        2010 Transit        
    ##  6 HD 149143 b       2005 Radial Velocity
    ##  7 HD 210702 b       2007 Radial Velocity
    ##  8 HIP 12961 b       2010 Radial Velocity
    ##  9 XO-5 b            2008 Transit        
    ## 10 HD 5608 b         2012 Radial Velocity
    ## # … with 3,449 more rows, and 18 more
    ## #   variables: pl_orbper <dbl>,
    ## #   pl_rade <dbl>, pl_bmasse <dbl>,
    ## #   pl_radj <dbl>, pl_bmassj <dbl>,
    ## #   pl_eqt <dbl>, pl_dens <dbl>,
    ## #   st_spectype <chr>, st_teff <dbl>,
    ## #   st_lum <dbl>, …

``` r
for (i in 1:length(metallicityData$pl_name)){
  if (metallicityData$pl_bmasse[i] >= 10){
    metallicityData$giantPlFlag[i] = "Giant"
  } else if(metallicityData$pl_bmasse[i] < 10){
    metallicityData$giantPlFlag[i] = "Sub-giant"
  }
  else {
    metallicityData$giantPlFlag[i] = NA
  }
}

metallicityHisto <- ggplot(metallicityData, aes(x = st_met))
metallicityHisto + geom_histogram(aes(y = ..density.., 
                                      fill = giantPlFlag),
                                  bins = 50, color = "red") +
  labs(x = "Stellar metallicity [Fe/H]",
       title = "The distribution of metallicity in giant and sub-giant planets",
       fill = "") + 
  geom_density(adjust = 0.5, alpha = 0.5, aes(fill = giantPlFlag))
```

![](README_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

An empirical cumulative distribution function affirms that, while 50% of
sub-giant (*M* &lt; 10*M*⊕) planets orbit a star with a metallicity
\[Fe/H\] = 0, 50% of giant planets (*M* &gt;  = 10*M*⊕) orbit stars with
a metallicity of \[Fe/H\] = 0.06.

``` r
# Render and label metallicity ECDF 
metallicityHisto + stat_ecdf(geom = "step", aes(color = giantPlFlag)) +
  labs(title="Empirical Cumulative Density Function \n Stellar Metallicity",
     y = "ECDF", x="[Fe/H]", color = "Planet category")
```

![](README_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
# Group data by planet status (giant/sub-giant) and calculate 
# mean/median metallicity. 
metallicityAverages <- metallicityData %>% group_by(giantPlFlag) %>%
  summarise(meanMetallicity = mean(st_met), medianMetallicity = median(st_met))
knitr::kable(metallicityAverages, 
             col.names = c("Classification", "Mean stellar metallicity [dex]",
                           "Median stellar metallicity [dex]"))
```

| Classification | Mean stellar metallicity \[dex\] | Median stellar metallicity \[dex\] |
|:---------------|---------------------------------:|-----------------------------------:|
| Giant          |                        0.0506846 |                               0.06 |
| Sub-giant      |                       -0.0131614 |                               0.00 |

### Mass-radius diagram

Planets with radii in the range 10 − 15*R*⊕ comprise 20% of our data
set. The bulk of the remaining observations - more than 60% - consists
of planets in the range of 0.1 − 5*R*⊕.

``` r
radiiFreq <- ggplot(annualDiscoveries, aes(x = pl_rade)) 
radiiFreq + geom_histogram(color = "#123456", fill = "#f7a22b", 
                           aes(y = (..count..)/sum(..count..)),
                           binwidth = 2) +
  labs(title="Distribution of planetary radii") +
  labs(x = TeX(r'(Radius $(R\oplus$))'), 
       y = "Frequency") +
  geom_density()
```

    ## Warning: Removed 7 rows containing non-finite
    ## values (stat_bin).

    ## Warning: Removed 7 rows containing non-finite
    ## values (stat_density).

![](README_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

By combining radii with the masses of planets, we can produce a
mass-radius diagram and calculate planetary densities. From this
diagram, it is also apparent that planetary radii tend to increase with
mass until approximately 1000*M*⊕, at which point gravity impels even
the hardiest planetary material to compress and constrains further
radial growth (Hoolst et al., 2019).

``` r
# New vector with temporary data
tempMassData <- exoplanetData 

# Scatter plot of masses/radii for discovered exoplanets
# Use LaTeX to denote the standard astronomical symbol for the Earth
tempMassScatter <- ggplot(tempMassData, aes(x = pl_bmasse, y = pl_rade))
tempMassScatter + geom_point(aes(col = pl_eqt, size = pl_dens), alpha = 0.6, position = "jitter") +
  scale_color_gradientn(colours = heat.colors(5)) +
  labs(x = TeX(r'(Planet mass $(log(M\oplus))]$))'), y = TeX(r'(Planet radius $(R\oplus)$)'),
       title = "Planetary mass-radius diagram", col = "Equilbrium temperature (K)",
       size = TeX(r'(Planet density $(g/cm^3)$)')) +
  # Label planets exceeding R = 25 (Earth), eliminate overlapping labels
  # using geom_text_repel from the ggrepel package
  geom_text_repel(aes(label = ifelse(pl_rade >= 25, pl_name,''))) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
   labels = scales::trans_format("log10", scales::math_format(10^.x)))
```

    ## Warning: Removed 94 rows containing missing
    ## values (geom_point).

    ## Warning: Removed 24 rows containing missing
    ## values (geom_text_repel).

![](README_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

### Exoplanet habitability

Among the 4501 exoplanets in the NASA Exoplanet Archive, which have the
potential to harbor life? The [Planetary Habitability
Laboratory](http://phl.upr.edu/projects/habitable-exoplanets-catalog/methods)
(PHL) attempts to answer this by narrowing conditions such that

1.  “The planet orbits an F, G, K, or M star.”
2.  “The planet orbits within the optimistic habitable zone defined by
    Kopparapu et al. (2014).”
3.  “The planet has a radius between 0.5 to 2.5 Earth radii or a minimum
    mass between 0.1 to 10 Earth masses.”

Moreover, as defined by Zsom et al., “an exoplanet is habitable if
liquid water, a crucial ingredient for life as we know it, is present on
its surface, and if the surface temperature and pressure are such that
complex organic molecules are stable” (Zsom et al., 2013). This this
end, we can use criteria from the PHL and planetary equilibrium
temperatures to compile our own list of potentially habitable
exoplanets. We can then compare this list against one compiled by the
PHL [Habitable Planets
Catalog](http://phl.upr.edu/projects/habitable-exoplanets-catalog).

First, we calculate the maxima and minima for habitable zone distances
and solar flux using formulae provided by Kopparapu et al. (Kopparapu et
al., 2014). The `hzFluxCalculator()` function takes an existing data
frame (such as the one produced by the `annualExoDiscoveries()`
function) and appends to it the columns *innerHZ*, *outerHZ*,
*innerFlux*, and *outerFlux*. Individual habitable zone values are
calculated using the `calculateHZ()` function and extended to larger
data sets via the `hzFluxCalculator()` function.

``` r
# Exoplanet data from the NASA Exoplanet Archive
planetData <- exoplanetData

# Calculate maxima/minima for habitabilize zone distances 
# and the solar flux incident on each exoplanet
planetData <- hzFluxCalculator(planetData)

head(planetData, n = 10)
```

    ## # A tibble: 10 × 25
    ##    pl_name    disc_year discoverymethod
    ##    <chr>          <int> <chr>          
    ##  1 OGLE-2016…      2020 Microlensing   
    ##  2 GJ 480 b        2020 Radial Velocity
    ##  3 Kepler-27…      2013 Transit        
    ##  4 Kepler-82…      2016 Transit        
    ##  5 K2-283 b        2018 Transit        
    ##  6 Kepler-47…      2016 Transit        
    ##  7 HAT-P-15 b      2010 Transit        
    ##  8 HD 149143…      2005 Radial Velocity
    ##  9 HD 210702…      2007 Radial Velocity
    ## 10 HIP 12961…      2010 Radial Velocity
    ## # … with 22 more variables:
    ## #   pl_orbper <dbl>, pl_rade <dbl>,
    ## #   pl_bmasse <dbl>, pl_radj <dbl>,
    ## #   pl_bmassj <dbl>, pl_eqt <dbl>,
    ## #   pl_dens <dbl>, st_spectype <chr>,
    ## #   st_teff <dbl>, st_lum <dbl>,
    ## #   pl_controv_flag <int>, …

Next, we supply this data frame to the `habitableExoFinder()` function
alongside our criteria for habitability. In the code chunk below, we use
effective temperatures in the range of 181 − 279 Kelvin and planet radii
in the range of 2.5 − 10*R*⊕.

``` r
# List habitable planets using "optimistic" parameters from
# Planetary Habitability Laboratory
listHabitablePlanets <- habitableExoFinder(planetData, minTemp = 181, maxTemp = 279,
                                           maxEarthMass = 10, maxEarthRadius = 2.5)
knitr::kable(listHabitablePlanets)
```

| pl\_name | pl\_eqt | spectralClass | pl\_bmasse | pl\_rade | pl\_orbeccen | pl\_orbsmax | innerHZ | outerHZ | innerFlux | outerFlux |
|:---------|--------:|:--------------|-----------:|---------:|-------------:|------------:|--------:|--------:|----------:|----------:|

Our combination of habitable zone distances, incident flux, effective
temperatures, planet masses, and planet radii yield 0 potentially
habitable exoplanets.

## References

<div id="refs" class="references csl-bib-body hanging-indent"
line-spacing="2">

<div id="ref-geosciences9030105" class="csl-entry">

Adibekyan, V. (2019). Heavy metal rules. I. Exoplanet incidence and
metallicity. *Geosciences*, *9*(3).
<https://doi.org/10.3390/geosciences9030105>

</div>

<div id="ref-pub.1041145028" class="csl-entry">

Hasegawa, Y., & Pudritz, R. E. (2014). PLANET TRAPS AND PLANETARY CORES:
ORIGINS OF THE PLANET-METALLICITY CORRELATION. *The Astrophysical
Journal*, *794*(1), 25. <https://doi.org/10.1088/0004-637x/794/1/25>

</div>

<div id="ref-doi:10.1080/23746149.2019.1630316" class="csl-entry">

Hoolst, T. V., Noack, L., & Rivoldini, A. (2019). Exoplanet interiors
and habitability. *Advances in Physics: X*, *4*(1), 1630316.
<https://doi.org/10.1080/23746149.2019.1630316>

</div>

<div id="ref-Kopparapu_2014" class="csl-entry">

Kopparapu, R. K., Ramirez, R. M., SchottelKotte, J., Kasting, J. F.,
Domagal-Goldman, S., & Eymet, V. (2014). HABITABLE ZONES AROUND
MAIN-SEQUENCE STARS: DEPENDENCE ON PLANETARY MASS. *The Astrophysical
Journal*, *787*(2), L29. <https://doi.org/10.1088/2041-8205/787/2/l29>

</div>

<div id="ref-pub.1058671435" class="csl-entry">

Kuchner, M. J. (2003). Volatile-rich earth-mass planets in the habitable
zone. *The Astrophysical Journal*, *596*(1), l105–l108.
<https://doi.org/10.1086/378397>

</div>

<div id="ref-cite-key" class="csl-entry">

service), S. (Online. (2016). *Methods of detecting exoplanets 1st
advanced school on exoplanetary science* (V. Bozza, L. Mancini, & A.
Sozzetti, Eds.). Cham : Springer International Publishing : Imprint:
Springer, 2016. <https://catalog.lib.ncsu.edu/catalog/NCSU3603337>

</div>

<div id="ref-pub.1021397281" class="csl-entry">

Zsom, A., Seager, S., Wit, J. de, & Stamenković, V. (2013). TOWARD THE
MINIMUM INNER EDGE DISTANCE OF THE HABITABLE ZONE. *The Astrophysical
Journal*, *778*(2), 109. <https://doi.org/10.1088/0004-637x/778/2/109>

</div>

</div>
