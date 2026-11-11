# Plot STRUCTURE like results on a map.

Plots the mean cluster assignment for each population on a map using the
`scatterpie` package alongside any additional simple feature objects
(see `sf` from the `sf` package). Assignments must be given in the
format provided by
[`plot_structure`](https://hemstrow.github.io/snpR/reference/plot_structure.md).
This function is a wrapper which sources code from github.

## Usage

``` r
plot_structure_map(
  assignments,
  k,
  facet,
  pop_coordinates,
  layers = NULL,
  pop_names = T,
  viridis.option = "viridis",
  alt.palette = NULL,
  radius_scale = 0.05,
  label_args = NULL,
  crop = FALSE,
  scale_bar = list(),
  compass = list(style = ggspatial::north_arrow_fancy_orienteering, location = "br"),
  ask = TRUE
)
```

## Arguments

- assignments:

  Structure like results, parsed in or generated via
  [`plot_structure`](https://hemstrow.github.io/snpR/reference/plot_structure.md),
  which generates the needed plot data.

- k:

  numeric. Value of K (number of clusters) to plot.

- facet:

  character. The facet by which data is broken down in the passed
  assignments.

- pop_coordinates:

  sf object, see the documentation for `sf` function from the `sf`
  package. sf object containing points/coordinates for each facet level.
  Must contain a column of data with population labels named identically
  to the provided facet (for example, named "pop" if "pop" is the
  provided facet.)

- layers:

  list of `ggplot2` layer objects, default NULL. Additional ggplot
  layers to be plotted in order, below the pie charts, such as maps with
  borders, temperatures, forest cover, etc. As a special note,
  `link[ggnewscale]{new_scale_fill}` can be used to add additional
  `fill` aesthetic layers without conflict from the resulting pie
  charts. See examples.

- pop_names:

  logical, default T. If true, facet level names will be displayed on
  the map.

- viridis.option:

  character, default "viridis". Viridis color scale option. See
  `scale_gradient` for details.

- alt.palette:

  character or NULL, default NULL. Optional palette of colors to use
  instead of viridis palette the pie charts.

- radius_scale:

  numeric 0-1, default 0.05. Scale for pie chart radii as a proportion
  of the total map space.

- label_args:

  list, default NULL. Named list of arguments passed to
  [`geom_label_repel`](https://ggrepel.slowkow.com/reference/geom_text_repel.html).
  For example, passing list(max.overlaps = 14) will add the max.overlaps
  argument to the function call.

- crop:

  logical, default F. If TRUE, will will crop the plot around the sample
  points. If false will show the full extent of the data, often set by
  any additional sf objects being plotted.

- scale_bar:

  list or NULL, default [`list()`](https://rdrr.io/r/base/list.html).
  Arguments passed to the
  [`annotation_scale`](https://paleolimbot.github.io/ggspatial/reference/annotation_scale.html)
  function from `ggspatial` to add a scale to the plot. If NULL, no
  scale added.

- compass:

  list or NULL, default
  `list(style = ggspatial::north_arrow_fancy_orienteering, location = "br")`.
  Arguments passed to
  [`annotation_north_arrow`](https://paleolimbot.github.io/ggspatial/reference/annotation_north_arrow.html)
  function from `ggspatial` to add a compass to the plot. If NULL, no
  compass added. Note that calls to alternative styles, like
  [`north_arrow_fancy_orienteering`](https://paleolimbot.github.io/ggspatial/reference/north_arrow_orienteering.html)
  cannot have `()` after them, like they would if called directly.

- ask:

  logical, default TRUE. Should the function ask for confirmation before
  sourcing github code?

## Details

Currently, this only works for simple, sample specific facets.
Coordinates for pie charts should be provided as an `sf` object, where
one column, named for the facet being plotted, provides the subfacet
level names matching those in the assignments. Additional sf objects can
be provided, which will also be plotted. Note that there is no need to
standardize the CRS across the objects, since each will be transformed
to match the sample coordinates.

## Examples

``` r
if (FALSE) { # \dontrun{
# get an sf of the sampling locations
lat_long <- data.frame(SMR = c(44.365931, -121.140420), 
                       CLF = c(44.267718, -121.255805), 
                       OPL = c(44.485958, -121.298360), 
                       ASP = c(43.891693, -121.448360), 
                       UPD = c(43.891755, -121.451600), 
                       PAL = c(43.714114, -121.272797)) # coords for points
lat_long <- t(lat_long)
colnames(lat_long) <- c("lat", "long")
lat_long <- as.data.frame(lat_long)
lat_long$pop <- rownames(lat_long)
psf <- sf::st_as_sf(as.data.frame(lat_long), coords = c("long", "lat"))
psf <- sf::`st_crs<-`(psf, "EPSG:4326")

# get the assignments (STRUCTURE-like results)
assignments <- plot_structure(stickSNPs, "pop", alpha = 1) 

# get a map of oregon as a background from the maps package. 
background <- rnaturalearth::ne_states(iso_a2 = "US", returnclass = "sp")
background <- sf::st_as_sf(background)
background <- background[background$name %in% "Oregon",]

# make the plot
plot_structure_map(assignments, k = 2, facet = "pop", pop_coordinates = psf,
                   layers = list(ggplot2::geom_sf(data = background, 
                                     ggplot2::aes(fill = "example")),
                   ggnewscale::new_scale_fill()), radius_scale = .2)
} # }
```
