## ----fig.show='hold', fig.width=7, fig.height=7-------------------------------
requireNamespace("png", quietly = TRUE)
img <- png::readPNG(paste0(getwd(),"/../inst/extdata/test/ALL_ALL_medres2.png"))
grid::grid.raster(img)

