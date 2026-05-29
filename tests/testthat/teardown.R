# close the NULL pdf device opened in setup.R; remove any straggler Rplots.pdf
# that slipped past it (e.g. opened from a forked future worker).
tryCatch(grDevices::dev.off(), error = function(e) NULL)
unlink("Rplots.pdf")
