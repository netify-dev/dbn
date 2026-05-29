# open a NULL pdf device so tests that call plot.* without an explicit
# device do not leave an Rplots.pdf behind. matched by teardown.R.
pdf(NULL)
