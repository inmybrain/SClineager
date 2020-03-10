pkgname <- "SClineager"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('SClineager')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("read_sclineager")
### * read_sclineager

flush(stderr()); flush(stdout())

### Name: read_sclineager
### Title: Read datsets for SClineager
### Aliases: read_sclineager

### ** Examples

# See the github page https://github.com/inmybrain/SClineager.



cleanEx()
nameEx("run_sclineager")
### * run_sclineager

flush(stderr()); flush(stdout())

### Name: run_sclineager
### Title: An wrapper function for running SClineager
### Aliases: run_sclineager

### ** Examples

# See the github page https://github.com/inmybrain/SClineager.



cleanEx()
nameEx("sclineager_internal")
### * sclineager_internal

flush(stderr()); flush(stdout())

### Name: sclineager_internal
### Title: The Gibbs sampler for SClineager
### Aliases: sclineager_internal

### ** Examples

# See the github page https://github.com/inmybrain/SClineager.



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
