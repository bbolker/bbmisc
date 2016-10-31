# celebrate

A silly package for audio enhancement of statistical model summaries

## Installation

- **Linux only**: install a recent version of the PortAudio library (e.g. for Debian/Ubuntu etc., `sudo apt-get install portaudio19-dev`) *before* installing (or re-installing) the `audio` package. (When you install from source you should see `checking for working PortAudio (API>=2)... yes` as the build messages scroll by ...)
- `devtools::install_github("bbolker/bbmisc",sub="celebrate")`
- `library(celebrate)`
- `?celebrate`

## To do

- audio is not working properly under VirtualBox on Ubuntu 14.04 (sigh). Try falling back to `/usr/bin/play` on Linux (available through `sox` package)?
- implement for `print.summary.glm` as well
