# celebrate

A silly package for audio enhancement of statistical model summaries

## Installation

- **Linux only**: Getting sound working may be a bit tricky. If you have a working version of `play` on your system (try `which play` from a Unix command prompt, or `system("which play")` from within R) (you can try `sudo apt-get install sox` to get it), this will be tried first. Otherwise, you can try installing a recent version of the PortAudio library (e.g. for Debian/Ubuntu etc., `sudo apt-get install portaudio19-dev`) *before* installing (or re-installing) the `audio` package. (When you install from source you should see `checking for working PortAudio (API>=2)... yes` as the build messages scroll by ...)
- `devtools::install_github("bbolker/bbmisc/celebrate")`
- `library(celebrate)`
- `?celebrate`

## To do

- implement for `print.summary.glm` as well
