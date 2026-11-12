.PHONY: all touch-sim

TAR_MAKE = Rscript -e 'targets::tar_make(script = here::here("forstmeier_sim", "_targets.R"), store = here::here("forstmeier_sim", "_targets"))'

all:
	$(TAR_MAKE)

## mark the existing simulation results as up to date without rerunning
## (see R/run_corrections_simulation.R); everything downstream still rebuilds
touch-sim:
	FORSTMEIER_TOUCH_SIM=1 $(TAR_MAKE)
