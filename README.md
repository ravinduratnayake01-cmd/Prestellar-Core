# Prestellar-Core

constants.f95          # Physical constants
prestel.f95            # Hydrodynamics solver
jeans_analysis.f95     # Jeans instability Analysis

main_acoustic.f95      # Acoustic wave test 
main_jeans.f95         # Jeans instability test 
main_freefall.f95      # Free-fall collapse test

Makefile               # Compilation rules

plot_acoustic.py       # Visualization
plot_jeans.py          
plot_freefall.py       


## Compiling

Compile the Acoustic Test:
	make acoustic
Compile Jeans Instability Test:
	make jeans
Compile Free-Fall Collapse Test:
	make

## Running Simulations
./acoustic_test.exe
./jeans_test.exe
./dense_core.exe

## Output Files

wave_analytic.dat           	# Acoustic wave snapshots
wave_t0.dat
wave_thalf.dat
wave_tT.dat
jeans_growth_a.dat       	# Jeans Growth/oscillation data
jeans_growth_b.dat
jeans_marginal.dat       	# Marginal stability case
freefall_initial.dat     	# Initial free-fall state
freefall_final.dat       	# Final collapse state


## Visualization

After a simulation run, generate plots with:
plot.py
plot_results.py
