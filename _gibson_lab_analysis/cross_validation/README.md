# Scripts ran & setup

In order to print out "x0" in the cross-validation forward simulation loop for MDSINE2, we ran the following scripts:
- simulation_data.Rmd
- forecasting_metrics.Rmd

In order to run the above scripts, we set up the R packages this way:
```bash
  Rscript -e "install.packages('fido', repos='https://cloud.r-project.org')"
  Rscript -e "devtools::install_local('tfPaper')"
  # skipped fido2 installation, due to a mismatch between fido1's exposed interfaces and the author's fido2 implementation.
  Rscript -e "devtools::install_github('gathanei/xyz')"
  Rscript -e "devtools::install_github('krisrs1128/mbtransfer')"
  Rscript -e "devtools::install_local('mdsine')"
  Rscript -e "mdsine::install_mdsine()"
```

# Notes

1) "tfPaper" has been modified to exclude fido2 (due to a failed installation, see section 2 above)
2) simulation_data.Rmd generates sim_input*.rda files in the home ~/Download directory. 
The forecasting_metrics.Rmd script assumes that these files have been moved/symlinked into ./tf_sim.
Alternatively, one can download the pre-generated files here: (https://www.dropbox.com/scl/fi/9phq3v0f58bykd8czfvch/tf_sim.tar.gz?rlkey=zsxau7ni9tbz4tdy3vta3wq1k&st=whml3ffy&dl=0)
