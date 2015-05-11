#/bin/bash
(cd Reconstruction; ./add_reconstruction_table.sh)
(cd CalibrationFiles; ./apply_calibration_files.sh)
(cd WMC_Geometry; ./apply_geometry.sh)
