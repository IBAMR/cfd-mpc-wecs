#!/bin/bash

# A script to download the required BEMIO files to parse ANSYS AQWA data and generate a state space representation of the WEC

# Read_AQWA.m
wget --no-check-certificate --content-disposition https://raw.githubusercontent.com/WEC-Sim/WEC-Sim/master/source/functions/BEMIO/Read_AQWA.m --output-document=Read_AQWA.m

# Normalize.m
wget --no-check-certificate --content-disposition https://raw.githubusercontent.com/WEC-Sim/WEC-Sim/master/source/functions/BEMIO/Normalize.m --output-document=Normalize.m 

# Plot_BEMIO.m
wget --no-check-certificate --content-disposition https://raw.githubusercontent.com/WEC-Sim/WEC-Sim/master/source/functions/BEMIO/Plot_BEMIO.m --output-document=Plot_BEMIO.m

# Radiation_IRF.m
wget --no-check-certificate --content-disposition https://raw.githubusercontent.com/WEC-Sim/WEC-Sim/master/source/functions/BEMIO/Radiation_IRF.m --output-document=Radiation_IRF.m 

# Radiation_IRF_SS.m
wget --no-check-certificate --content-disposition https://raw.githubusercontent.com/WEC-Sim/WEC-Sim/master/source/functions/BEMIO/Radiation_IRF_SS.m --output-document=Radiation_IRF_SS.m

# Excitation_IRF.m
wget --no-check-certificate --content-disposition https://raw.githubusercontent.com/WEC-Sim/WEC-Sim/master/source/functions/BEMIO/Excitation_IRF.m --output-document=Excitation_IRF.m


