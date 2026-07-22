#!/bin/bash

#Path to executables
EXEC_PATH="/home/kelsey/GAPS_TrackerEnergyDep/build/"
#Path to data
DATA_PATH="/home/kelsey/simulations/simdat/flight/251221/FPSI/starlink251221_"
#TOF Factor
TOF_FACTOR="1.4"
#TKR Factor
TKR_FACTOR="1.15"
#MainloopScaleFactor
MSCALE="5"

echo "=== Data Quality Plots Incoming! ==="
echo "Creating beta directory if it doesn't exist..."
if [ ! -d "beta" ]; then
    mkdir beta
    echo "beta directory created"
else
    echo "beta directory already exists"
fi

echo "Beta proxy..."
"$EXEC_PATH"RecEdepBeta -i "$DATA_PATH" -l 0.2 -u 1.2 -r 2 -f "$TOF_FACTOR" -k "$TKR_FACTOR" -o beta
echo "Other betas..."
"$EXEC_PATH"BetaHisto -i "$DATA_PATH" -l 0.2 -u 1.8 -r 2 -o beta

echo "Angle corrected energy depositions vs beta plots..."
"$EXEC_PATH"PartsRecEdepvB -i "$DATA_PATH" -r 2 -l 0.2 -u 1.8 -b 25 -t 10

echo "Occupancy pdf..."
"$EXEC_PATH"/HoursOccu -i "$DATA_PATH" -m "$MSCALE"
