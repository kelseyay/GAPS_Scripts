#!/bin/bash

#Path to executables
EXEC_PATH="/home/kelsey/GAPS_TrackerEnergyDep/build/"
#Path to data
DATA_PATH="/home/kelsey/simulations/simdat/ground/251204/25.10/et"
#post-fix (not prefix for the data)
PFX="_251204_v25.10"
#TOF Factor
TOF_FACTOR="1.50"
#TKR Factor
TKR_FACTOR="1.00"
#MainloopScaleFactor
MSCALE="5"
#Trigger (usually either 2 or 0!)
TRG="0"


echo "=== Data Quality Plots Incoming! ==="
echo "Creating beta directory if it doesn't exist..."
if [ ! -d "beta" ]; then
    mkdir beta
    echo "beta directory created"
else
    echo "beta directory already exists"
fi

echo "Beta proxy..."
"$EXEC_PATH"RecEdepBeta -i "$DATA_PATH" -l 0.2 -u 1.2 -r "$TRG" -f "$TOF_FACTOR" -k "$TKR_FACTOR" -o beta -e "$PFX"
echo "Other betas..."
"$EXEC_PATH"BetaHisto -i "$DATA_PATH" -l 0.2 -u 1.8 -r "$TRG" -o beta -e "$PFX"

echo "Angle corrected energy depositions vs beta plots..."
"$EXEC_PATH"PartsRecEdepvB -i "$DATA_PATH" -r "$TRG" -l 0.2 -u 1.8 -b 25 -t 10 -e "$PFX"

echo "Occupancy pdf..."
"$EXEC_PATH"/HoursOccu -i "$DATA_PATH" -m "$MSCALE"
