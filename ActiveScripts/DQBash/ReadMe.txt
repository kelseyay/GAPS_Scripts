Attempt 1 at making a bash script that uses executables from GAPS_TrackerEdep to make all of the DQ plots on the DQ page: https://gaps1.astro.ucla.edu/wiki/gaps/index.php?title=Data_Quality_Version_Checker#plots 

Note: First you must find the gain factors for the TOF and TKR for a given dataset. This should be done and checked manually (until I fix the max bin so that it's fitted instead o-o;;). This is done with the DZEdep3 executable in the following way:
	Run yes TOF, no TKR, both factors = 1 to find the TOF factor. Run yes tof, no TKR with the factor to see if it worked. This plot should be saved and added to DQ. Repeat with TKR. Repeat with both factors both systems on. 

These are the commands to run for the TOF and TKR factors:

(TOF only)
/home/kelsey/GAPS_TrackerEnergyDep/build/DZedep3 -i /home/kelsey/simulations/simdat/flight/251221/FPSI/starlink251221_0 -s 0 -t 1 -f 1 -k 1 -r 2 -l 0.9 -u 0.99 -o test/
Re-run with the correct TOF factor

(TKR only)
/home/kelsey/GAPS_TrackerEnergyDep/build/DZedep3 -i /home/kelsey/simulations/simdat/flight/251221/FPSI/starlink251221_0 -s 1 -t 0 -f <TOF FACTOR> -k 1 -r 2 -l 0.9 -u 0.99 -o test/
Re-run with the correct TKR factor

(BOTH)
/home/kelsey/GAPS_TrackerEnergyDep/build/DZedep3 -i /home/kelsey/simulations/simdat/flight/251221/FPSI/starlink251221_0 -s 1 -t 0 -f <TOF FACTOR> -k <TKR FACTOR> -r 2 -l 0.9 -u 0.99 -o test/
	
	Feed the bash script the TOF and TKR figures and I want something to just execute, automatically, all of the scripts needed for populating the 
	
	Bash script is called DQbash.sh
	Load up the TOF and TKR variables in the bash script as well as the path to the data files!
	
	Run 
	chmod +x DQbash.sh
	
	Then
	./DQbash.sh
	
	
If you want to keep the old name as part of the new name, you can pass the original filename as a variable

file="old_name.txt"
mv "$file" "${file%.txt}_backup.txt"

The easiest way to update the wiki is to copy the 26.03 gallery. Change the upload name to [yadda yadda][v26.XX][F or G for flight or ground]

Upload all of the 9 files in quick succession
