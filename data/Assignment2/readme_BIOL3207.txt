"Data, analysis script and associated files for:"			
"Clark TD, Raby GD, Roche DG, Binning SA, Speers-Roesch B, Jutfelt F, Sundin J (in prep) Ocean acidification does not impair the behaviour of coral reef fishes"			
 			
"Data collected by TDC, GDR, DGR, SAB, BSR, FJ and JS. Please refer to the manuscript for data collection methods and statistical analyses. For questions or to notify the authors if any errors are identified in the data, please contact Tim Clark (t.clark@deakin.edu.au), Graham Raby (graham.d.raby@gmail.com), and/or Dominique Roche (dominique.roche@mail.mcgill.ca)."			
			
			
			
			
##################			
# DATA files			
##################			

Note: Fish standard length (SL) data were measured directly in some experiments (e.g., 2016) and estimated from regression relationships in others.  
For example, the Acanthochromis polyacanthus used in choice flume experiments in 2015 arrived in the lab at a very uniform size, and batches of fish 
were sacrificed from holding tanks at various timepoints throughout the study.  Regressions of SL vs. time were then used to estimate the SL of each 
of the fish used in experimental trials, as we did not want to risk harming these fragile fish by taking length measurements.  Fish SL in the 2014 flume 
data was measured from the video footage and then a post hoc calibration relationship was formulated for measured SL and footage-derived SL, to enable 
an estimate of the former from the latter for each of the experimental fish.  Excess decimal places have been retained in the dataset in instances where 
SL has been calculated rather than measured. This also explains the repeat values of SL for different individuals with excess decimal places.
			
### OA_flumedat_20190302.csv ###			
			
columnHeading		description	
-------------		-----------	
loc			Location, and year, where the data were collected. AIMS = Australian Institute of Marine Science; LIRS = Lizard Island Research Station
species			Species name: Acantho = Acanthochromis; Ambon = Pomacentrus amboinensis; Chromis = Chromis atripectoralis; Humbug = Dascyllus aruanus; Lemon = Pomacentrus moluccensis	
treatment		Elevated CO2 [CO2] (850-1,050 µatm) or control [Control] (400 - 450 µatm) groups
animal_id		Fish identity
SL			Standard length of the fish in mm
size			Size grouping of the fish, separated at 15 mm standard length into 'big' or 'small'
attraction		Percent of time the fish spent in the side of the flume containing the predator chemical cue, mean of all minutes for which side preference was calculated
comment			Comment with notes on the origin of the data
			
*** Data used to generate Fig. 1 and Fig. 2 in the manuscript			
			
			
			
### OA_activitydat_20190302.csv ###			
			
columnHeading		description	
-------------		-----------	
loc			Location, and year, where the data were collected. AIMS = Australian Institute of Marine Science; LIRS = Lizard Island Research Station
species			Species name: acantho = Acanthochromis; Ambon = Pomacentrus amboinensis; Chromis = Chromis atripectoralis; Humbug = Dascyllus aruanus; Lemon = Pomacentrus moluccensis	
treatment		Elevated CO2 [CO2] (850-1,050 µatm) or control [Control] (400 - 450 µatm) groups
animal_id		Fish identity
sl			Standard length of the fish in mm
size			Size grouping of the fish, separated at 15 mm standard length into 'big' or 'small'
activity		Number of seconds the fish was active per minute, averaged across the duration of the trial
comment			Comment with notes on the origin of the data
			
*** Data used to generate Fig. 4 in the manuscript			
			
			
			
### LIRS 2014 lat data.csv ###			
			
columnHeading		description	
-------------		-----------	
treatment		"Elevated CO2 [CO2] (850-1,050 µatm) or control [Control] (400 - 450 µatm) groups"	
species			Species name: Acanthochromis = Acanthochromis; Ambon = Pomacentrus amboinensis; Chromis = Chromis atripectoralis; Humbug = Dascyllus aruanus; Lemon = Pomacentrus moluccensis
fishID			Fish identity
daysExposed		"Numbers of days in elevated CO2 conditions (850-1,050 µatm)"	
daysUntilExp		Numbers of days in experimental conditions (Control or CO2)	
tank			Tank number
tankFishIdentifier	Unique fish identifier		
left			Number of turns to the left (out of 10)
right			Number of turns to the right (out of 10)
max			Maximum number oif turns to one side (left OR right)
total			Total number of turns
relLat			Relative lateralisation score
absLat			Absolute lateralisation score
number			Trial number (not in chronological order)
score			Did the fish turn right [1] or left [0] in the trial (i.e. row)?
			
*** Data used to generate Fig. 5 in the manuscript			
			
			
			
### AIMS 2015 lat data.csv ###			
			
columnHeading		description	
-------------		-----------	
treatment		"Elevated CO2 [CO2] (850-1,050 µatm) or control [Control] (400 - 450 µatm) groups"	
population		"Collected in the wild [Wild] or aquarium reared at Reef HQ in Cairns, QLD [RHQ]"	
fishID			Fish identity
experimenter		Name of experimenter	
number			Trial number (not in chronological order)
score_offset		Did the fish turn right [1] or left [0] in the trial (i.e. row) when the barrier at one end of the detour test runway was offset?	
score			Did the fish turn right [1] or left [0] in the trial (i.e. row)?
			
*** Data used to generate Fig. 5 in the manuscript			
			


### Time_in_cue_2+2.csv ###			
			
columnHeading		description	
-------------		-----------			
Location		Location, and year, where the data were collected. AIMS = Australian Institute of Marine Science; LIRS = Lizard Island Research Station
Species			Abbreviated species name; see paper for details
Treatment		Elevated CO2 or control groups
Fish_ID			Individual number associated with each fish
SL_mm			Standard length of the fish in mm
Size_group		Categorical grouping of each fish into "small", "medium" or "big" grouping based on standard length
Time_in_cue		Percentage of time spent in the size of the choice flume containing predator chemical cue
Comment			Note on which sections of the data file were taken for the analysis (equal to 2 + 2 minutes of data [i.e. before and after switching sides of the predator chemical cue in the choice flume])
			
*** Data used to generate Fig. 3 in the manuscript		



### Water_chemistry_data.csv ###			
			
columnHeading		description	
-------------		-----------
location		Location where the measurements were carried out: LIRS = Lizard Island Research Station; AIMS = Australian Institute of Marine Science
year			Year when the measurements were carried out.
date			Date when the measurements were carried out. MONTH_day
time			Time when the measurements were carried out.
treatment		Treatment group: control = present-day CO2 concentration; CO2 = end of century CO2 concentration. See the paper for numerical concentrations.
variable		Type of measurement: temperature (degrees Celsius); pCO2 (uatm); alkalinity (umol kg-1)
value			Value of the measurement (see 'variable' for units of measurement)

*** Data not used to generate any figure
	


			
##################			
# SCRIPT files			
##################			
			
			
### OA_choice flume analysis script.R ###			
			
Annotated R script to examine the effects of elevated pH on attraction to odour cues and generate associated figures.			
			
			
			
### OA_activity analysis script.R ###			
			
Annotated R script to examine the effects of elevated pH on activity levels and generate associated figures.			
			
			
			
### OA_lateralisation_script.R ###			
			
Annotated R script to examine the effects of elevated pH on individual- and population-level lateralisation and generate associated figures.			



### OA_bootstrapping_script.R ###

Annotated R script for bootstrapping simulations to compare the data collected in this study with those of previous studies (see Electronic Supplementary material).  
