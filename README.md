# TOROS-Field-Aligner
This function takes in two fields of detected stars and matches their photutils DAOStarFinder starIDs together.

## What code is doing
The user inputs 2 csv files containing photometric data of stars taken from two different fields. First, the pixel deviation between fields is computed by finding the brightest stars in both fields and matching them. This process is run
until the pixel deviation stablizes at some value, at which the deviation is subtracted from the *comparator* field so that stars are approximately all at the same pixel coordinates between frames. Next, the two fields are sent into the
star matching algorithms after several steps of filtering. The starmatcher picks a star from the *comparator* field and iterates through all stars in the *master* field calculating the distance between the *cstar* and *mstar*. The nearest *mstar* is matched with the *cstar's* starID, and if there is a tie then the magnitude difference of the stars will be used to break it. This process is repeated by picking an *mstar* and comparing it to all *cstars* using the same method as before. The only stars that are retained are those that matched to the same starIDs between the *master* field and the *comparator* field in both iteration processes. This matched id list is then returned to the python code, where a new astropy table is created and written to disk.

## Setup 
***1.*** Place dllmain.cpp, setup.py, and starmatchingalg.cpp into same folder. \
***2.*** In Native tools cmd, run the following script: [path to folder you made in step one]>python setup.py build_ext --inplace \
***3.*** Before running TOROSFieldAligner in the .py file, you must run sys.path.append([path to folder you made in step one]) so that python can source starmatcher.

## Use
This section will explain the function parameters. (**function_parameter** [data type]: Explanation of parameter.) \
\
\
**master_file** [str]: The field whose starIDs are prefered. Other fields will have their starIDs matched to the IDs in the master. This is the field that should stay constant if you use this function many times. \
\
**comparator_file** [str]: The field whose starIDs get overwritten by the master IDs. \
\
**alligned_file** [str]: The file name you want for the new astropy table that gets written containing all matched stars (ignore my typo I dont wanna rewrite everything ;-;) \
\
**dist** [float]: The maximum distance to be considered by the starmatcher when iterating through stars. (keep this number low) \
\
**pixel_cutoff** [float]: Similar to dist but applied during the offset calibration. \
\
**snr_threshold** [int]: Minimum SNR you want a star to have. Assumes poisson noise (ie SNR ~ flux / sqrt(flux)). This is a filter applied before using starmatcher. 

### Versions
I've included two versions of the program. TOROSFieldAligner.py uses C++ for the actual matching. In the event the C++ stuff does not work, I have also included a version of the program (TOROSFieldAlignerSlow.py) that processes everything through python. It is much slower, but has the benefit of not needing extra installation steps to use.

