ReadMe:

RunIUPred2
----------

--------------------------- #
This script is designed to: #
--------------------------- #

	1) Extract protein sequences from a MySQL database
	2) Run the disorder predictor IUPred2 on each protein sequence
	3) (optional) Run IUPred's new function Anchor on each protein sequence
	4) Upload the raw disorder prediction scores and the mean disorder prediction score to the MySQL database where the proteins were extracted from



--------------------------- #
       Dependencies:        #
--------------------------- #

The user will need to run the script with Python3 with the following modules:

	- BioPython: https://anaconda.org/conda-forge/biopython
	- mysql.connector : https://pynative.com/install-mysql-connector-python/

To get Python3, anaconda can be downloaded from: https://www.anaconda.com/download/


The user will also need to have the IUPred2 executables. These can be accessed here: https://iupred2a.elte.hu/download

(Select IUPred2A to download)



