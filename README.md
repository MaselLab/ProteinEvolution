# ProteinEvolution

This  is a common repository for submitting communal scripts used in studying protein evolution. 

## Notes About This Repository
### Submitting Communal Scripts
All scripts are communal, so if anyone creates something that’s more than a single-use, user-specific script that would be of value to the rest of the group (such as Scott’s IUPred script or Jason’s), it should be uploaded to this github repository. As a result, no one should be reinventing the wheel when it comes to code! Since these scripts are communal (and this is good practice for your personal scripts as well), the name of the script should be concise and descriptive and should give the user an idea of what it does. Scripts titled things like Test.py or PleaseWork.r are not good choices. 

### Edits
If edits need to be made to a script already on the site, edit it (though it might be a good idea to discuss it with the author first) and commit those changes so that everyone has access to the most recent versions of everything. If editing a script, add your name to the header of the code and the date that the edits were made.

### ReadMe
There is a ReadMe associated with each directory in the Github repository. This readme should contain the names of all the scripts in that directory in alphabetical order. The ReadMe file should have the extension .md so GitHub can read the fancy [syntax](https://help.github.com/articles/basic-writing-and-formatting-syntax/) that you can use in Github. The ReadMe should give an overview of what each script does, how to run it, dependencies for the script and how to install them, the inputs, and anything else the user should know to run your code. 

Individual ReadMe’s may also helpful to submit with each script depending on the complexity of the code. If submitting a ReadMe with your script, here is a [good template](https://gist.github.com/PurpleBooth/109311bb0361f32d87a2) for how to create one.

### Commenting Your Code
Commenting your code is important! This is not just for the sanity of those using your code, this is also for your own sanity. There are a lot of resources out there to help do this effectively, so make use of them! These comments shouldn’t simply restate what the code is doing (which either forces you to recomment the code every time you make an edit, or makes your comments invalid after any change), but should facilitate understanding. 

Another important point is to use descriptive variable names to help your code narrate the story. Code that looks like:

```
for thing in Thing:
	Thing2.doThingTest(thing)
```

Doesn’t tell us anything and obscures the meaning of the code. Renaming these variables

```
For fish in wheelbarrow:
	pengin.eat(fish)
```

Tells us a lot more about what’s going on.

-------
## Overview of MySQL Tables

For work on the Templeton grant, there are several important data tables on the MySQL database. Here I'll give a rundown these tables, including: where and how the information is stored, the names of the important tables, and where to look for specific data.

As a general overview, the tables contain the protein-coding sequences, their translations, Pfam annotations, and various metrics for over 400 fully-sequenced species. 

### Coding Tables 

There are two tables where all protein-coding nucleotide sequences are stored. These tables are specifically:

* EnsemblGenomes_Coding_Complete
* NCBIGenomes_Coding_Complete

## Protein Tables

Each coding table's sequences have been translated and annotated with the Pfam domains and are stored in separate data tables:

* EnsemblGenomes_Protein_Complete
* NCBIGenomes_Protein_Complete

## Metrics Tables

All metrics for the full protein sequences are stored in their own tables. These metrics include

* Mean ISD (IUPred2, cysteine's included)
* Tango Scores:
  * Number of aggregation-prone regions
  * Density of aggregation-prone regions
  * Number of amino acids in aggregation-prone regions
  * Density of amino acids in aggregation-prone regions
* Percent amino acids composition
  * (ordered A,R,N,D,C,E,Q,G,H,O,I,L,K,M,F,P,U,S,T,W,Y,V)
* Normalized Dispersion Index
  * Trunc and AllFrames are both stored
  * Option 'FILVM' selected from Jason's script
  * Window size = 6
  
These metrics were originally stored in the protein table, but to reduce the size of those tables (which were monsterous), and to improve organization, they are now stored in their own tables

## Pfam metrics table

There is an additional table that stores each Pfam/Gene combination as its own row
