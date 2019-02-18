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

----------


## Directories in this Repository

### Genomic Data Collection Scripts
This directory contains all scripts that are used in collecting sequence data and their annotations, e.g. the script used in collecting the genomes from the NCBI databases and assigning Pfam annotations to them for the Templeton grant. 

### Phylostratigraphy Dating Species and Pfams
All scripts and files involved in dating domains/genes go here. A master phylogenetic tree is kept here in Newick Format for all species used in the work for the Templeton grant.

### Protein Sequence Metrics Scripts
All scripts used in generating metrics for proteins and their domains go in this directory. This directory is broken into two subdirectories: one for full protein analysis, one for domain analysis. Metrics include, but are not limited to, clustering, intrinsic structural disorder, aggregation propensity, and amino acid composition.

### Simple MySQL Manipulating Scripts
All scripts that are used to perform basic operations on MySQL databases go here. These include scripts that create backups of tables, add indices to table rows, and delete specific entries. 

---------------

## Dependancies 

### Python3
This repository contains many scripts that require Python3 to run. Fusion doesn't allow for the global installation of Python3, so the user should download [anaconda3](https://www.anaconda.com/distribution/), unpack/make/configure the files, and then should specify the path to the python3 executable (located in bin) when running code. 

Python3 also has an executable in its bin called conda which allows the user to download python packages from [the anaconda cloud](https://anaconda.org/). Examples of imporant modules that are available to download using conda are the following:

#### [mysql.connector](https://anaconda.org/anaconda/mysql-connector-python)
Allows user to connect to mysql

#### [ftputil](https://anaconda.org/conda-forge/ftputil)
Allows user to connect to an online ftp site

#### [matplotlib](https://anaconda.org/conda-forge/matplotlib)
Package for all your plotting need

#### [datetime](https://anaconda.org/trentonoliphant/datetime)
Basic time-keeping module


