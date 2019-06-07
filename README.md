# ProteinEvolution

This is a repository for submitting communal scripts used in studying protein evolution. 

## Notes About This Repository
### Submitting Communal Scripts
All scripts are communal, so if anyone creates something that’s more than a single-use, user-specific script that would be of value to the rest of the group (such as Scott’s IUPred script or Jason’s), it should be uploaded to this github repository. As a result, no one should be reinventing the wheel when it comes to code! Since these scripts are communal (and this is good practice for your personal scripts as well), the name of the script should be concise and descriptive and should give the user an idea of what it does. Scripts titled things like Test.py or PleaseWork.r are not good choices. 

### Edits
If edits need to be made to a script already on the site, edit it (though it might be a good idea to discuss it with the author first) and commit those changes so that everyone has access to the most recent versions of everything. If editing a script, add your name and the changes made. 

There are various ways to make edits, one of which is using ssh from your terminal. To begin, you should [set up an ssh key](https://help.github.com/en/articles/connecting-to-github-with-ssh) to connect to Github. You can then [clone the repository](https://confluence.atlassian.com/bitbucket/clone-a-repository-223217891.html) you want to edit, edit the scripts you want to make changes to, [add, commit, and push your changes](https://dev.to/juni/git-and-github---must-know-commands-to-make-your-first-commit-333c). 

When committing your changes, you have the option to add a message. This is where you should leave a brief comment on the changes you made. For example:

```
git commit -m "Added more fish to penguin-feeding function -- Sara Willis"
```

To update your local machine with the most recent version of the cloned Github repository, use the command

```
git pull origin master
```

### ReadMe
There is a ReadMe associated with each directory in the Github repository. This readme should contain the names of all the scripts in that directory in alphabetical order. The ReadMe file should have the extension .md so GitHub can read the fancy [syntax](https://help.github.com/articles/basic-writing-and-formatting-syntax/) available. The ReadMe should give an overview of what each script does, how to run it, dependencies for the script and how to install them, the inputs, and anything else the user should know to run your code. 

Individual ReadMe’s may also helpful to submit with individual scripts depending on the complexity of the code. If submitting a ReadMe with your script, here is a [good template](https://gist.github.com/PurpleBooth/109311bb0361f32d87a2) for how to create one.

### Commenting Your Code
Commenting your code is important! This is not just for the sanity of those using your code, this is also for your own. There are a lot of resources out there to help do this effectively, so make use of them! These comments shouldn’t simply restate what the code is doing (which either forces you to recomment the code every time you make an edit, or makes your comments invalid after any change), but should facilitate understanding. 

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

is much more informative


## Directories in this Repository

### Data Tables Creation and Updates
These scripts are used to either update existing values in a MySQL table or to insert values into a new table based from an old table. For example, averaging over particular values already stored in a table and saving the average to a new column, creating indices to map between tables, or creating a new table using preexisting rows from an old table.

### Figures 
Scripts used to generate figures. This is where the massive BoxAndWhiskerPlots_LinearModelSlopes_MetricsVsAge.py script is stored which creates box and whiskers plots as well as csv slope files for the user.

### Genomic Data Collection Scripts
This directory contains all scripts that are used in collecting sequence data and their annotations, e.g. the script used in collecting the genomes from the NCBI databases and assigning Pfam annotations to them for the Templeton grant. 

### HomologyDetection_AbandonedHMMERPipeline
This directory contains all scripts that were written with the intention of creating homology groups using HMMER. Our initial goal was being able to use HMMER's sensitive hmomology detection algorithms to partition the mouse genome into homology groups with each gene segment being represented exactly once. This pipeline turned out to be far more complex than we had anticipated and so we switched to using Pfam annotations to speed up our research. This pipeline is on hiatus pending a grant proposal.

### Linear Models
Where scripts are stored that run linear models on our data. This is a good repository for R files.

### Metrics
All scripts that are used to calculate specific metrics go here. These are scripts that do things like run IUPred and Tango, calculate the mean index of dispersion, calculate the CAI of a genome, calculate the GC% content of a genome, calculate the percent amino acid composition of proteins, etc.

### Phylostratigraphy Dating Species and Pfams
All scripts and files involved in dating domains/genes go here. A master phylogenetic tree is kept here in Newick Format for all species used in the work for the Templeton grant.

### Sequence Alignments
Any scripts that generate multiple sequence alignments.

### Simple MySQL Manipulating Scripts
All scripts that are used to perform basic operations on MySQL databases go here. These include scripts that create backups of tables, add indices to table rows, and delete specific entries. 

### Useful Submodules
These are scripts that contain small submodules that are useful for work in Bioinformatics. Used for things like generating fasta files, running InterProScan, calculating the optimal value of lambda for a Box-Cox transform given a dataset, etc.

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

#### [biopython](https://anaconda.org/anaconda/biopython)
This is an excellent tool for bioinformatics and I highly suggest you get it!

### R
There are scripts in this repository that make use of [R](https://www.r-project.org/). 

For any help on installing packages for R on Fusion where you don't have root access, a quick tutorial can be found [here](https://cmdlinetips.com/2012/05/how-to-install-a-r-package-locally-and-load-it-easily/). The best way to get packages to load in this way can be found in the [CRAN package repository](https://cran.r-project.org/web/packages/).
