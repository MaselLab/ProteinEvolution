import os, json, sys, csv, requests, re, mysql.connector, datetime, copy, matplotlib
from scipy import stats
from collections import defaultdict
import numpy as np
matplotlib.use('agg')
import matplotlib.pyplot as plt

'''
Author: Jenny James
Uploaded: 11th March 2019  

This function groups species in the dataset into their kingdoms, and return a dictionary of kingdoms and species UIDs

'''

def KingdomUID():

    Kingdom_SpeciesUID = defaultdict(list)

    ### NCBI get species category based on ftp site
    NCBI_genome_categories = defaultdict(list)

    dir = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/"
    kingdoms = ("fungi", "plant", "archaea", "protozoa", "invertebrate", "vertebrate_mammalian", "vertebrate_other")
    for kingdom in kingdoms:
        r = requests.get(dir+kingdom)
        if not r.ok:
          r.raise_for_status()
          sys.exit() 
        NCBI_species = re.findall(r'>(.+)/</a>', r.text)
        for x in NCBI_species:
            NCBI_genome_categories[kingdom].append(x)
    for x in NCBI_genome_categories['vertebrate_mammalian']:
        NCBI_genome_categories['vertebrate'].append(x)
    for x in NCBI_genome_categories['vertebrate_other']:
        NCBI_genome_categories['vertebrate'].append(x)
    del NCBI_genome_categories['vertebrate_mammalian']
    del NCBI_genome_categories['vertebrate_other']


    # This log into MySQL to access Ensembl records and UID
    # input user data
    Database = ''
    User = ''
    Host = ''
    Password = ''

    # This sets up a connection with MySQL using the information you provided above
    cnx = mysql.connector.connect(user = User,
                                password = Password,
                                host = Host,
                                database = Database)

    # a cursor 
    mycursor = cnx.cursor(buffered = True)

    # We start by extracting the speciesUID, NewickSpeciesName, and Source Database from our list of species
    mycursor.execute("SELECT SpeciesUID,NewickSpeciesName,SourceDatabase FROM SpeciesList")
    # fetchall() pulls the results out of our cursor as a list of tuples
    SpeciesResults = mycursor.fetchall()
    # search through the species we've extracted
    for speciesEntry in SpeciesResults:
        # We assign variable names for ease of use
        SpeciesUID,SpeciesName,Source = speciesEntry
        if Source == 'NCBI':
            ### use the NCBI dictionary
            for k, v in NCBI_genome_categories.items():
                if SpeciesName in v:
                    Kingdom_SpeciesUID[k].append(SpeciesUID) 
        ### else: Ensembl source, use ensembl sourceID
        elif Source == 'V93':
            Kingdom_SpeciesUID['vertebrate'].append(SpeciesUID)
        elif Source == 'EnsemblPlantV40':
            Kingdom_SpeciesUID['plant'].append(SpeciesUID)
        elif Source == 'EnsemblMetazoaV40':
            Kingdom_SpeciesUID['invertebrate'].append(SpeciesUID)
        elif Source == 'EnsemblProtistV40':
            Kingdom_SpeciesUID['protozoa'].append(SpeciesUID)
        elif Source == 'EnsemblFungiV40':
            Kingdom_SpeciesUID['fungi'].append(SpeciesUID)    
        
    return Kingdom_SpeciesUID



