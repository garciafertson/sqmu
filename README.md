This program is a collection of utlities to analyze Squeezemeta results.


INSTALLATION

Clone the git repository into your computer:

$ git clone https://github.com/garciafertson/sqmu.git


Create and activate a conda enviroment using the provided yml file:

$ conda env create -f conda_env.yml
$ conda activate sqmu


Install keggrest package using pip:

$ pip install keggrest


Add a sym link into the enviroment's bin folder (alernatively add the project's
bin folder in your $PATH)

$ cd $CONDA_PREFIX
$ ln -s /path/to/project/sqmu/bin/sqmu.py bin/.



USER GUIDE

This program performs a simple exploratoty analysis on Squeezemeta v1.3 results. It includes three analysis types focused on KO and taxonomy annotation provided by SQM results.

For more information please type:

    sqmu.py tax_module -h

    sqmu.py tax_network -h

    sqmu.py graftm_abundance -h
