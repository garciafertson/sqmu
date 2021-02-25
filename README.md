This program is a collection of utlities to analyze Squeezemeta results.


INSTALLATION

Clone the git repositorie into your computer:

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


