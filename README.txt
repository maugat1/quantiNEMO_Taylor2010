Malcolm Augat <ma5ke@virginia.edu>

Data, code, and files for my Fall 2010 semester first-year rotation in the Taylor lab. The aim of this project is to model inbreeding and selection in a metapopulation, and the work of the project is all code and thought-based. The modeling framework I am using (and occasionally extending) is quantiNEMO.

Roughly, the major organization of this project is as such:
    /               This (toplevel) directory
    /README.txt     This file
    /src/           The (modified) quantiNEMO source code
    /doc/           Documentation for quantiNEMO and other things
    /results/       The results and output from runs
    /pres/          My presentations of the project 

TODO:
    - test colonization code
    - allow different colonization models when colonization is different than
        migration (only try to colonize extinct patches or survival prob for
        colonizers going to extant patch)
    - change patch extinction (and modify colonization stuff above) to allow
        group selection
