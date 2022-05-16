for example usage:

* step 1.  
    cd snowglobe  
    ./supernova.pl argon ar40kt  
    cd ..  

* step 2.  
    cd src  
    make && make install  
    cd ..  

* step 3.  
    ./supernova argon  


As with flux parameters and energy binning information, see source/src/supernova.cpp 
line20 to line36.  
Adjustement is desired through little effort for future plan.  
Enviroment variable GLB_DIR should be set to the build path of ``GLOBE``, my ``GLOBE`` 
version is 3.2.18

**WARRNING**:  
the code is for private usage and I do not follow the LICENSE prescripted by ``SNOWGLOBE`` and
``GLOBE`` author! If someone wants to use it for constructing public available tools, please make
sure to credit to the original ``SNOWGLOBE`` authors and obey the LICENSE rule given by them!
