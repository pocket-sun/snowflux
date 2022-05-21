for example usage:

* step 1.  
    cd src  
    make && make install  
    cd ..  

* step 2.  
    ./supernova

main.cpp is an example file for usage.
Enviroment variable GLB_DIR should be set to the build path of ``GLOBE``, my ``GLOBE`` 
version is 3.2.18.  
To create the ``.so`` shared lib, replace second line of step 1. with ``make shared``.

**WARRNING**:  
the code is for private usage and I do not follow the LICENSE prescripted by ``SNOWGLOBE`` and
``GLOBE`` authors! If someone wants to use it for constructing public available tools, please make
sure to credit to the ``SNOWGLOBE`` and ``GLOBE`` authors and obey the LICENSE rule given by them!

CLS, 5/21/2022
