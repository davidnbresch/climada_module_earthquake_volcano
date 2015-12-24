climada_module_earthquake_volcano
=================================
formerly named climada_module_eq_global until 20151224

Global earthquake model. See eq_global/doc for documentation.

In order to grant core climada access to additional modules, create a folder ‘modules’ in the core climada folder and copy/move any additional modules into climada/modules, without 'climada_module_' in the filename. 

E.g. if the addition module is named climada_module_MODULE_NAME, we should have
.../climada the core climada, with sub-folders as
.../climada/code
.../climada/data
.../climada/docs
and then
.../climada/modules/MODULE_NAME with contents such as 
.../climada/modules/MODULE_NAME/code
.../climada/modules/MODULE_NAME/data
.../climada/modules/MODULE_NAME/docs
this way, climada sources all modules' code upon startup

see climada/docs/climada_manual.pdf to get started

copyright (c) 2014, David N. Bresch, david.bresch@gmail.com all rights reserved.
