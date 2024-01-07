 ----------------------------------------------------------
	                Chocked Class 

 - written using the python package skysurvey By Mickael Rigault 
   (https://skysurvey.readthedocs.io/en/latest/#skysurvey)

 - Thanks to Jason for sharing codes of the afterglow model
					    	               
 ----------------------------------------------------------

How to run the codes:

From terminal launch python scripts:


1) chockedday.py

In the python script chockedday.py the Target is generated as described in the
model presented in Rabinak&Waxman 2011 doi:10.1088/0004-637X/728/1/63				

Intensity f(nu, t) depends on several parameters fixed in the code
Is it possible to change the Mass, the Energy and the Radius of the source from the code


2) chocked_class.py

Using chockedday.py this script defines a Class of the target that will be used in the last script 


3) UVchocked.py

This script put together the target (Chocked) and the survey (ULTRASAT) creating a DATASET that would have collected observing this target with this survey





