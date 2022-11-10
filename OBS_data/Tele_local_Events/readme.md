**Folder information and work flow**

 **Data** 

    1. This folder contains the 2 tele events and 2 local events with time errors at station H36
 
 **event_lst.txt**
 
    1. This files contain the event information with magnitude more than 6 using *wilber3 (IRIS)*
  
 **get_Data_change_header_IRIS.sh**
 
    1. it will help to extract the events information from event_lst.txt (previos file)
    
    2. in step 2 , data will be copied for each station from directory with specfic day
    
    3. we will change the header information, (insert events informations, triggering time etc)
    
    4. cut the event data with an appropriate length to reduce the size and avoiding mulitple evnets.
    
    
**caltime-taup.sh**

    1. This script will help to calculate the theoratical arrival times on each station.
  
    2. Each station must have the events information (including latitude, longitude, depth and distance from station)
    
    3. Store the arrival and events information in respective events waveforms folder with file name *phase.arr*
    
**plt.eqprf.sh**

    1. This script will helps to plot events profiles with both predicted and observed travel times
    
    2. Both predicted file (phase.arr) and observed file (time_picks.txt) must be in the event wf folder
    
    3. adjust the scale based on events distance and arrival time. 
    
 **get_ppicks_for_each_Stations**
    
    1. This script will help to calcuate the differential time predicted and onserved arrival times for specific station and
    save in output to further use in plots folder.
    
    2. user need to define the path and station in script
   
**Plots**

    1. In this folder, *dat files are the differential times of each stations for all events.
    
    2. local referes to local events and tele referes to teleseismic earthquakes.
    
    3. Differential times is script to plots these diff.time w.r.t time and time error
    
    4. jpeg file is final output of waveforms inspection test

**misc.scripts**
    not complete yet

    
 
