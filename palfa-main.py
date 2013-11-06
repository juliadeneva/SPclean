#!/usr/bin/python
# Steps for accessing SQL Server from linux machines using python:
# 1) Need python 2.4 installed
# 2) Download freetds source from http://www.freetds.org/
# 3) ./configure --with-tdsver=8.0 --enable-msdblib 
# 4) make, make install
# 5) Download pymssql from http://pymssql.sourceforge.net/
# 6) Run "python setup.py install" for installation
# 7) Test connection with this script

import pymssql
import sys
import os
from glob import glob
from subprocess import *
#import sm

mjd_low = 54147
mjd_high = 54200

smfile = 'giant.alfa.sm'
SIZEMAX = 400000000
mode = 'match'

fskip = open("skip.txt","r")
badpointings = fskip.read();
# print badpointings
# sys.exit()

mjd = mjd_low
while mjd <= mjd_high:
    
    # Open direct connection to the database
    print "Opening DB connection"
    con = pymssql.connect(host='scidata2.tc.cornell.edu',user='deneva',password='@astro.cornell.edu',database='arecibo')
    cur = con.cursor()
    
    # Get the pointing names for this MJD
    query = "SELECT source_name FROM Headers WHERE source_name LIKE 'G%' AND beam_id = 0 AND timestamp_mjd > "+str(mjd)+" AND timestamp_mjd < "+str(mjd+1)
    ret=cur.execute(query)

    print "Closing DB connection"
    #con.commit()
    con.close()
    
    if cur.rowcount == 0:
        print "No SP plots for MJD "+str(mjd)
        mjd = mjd + 1
        continue

    # Make a sub-directory for this mjd
    print "\n*** Working on mjd "+str(mjd)+"; total pointings: "+str(cur.rowcount)+" ***\n"
    if not os.path.exists(str(mjd)):
        pr = Popen("mkdir "+str(mjd), shell=True, stdout=PIPE)
        pr.wait()
    os.chdir(str(mjd))
    # Copy the delays and dmlist files into current directory
    pr = Popen("cp ../delays .", shell=True, stdout=PIPE)
    pr = Popen("cp ../dmlist_palfa_ctc .", shell=True, stdout=PIPE)
    pr = Popen("cp ../giant.alfa.sm .", shell=True, stdout=PIPE)
        
    # Handle each set of 7 files separately
    pointings = cur.fetchall()
    for row in pointings:
        spfiles = ""
        pointing = row[0]
        print "\n*** Working on pointing "+row[0]+" ***\n"

        # If in skip list (i.e. there's something pathological about it), skip
        if badpointings.find(row[0]) > -1:
            print "Pointing "+row[0]+" in skip list, skipping."
            continue

        # If plots already made, skip
        if (glob('*'+row[0]+'*.gp7.*') and glob('*'+row[0]+'*.gp7clean.*')) or glob('*'+row[0]+'*.gpclean.*'):
            print "Plots already made for this pointing, skipping."
            continue

        # Open direct connection to the database
        print "Opening DB connection"
        con = pymssql.connect(host='scidata2.tc.cornell.edu',user='deneva',password='@astro.cornell.edu',database='arecibo')
        cur = con.cursor()
        
        # Get zip filenames and sizes for this pointing
        query1 = "SELECT source_name,beam_id,timestamp_mjd,observation_time,zip_sift_filename,zip_sift_filesize FROM Headers, Pulses_files WHERE source_name LIKE '"+row[0]+"' AND timestamp_mjd > "+str(mjd)+" AND timestamp_mjd < "+str(mjd+1)+" AND Headers.header_id = Pulses_files.header_id ORDER BY beam_id"
        query2 = "SELECT source_name,beam_id,timestamp_mjd,observation_time,zip_sift_filename,zip_sift_filesize,zip_sift_filedata FROM Headers, Pulses_files WHERE source_name LIKE '"+row[0]+"' AND timestamp_mjd > "+str(mjd)+" AND timestamp_mjd < "+str(mjd+1)+" AND Headers.header_id = Pulses_files.header_id"
                
        if mode == 'sift':
            pass
        elif mode == 'match':
            query1 = query1.replace('sift','best')
            query2 = query2.replace('sift','best')
        else:
            print 'Mode should be "match" or "sift"; mode found: '+mode
            sys.exit()

        try:
            ret=cur.execute(query1)
        except pymssql.DatabaseError:
            print "Caught a DatabaseError, skipping current query."
            print "Query: "+query1
            continue
        except:
            print "Unexpected error: ",sys.exc_info()[0]
            print "Query: "+query1
            raise
        
        # If more than 7 beams satisfy query, it must have been processed more than once, for now skip
        if cur.rowcount > 7:
            print "Skipping pointing "+row[0]+" from MJD "+str(mjd)+": "+str(cur.rowcount)+" beams returned"
            continue
        
        beam_conditions = []
        stampfile = row[0]+".tstamps"
        timestamps = file(stampfile,"w")
        for row in cur.fetchall():
            if row[5] > SIZEMAX:
                print "Skipping file "+row[4]+" with size "+str(row[5])
                beam_conditions.append("beam_id != "+str(row[1]))
                spfiles = spfiles+"dummy "
                timestamps.write("0\n")
                continue
            elif not row[5]:
                print "Pulse file missing."
                beam_conditions.append("beam_id != "+str(row[1]))
                spfiles = spfiles+"dummy "
                timestamps.write("0\n")
                continue
            else:
                spfiles = spfiles+row[4]+" "
                timestamps.write(str(row[2])+"\n")
        timestamps.close()

        if len(beam_conditions):
            beam_conditions = ' AND '.join(beam_conditions)
            query2 = query2+" AND "+beam_conditions
        query2 = query2+" ORDER BY beam_id"
        #beam_conditions = '('+beam_conditions+')'
        
        #print beam_conditions
        
        try:
            ret=cur.execute(query2)
        except pymssql.DatabaseError:
            print "Caught a DatabaseError, skipping current query."
            print "Query: "+query2
            continue
        except:
            print "Unexpected error: ",sys.exc_info()[0]
            print "Query: "+query2
            raise
        
        # print cur.rowcount
        # print cur.fetchall()

        # Save zip files to disk and relative timestamps to file
        for row in cur.fetchall():
            #print "Writing file "+row[4]
            f_out=open(row[4],'wb')
            f_out.write(row[6])
            f_out.close()

        print "Closing DB connection"
        #con.commit()
        con.close()

        # Unzip all the SP files and delete the zipped ones
        print "Unzipping files..."
        zipfiles = glob('*.zip')
        for zipfile in zipfiles:
            pr = Popen("unzip "+zipfile+" -x 'spectra*'", shell=True, stdout=PIPE)
            pr.wait()
        pr = Popen("\\rm *.zip", shell=True, stdout=PIPE)
        pr.wait()

        if mode == 'sift':
            spfiles1 = spfiles.replace("zip","dat")
            spfiles2 = spfiles.replace("zip","dat2")
        elif mode == 'match':
            if glob('*match.best'):
                spfiles1 = spfiles.replace(".best.zip",".match.best")
                spfiles2 = spfiles.replace(".best.zip",".match.best2")
            else:
                spfiles1 = spfiles.replace(".zip","")
                spfiles2 = spfiles.replace(".zip","2")

        #print "spfiles: "+spfiles
        #print "spfiles1: "+spfiles1
        #print "spfiles2: "+spfiles2
        plotfile1 = spfiles1.split(" ")[0]
        if mode == 'sift':
            plotfile1 = plotfile1.split(".pulse_0")[0]
            plotfile1 = plotfile1+".gp7.cluster.gif"
            plotfile2 = plotfile1.replace("gp7","gp7clean")
        elif mode == 'match':
            plotfile1 = plotfile1.split(".best")[0]
            plotfile1 = plotfile1+".gp7.match.gif"
            plotfile2 = plotfile1.replace("gp7","gp7clean")
        # print plotfile

        # Clean SP files
        if mode == 'sift':
            pr = Popen("../spclean "+spfiles1+" "+stampfile+" 1", shell=True, stdout=PIPE)
        elif mode == 'match':
            pr = Popen("../spclean "+spfiles1+" "+stampfile+" 0", shell=True, stdout=PIPE)
        pr.wait()
        
        # Run Supermongo and make plots from the cleaned SP files
        smcomm = []
        smcomm.append('define smfile "'+smfile+'"\n')
        smcomm.append('macro read "$!smfile"\n')
        #smcomm.append('device gif '+plotfile1+'\n')
        #smcomm.append('sevenbeams '+spfiles1+' dmlist_palfa_ctc 0.000064 \n')
        #smcomm.append('hardcopy \n')
        smcomm.append('device gif '+plotfile2+'\n')
        smcomm.append('sevenbeams '+spfiles2+' dmlist_palfa_ctc 0.000064 \n')
        smcomm.append('hardcopy \n')

        spfiles2 = spfiles2.split(" ")
        for spfile in spfiles2:
            if mode == 'sift':
                plotfile = spfile.split(".pulse_0")[0]
                plotfile = plotfile+".gpclean.cluster.gif"
            elif mode == 'match':
                plotfile = spfile.split(".best")[0]
                plotfile = plotfile+".gpclean.match.gif"
            smcomm.append('device gif '+plotfile+'\n')
            smcomm.append('singlebeam '+spfile+' dmlist_palfa_ctc 0.000064 \n')
            smcomm.append('hardcopy \n')
            
        smcomm.append('quit \n')
        smgpout = file(pointing+'.smgp.out','w')
        smgperr = file(pointing+'.smgp.err','w')
        SMG = Popen('sm',shell=True,stdin=PIPE,stdout=smgpout,stderr=smgperr)
        SMG.stdin.writelines(smcomm)
        SMG.wait()
        
        # Clean up
        pr = Popen("\\rm *dat* *best* dummy*",shell=True, stdout=PIPE)
        pr.wait()
        #sys.exit()
    # End for loop over pointings

    # Clean up
    pr = Popen("\\rm delays dmlist* *.sm",shell=True, stdout=PIPE)
    pr.wait()
    os.chdir("..")
    
    mjd = mjd + 1
# End while loop over MJDs
