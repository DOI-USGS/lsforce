import sqlite3 as lite
import os
import stat
import glob
import shutil
import numpy as np
import subprocess
#from reviewData import reviewData

def unique_list(seq):  # make a list only contain unique values and keep their order
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]

def setup(event_id, modelfile, stacodes, samplerate, duration, T0, stfilename=None,
          mainfolder='/Users/kallstadt/LSseis/LSPy/INVERSION_FILES',
          database='/Users/kallstadt/LSseis/landslideDatabase/lsseis.db'):
    """
    This function sets up the folder structure and creates some files
    needed for computing Green's functions using CPS

    INPUTS
    event_id = SQLite3 event_id of landslide event
    modelfile = full path to model file (CPS format) to use
    stacodes = list of station ids (from tr.id - NET.STA.LOC.CHAN)
    samplerate = samplerate to use in Green's functions in samples per second
    duration = desired duration in seconds
    (will be converted to number of samples, nearest power of 2)
    T0 = start time of records (negative number recommended) in seconds before impulse time
    mainfolder = full path to main folder (INVERSION_FILES suggested)
    stfilename = name given to a particular run to save these greens functions separately from other
        greens functions

    OUTPUTS
    T0.txt - tiny text file that says what T0 was (in seconds)
    stadistlist.txt - list of stations and corresponding number in dist file
    (for renaming Green's functions)
    Folder structure
    dist.txt - dist file needed for
    modelfile - moves it to the folder
    CPScommands.sh - bash shell script to run
    sacodir - full path of sac directory to put orignal waveforms into
    """

    def nextpow2(val):
        import math
        temp = math.floor(math.log(val, 2))
        return int(math.pow(2, temp+1))
        pass

    #make folders
    if stfilename is not None:
        stfilename = stfilename+'_'
    else:
        stfilename = ''
    evdir = ('%s/EV%s') % (mainfolder, event_id)
    moddir = ('%s/%s') % (evdir, stfilename+modelfile.split('/')[-1].split('.')[0],)
    sacodir = ('%s/%s') % (moddir, 'sacorig')
    sacdir = ('%s/%s') % (moddir, 'sacdata')
    if not os.path.exists(evdir):
        os.mkdir(evdir)
    #except Exception as e:
    #    print(e)
    if not os.path.exists(moddir):
        os.mkdir(moddir)
    if not os.path.exists(sacdir):
        os.mkdir(sacdir)  # folder for renamed sac files to go into
    if not os.path.exists(sacodir):
        os.mkdir(sacodir)  # folder to keep original sac files

    #write T0 file
    f = open(moddir+'/T0.txt', 'w')
    f.write(('%3.2f') % T0)
    f.close()

    #make sure there is only one occurrence of each station in list (ignore channels)
    stacods = [lis[0]+'.'+lis[1] for lis in [stacod.split('.') for stacod in stacodes]]
    stacods = unique_list(stacods)
    stas = [stacod.split('.')[1] for stacod in stacods]
    nets = [stacod.split('.')[0] for stacod in stacods]
    #locs = [stacod.split('.')[2] for stacod in stacods]

    #get distance of each station from database
    connection = None
    connection = lite.connect(database)
    connection.text_factory = str
    dist = []
    stations = []
    with connection:
        cursor = connection.cursor()
        for k in np.arange(len(stas)):
            try:
                cursor_output = (cursor.execute(
                    'SELECT stasource_radius_km FROM sta_nearby, stations ON sta_nearby.station_id = stations.Sid WHERE event_id = ? AND Name = ? AND Network = ?', (event_id, stas[k], nets[k])))
                retrieved_data = cursor_output.fetchall()
                dist.append(retrieved_data[0][0])
                stations.append(stas[k])
            except Exception as e:
                print(e)

    stations, dist = list(zip(*sorted(zip(stations, dist))))

    #write stadistlist.txt
    f = open(moddir+'/stadistlist.txt', 'w')
    for sta, dis in zip(stations, dist):
        f.write(('%s\t%5.1f\n') % (sta, dis))
    f.close()

    #write dist file in free format
    #figure out how many samples
    samples = nextpow2(duration*samplerate)
    f = open(moddir+'/dist', 'w')
    for dis in dist:
        f.write(('%0.1f %0.2f %i %i 0\n') % (dis, 1./samplerate, samples, T0))
    f.close()

    #move copy of modelfile to current directory
    shutil.copy2(modelfile, moddir+'/'+modelfile.split('/')[-1])

    #write shell script to run Green's functions
    f = open(moddir+'/CPScommands.sh', 'w')
    f.write(("""#!/bin/bash
rm %s/*.sac
rm %s/*.sac
hprep96 -HR 0. -HS 0. -M %s -d %s -R -EXF
hspec96 > hspec96.out
hpulse96 -d %s -V -OD -p > Green
f96tosac Green
cp *.sac %s/.
mv *.sac %s/.
""") % (sacodir, sacdir, modelfile, moddir+'/dist', moddir+'/dist', sacdir, sacodir))
    os.chmod(moddir+'/CPScommands.sh', stat.S_IRWXU)

    return moddir


def setup_orphan(event_id, modelfile, stacodes, dists, samplerate, duration, T0, stfilename=None,
                 mainfolder='/Users/kallstadt/LSseis/LSPy/INVERSION_FILES'):
    """
    This function sets up the folder structure and creates some files
    needed for computing Green's functions using CPS for events without event_id's (not in database)

    INPUTS
    event_id = short name of event
    modelfile = full path to model file (CPS format) to use
    stacodes = list of station ids (from tr.id - NET.STA.LOC.CHAN)
    dists = list of source to station distances, same length as stacodes
    samplerate = samplerate to use in Green's functions in samples per second
    duration = desired duration in seconds
    (will be converted to number of samples, nearest power of 2)
    T0 = start time of records (negative number recommended) in seconds before impulse time
    mainfolder = full path to main folder (INVERSION_FILES suggested)
    stfilename = name given to a particular run to save these greens functions separately from other greens functions

    OUTPUTS
    T0.txt - tiny text file that says what T0 was (in seconds)
    stadistlist.txt - list of stations and corresponding number in dist file
    (for renaming Green's functions)
    Folder structure
    dist.txt - dist file needed for
    modelfile - moves it to the folder
    CPScommands.sh - bash shell script to run
    sacodir - full path of sac directory to put orignal waveforms into
    """

    def nextpow2(val):
        import math
        temp = math.floor(math.log(val, 2))
        return int(math.pow(2, temp+1))
        pass

    #make folders
    if stfilename is not None:
        stfilename = stfilename+'_'
    else:
        stfilename = ''
    evdir = ('%s/EV%s') % (mainfolder, event_id)
    moddir = ('%s/%s') % (evdir, stfilename+modelfile.split('/')[-1].split('.')[0],)
    sacodir = ('%s/%s') % (moddir, 'sacorig')
    sacdir = ('%s/%s') % (moddir, 'sacdata')
    try:
        os.mkdir(evdir)
    except Exception as e:
        print(e)
    try:
        os.mkdir(moddir)
    except Exception as e:
        print(e)
    try:
        os.mkdir(sacdir)  # folder for renamed sac files to go into
    except Exception as e:
        print(e)
    try:
        os.mkdir(sacodir)  # folder to keep original sac files
    except Exception as e:
        print(e)

    #write T0 file
    f = open(moddir+'/T0.txt', 'w')
    f.write(('%3.2f') % T0)
    f.close()

    #make sure there is only one occurrence of each station in list (ignore channels)
    stacods = [lis[0]+'.'+lis[1] for lis in [stacod.split('.') for stacod in stacodes]]
    stacods = unique_list(stacods)
    stas = [stacod.split('.')[1] for stacod in stacods]

    stations, dist = list(zip(*sorted(zip(stas, dists))))

    #write stadistlist.txt
    f = open(moddir+'/stadistlist.txt', 'w')
    for sta, dis in zip(stations, dist):
        f.write(('%s\t%5.1f\n') % (sta, dis))
    f.close()

    #write dist file in free format
    #figure out how many samples
    samples = nextpow2(duration*samplerate)
    f = open(moddir+'/dist', 'w')
    for dis in dist:
        f.write(('%0.1f %0.2f %i %i 0\n') % (dis, 1./samplerate, samples, T0))
    f.close()

    #move copy of modelfile to current directory
    shutil.copy2(modelfile, moddir+'/'+modelfile.split('/')[-1])

    #write shell script to run Green's functions
    f = open(moddir+'/CPScommands.sh', 'w')
    f.write(("""#!/bin/bash
rm %s/*.sac
rm %s/*.sac
hprep96 -HR 0. -HS 0. -M %s -d %s -R -EXF
hspec96 > hspec96.out
hpulse96 -d %s -V -OD -p > Green
f96tosac Green
cp *.sac %s/.
mv *.sac %s/.
""") % (sacodir, sacdir, modelfile, moddir+'/dist', moddir+'/dist', sacdir, sacodir))
    os.chmod(moddir+'/CPScommands.sh', stat.S_IRWXU)

    return moddir


def compute_greens(event_id, filepath, shellscript='CPScommands.sh'):
    """
    This function actually runs the codes, and names the sac files properly
    INPUT
    filepath - file path to main folder for the current event to run
    shellscript - name of shell script to run
    """
    currentdir = os.getcwd()
    os.chdir(filepath)
    cmd = os.path.join(filepath, shellscript)
    proc = subprocess.Popen(cmd,
                            shell=True,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE
                            )
    stdout, stderr = proc.communicate()
    retcode = proc.returncode
    if retcode == 0:
        retcode = True
    else:
        retcode = False

    #load in stadistlist
    f = open(filepath+'/stadistlist.txt', 'r')
    lines = f.readlines()
    lines = [line.split('\t') for line in lines]
    temp2 = [(line[0], float(line[1].split('\n')[0])) for line in lines]
    sta, dist = list(zip(*temp2))

    #copy and rename files
    files = [os.path.basename(x) for x in glob.glob(filepath+'/sacdata/*.sac')]
    files.sort()
    for file1 in files:
        #get number of event
        indx = int(file1[1:4])-1
        GFt = file1[6:9]
        #rename
        newname = ('GF_ev%i_%s_%1.0fkm_%s.sac') % (event_id, sta[indx], dist[indx], GFt)
        os.rename(filepath+'/sacdata/'+file1, filepath+'/sacdata/'+newname)
    os.chdir(currentdir)
    return (retcode, stdout, str(stderr))
