import sys
import os
import glob
import sqlite3
from optparse import OptionParser

from matplotlib import use
use('Agg')
import pylab

from glue.lal import CacheEntry
from glue.ligolw import lsctables

from pylal import SnglInspiralUtils
from pylal import git_version


def parse_command_line():
    """
    Parse the command line, return an options object
    """
    parser = OptionParser(
        version     = "Name: %%prog\n%s" % git_version.verbose_msg,
        usage       = "%prog [--run-chisq|--assemble-results] --ifo ifo --basedir dir --gps-start-time start --gps-end-time end --snr-threshold snr",
        description = ""
    )

    parser.add_option("-r", "--run-chisq",   action = "store_true", help = "Construct trigger, sub, and cache files to calculate the chisq values of each trigger exceeding threshold")
    parser.add_option("-a", "--assemble-results", action = "store_true", help = "Create snr vs time plots for unclustered and 100-millisecond clustered triggers, and new snr vs time plot for clustered triggers.  Print a report in wiki-friendly format.")
    parser.add_option("-i", "--ifo",     metavar = "ifo")
    parser.add_option("-b", "--basedir", metavar = "dir", help = "base directory within which to create analysis directories and store results")
    parser.add_option("-u", "--baseurl", metavar = "url", help = "base url for images")
    parser.add_option("-s", "--gps-start-time", metavar = "start", help = "Earliest time to look for triggers")
    parser.add_option("-e", "--gps-end-time",   metavar = "end",   help = "Latest time to look for triggers")

    parser.add_option("-t", "--snr-threshold",  metavar = "threshold", help = "Lower limit of SNRs to examine")

    
    options, others = parser.parse_args()

    if not (options.run_chisq or options.assemble_results):
        raise ValueError, "Missing required argument, one of --run-chisq or --assemble-results"

    if (options.run_chisq and options.assemble_results):
        raise ValueError, "Please specify only one of --run-chisq or --assemble-results"

    if not options.gps_start_time:
        raise ValueError, "Missing required argument --gps-start-time"

    if not options.gps_end_time:
        raise ValueError, "Missing required argument --gps-end-time"

    if not options.snr_threshold:
        raise ValueError, "Missing required argument --snr-threshold"

    if not options.ifo:
        raise ValueError, "Missing required argument --ifo"
    
    return options



def times_to_dirs(start, end):
    sec_in_day = 60 * 60 * 24

    return ['%s/%s/' % (IHOPE_DAILY_DIR,tconvert(x, fmt='%Y%m/%Y%m%d')) for x in range(start, end, sec_in_day)]



def get_new_snr(trigger,index=6.0):
    rchisq = trigger.chisq/(2 * trigger.chisq_dof - 2)
    nhigh  = 2.

    if rchisq > 1.:
        return trigger.snr / ((1 + rchisq**(index/nhigh))/2)**(1./index)
    else:
        return trigger.snr


def plot_triggers(trig_file, center_time, center_time_ns, title, dest, new_snr = False):
    start_time  = center_time - 10
    end_time    = center_time + 10


    triggers = SnglInspiralUtils.ReadSnglInspiralFromFiles([trig_file])
    triggers = [t for t in triggers if t.end_time >= start_time and t.end_time < end_time]

    xs = [t.end_time + t.end_time_ns * 1.0e-9 - center_time for t in triggers]
    ys = new_snr and [get_new_snr(t) for t in triggers] or [t.snr for t in triggers]
    zs = [t.tau0 for t in triggers]

    pylab.figure()
    pylab.scatter(xs, ys, c=zs, edgecolor='none')

    if max(ys) > 10:
        pylab.yscale('log')

    pylab.ylim(min(ys), max(ys))

    pylab.xlabel('Time (s since %.3f)' % center_time)
    pylab.ylabel(new_snr and 'New SNR' or 'SNR')
    pylab.title(title)

    cb = pylab.colorbar()
    cb.ax.set_ylabel('Tau 0')

    pylab.savefig(dest)

    loudest_snr = max(ys)
    pos = 0
    for i in range(len(ys)):
        if ys[i] == loudest_snr:
            pos = i

    # Find the original trigger so we can report on it's chisq, etc
    orig = [t for t in triggers if t.end_time == center_time and t.end_time_ns == center_time_ns][0]

    return orig.chisq, get_new_snr(orig), triggers[pos].end_time, triggers[pos].end_time_ns, triggers[pos].snr


def build_command_line(dag_file, job_id):
    f    = open (dag_file,'r')
    job  = None
    args = None
    prms = {}

    for l in f:
        if l.startswith('#'):
            continue

        l    = l[:-1]
        tmp1 = l.split(' ')

        if len(tmp1) < 2:
            continue

        curr_job_id = tmp1[1]

        if job_id == curr_job_id:
            if l.startswith('JOB'):
                cmd_f = open(l[37:].strip(),'r')
                for l2 in cmd_f:
                    l2 = l2[:-1]
                    if l2.startswith('executable'):
                        job = l2[13:]
                    if l2.startswith('arguments'):
                        job += ' '
                        job += l2[12:]
                        job  = job.replace('$','%')
                        job  = job.replace(')',')s')
                        job  = job.replace('"','')

                cmd_f.close()
                
            if l.startswith('VARS'):
                args = l[38:].split('" ')
                for arg in args:
                    name,value = arg.split('=')
                    value      = value.replace('"','')
                    prms[name] = value
                    
            if args != None and job != None:
                return job % prms
    
def tconvert(tme, fmt=None):
    cmd = str(tme)

    if fmt:
        cmd = '-f ' + fmt + ' ' + cmd

    return os.popen('tconvert ' + cmd).next().strip()

IHOPE_DAILY_DIR='/archive/home/cbc/ihope_daily/'
        

def file_with_trigger(dir, ifo, clustering, t_end_time):
    for filename in glob.glob('%s/%s-INSPIRAL_%s-*xml.gz' % (dir, ifo, clustering)):
        chunks = filename.split('-')
        tme    = int(chunks[2])
        dur    = int(chunks[3][:-7])

        if t_end_time >= tme and t_end_time < (tme + dur):
            return filename

    return None


def make_chisq_dir(basedir, ifo, trigger_time):
    start = trigger_time - 10
    end   = trigger_time + 10

    orig_dir = os.getcwd()
    newdir   = basedir + '/' + str(trigger_time)
    os.makedirs(newdir)
    os.chdir(newdir)

    dir = IHOPE_DAILY_DIR + tconvert(start, fmt='%Y%m/%Y%m%d')

    for filename in glob.glob('%s/%s-INSPIRAL_UNCLUSTERED-*xml.gz' % (dir, ifo)):
        chunks = filename.split('-')
        tme    = int(chunks[2])
        dur    = int(chunks[3][:-7])

        if start >= tme and start < (tme + dur):
            # Find the job in the dag that makes the unclustered triggers for this time
            os.chdir(dir)

            timearg = 'macrotrigstarttime="%d"' % tme
            f = [l for l in open('%s/daily.dag' % dir) if l.find('UNCLUSTERED_00') > 0 and l.find(timearg) > 0]

            if f == []:
                print "No line in dag to create triggers at time %d" % tme
                sys.exit(0)

            jobid = f[0].split(' ')[1]
            job   = build_command_line('daily.dag',jobid)

            # Change the user tag to CHISQ
            job = job.replace('UNCLUSTERED_00','CHISQ')

            # Read LDAS-STRAIN instead of DMT-STRAIN, which may have moved to tape
            job = job.replace('DMT-','LDAS-')

            # Turn on chisq
            job = job.replace('chisq-bins 0','chisq-bins 16')

            # Look at the full template bank, not just the 00 subbank
            job = job.replace('_00','')

            # and look for it in the right place
            job = job.replace('bank-file ','bank-file %s/' % dir)

            # change the frame cache to the one we'll create
            pos1 = job.find('frame-cache ')
            pos2 = job.find(' ', pos1 + len('frame-cache '))
            job  = job[:pos1] + 'frame-cache ldas.cache ' + job[pos2:]
            
            # get rid of the executable name
            pos = job.find(' ')
            job = job[pos:]

            # Add the final argument
            job = job + ' --trigger-file triggers.xml.gz'

            os.chdir(orig_dir)
            os.chdir(newdir)

            # Create a submit file.
            newsub = open('chisq.sub','w')
            print >>newsub,"""Executable = /archive/home/cbc/opt/s6b/latest/bin/lalapps_populate_chisq
Args       = %s
getenv     = True
Universe   = local
output     = chisq.$(cluster).$(process).out
error      = chisq.$(cluster).$(process).error
Log        = chisq.$(cluster).$(process).log
Queue""" % job

            newsub.close()

            # Next, make a trigger file with just the triggers at the time of interest
            target = dbtables.get_connection_filename('triggers.db')
            connection = ligolw_sqlite.setup(target)

            ligolw_sqlite.insert(connection, ['%s/%s-INSPIRAL_100MILLISEC_CLUSTERED-%d-%d.xml.gz' % (dir, ifo, tme, dur)])

            cursor = connection.cursor()
            cursor.execute('delete from sngl_inspiral where end_time < %d or end_time > %d' % (start, end))

            ligolw_sqlite.extract(connection, 'triggers.xml.gz')

            os.remove('triggers.db')

            # Make a new cache, since the robots have all the original DMT data!
            os.system('ligo_data_find --observatory %s --url-type file --server ldr.ligo.caltech.edu --gps-start-time %d --gps-end-time %d --output ./ldas.cache --lal-cache --type %s_LDAS_C02_L2' % (ifo[0], (start - 2048), (start + 2048), ifo))

            os.chdir(orig_dir)


if __name__ == '__main__':
    options = parse_command_line()

    start_time = int(options.gps_start_time)
    end_time   = int(options.gps_end_time)
    threshold  = int(options.snr_threshold)
    ifo        = options.ifo
    basedir    = options.basedir

    dirs  = times_to_dirs(start_time, end_time)
    trigs = []

    for d in dirs:
        f_in = open('%s/%s-0-INSPIRAL_16SEC_CLUSTERED.csv' % (d, ifo))

        for l in f_in:
            trig = l.strip().split(',')

            if float(trig[3]) > threshold:
                trigs.append(  (d, int(trig[0]), int(trig[1]), float(trig[3]))  )


    if options.run_chisq:
        # Can't load these at the top, if we do then
        # calls to ReadSnglInspiral will fail
        from glue.ligolw import dbtables
        from glue.ligolw.utils import ligolw_sqlite

        # Remove the constraint that event_ids must be unique,
        # since first inspiral triggers all have event_id = 0
        del lsctables.SnglInspiralTable.constraints
        for dir, t_end_time, t_end_time_ns, snr in trigs:
            make_chisq_dir(basedir, ifo, t_end_time)

    if options.assemble_results:
        print '|| OL end time || OL end time ns || OL SNR || NL end time || NL end time ns || NL new SNR || Vetoed at cat ||',

        if options.baseurl:
            print '  Unclustered || 100 millisec clustered || New SNR ||'
        else:
            print

        cwd = os.getcwd()

        for dir, t_end_time, t_end_time_ns, snr in trigs:
            chisq_dir = '%s/%d' % (basedir, t_end_time)

            os.chdir(cwd)
            os.chdir(chisq_dir)

            unclustered_trig_file = file_with_trigger(dir, ifo, 'UNCLUSTERED', t_end_time)
            clustered_trig_file   = file_with_trigger(dir, ifo, '100MILLISEC_CLUSTERED', t_end_time)
            chisq_trig_file       = glob.glob('*CHISQ*')

            if chisq_trig_file == []:
                print >>sys.stderr,"Job in directory %s did not succeed" % chisq_dir
                continue

            chisq_trig_file = chisq_trig_file[0]

            plot_triggers(unclustered_trig_file, t_end_time, t_end_time_ns, "Unclustered", '%d_unclustered.png' % t_end_time)
            plot_triggers(clustered_trig_file, t_end_time, t_end_time_ns, "100 Millisec clustered", '%d_clustered.png' % t_end_time)
            original_chisq, original_new_snr, loudest_t, loudest_t_ns, loudest_new_snr = plot_triggers(chisq_trig_file, t_end_time, t_end_time_ns, "100 Millisec clustered, new snr", '%d_chisq.png' % t_end_time, new_snr = True)

            vetoed = 'NA'

            for cat in [1,2,4]:
                f_in  = open('%s/%s-%d-INSPIRAL_16SEC_CLUSTERED.csv' % (dir, ifo, cat))
                found = False

                for l in f_in:
                    trig = l.strip().split(',')
                    if int(trig[0]) == t_end_time and int(trig[1]) == t_end_time_ns:
                        found = True
                        break

                if not found:
                    vetoed = str(cat)
                    break
            
            print '|| %d.%d || %.2f || %.2f || %.2f || %d.%d || %.2f || %s ||' % (t_end_time, t_end_time_ns, snr, original_chisq, original_new_snr, loudest_t, loudest_t_ns, loudest_new_snr, vetoed),

            if options.baseurl:
                u = options.baseurl
                t = t_end_time

                print ' [[%s/%d_unclustered.png|img]] || [[%s/%d_clustered.png|img]] || [[%s/%d_chisq.png|img]]  ||' % (u, t, u, t, u, t)
            else:
                print



