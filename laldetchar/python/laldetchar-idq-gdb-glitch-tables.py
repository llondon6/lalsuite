# Copyright (C) 2014 Reed Essick, Ruslan Vaulin
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

description = \
    """The program  generates a summary of iDQ output during a short time period. Its output is written into a new xml file(s)"""

from ConfigParser import SafeConfigParser

from ligo.gracedb.rest import GraceDb
from laldetchar.idq import idq
from laldetchar.idq import idq_tables
from laldetchar.idq import idq_tables_dbutils

from glue.ligolw import ligolw
from glue.ligolw.utils import ligolw_add
from glue.ligolw import utils as ligolw_utils
from glue.ligolw import lsctables
from glue.ligolw import table
from glue.ligolw.utils import ligolw_sqlite
from glue import segments
from glue.ligolw import dbtables
from glue import lal
import sqlite3

import sys

from laldetchar import git_version

__author__ = 'Reed Essick <reed.essick@ligo.org>, Ruslan Vaulin <ruslan.vaulin@ligo.org>'
__version__ = git_version.id
__date__ = git_version.date

##################################################
# main
##################################################
from optparse import OptionParser
parser = OptionParser(version='Name: %%prog\n%s'% git_version.verbose_msg,
    usage='%prog [options]',
    description=description)

parser.add_option('-v',
    '--verbose',
    default=False,
    action='store_true')

parser.add_option(
    '-c', '--config',
    default='idq.ini',
    type='string',
    )

parser.add_option(
    '-s',
    '--gps-start',
    dest='start',
    default=0,
    type='float',
    help='the gps start time of the time range of interest',
    )

parser.add_option(
    '-e',
    '--gps-end',
    dest='end',
    default=0,
    type='float',
    help='the gps end time of the time range of interest',
    )

parser.add_option('-g',
    '--gracedb-id',
    default=None,
    type='string',
    help='GraceDB ID')

parser.add_option('',
    '--skip-gracedb-upload',
    default=False,
    action='store_true',
    help='skip steps involving communication with GraceDB. Automatically set to True if --gracedb-id==None')

#parser.add_option('',
#    '--ifo',
#    type='string',
#    help='the ifo for which predictions were made')

parser.add_option('-C',
    '--classifier',
    default='ovl',
    type='string',
    help='the classifier that was used to generate the data. Default="ovl"')

#parser.add_option('-i',
#    '--input-dir',
#    default='./',
#    type='string',
#    help='the directory through which is searched for relevant *glitch*.xml files. Assumes directory structure generated by laldetchar-idq-realtime.py'
#    )

#parser.add_option('-o',
#    '--output-dir',
#    default='.',
#    type='string',
#    help='the output directory')

#parser.add_option('-t',
#    '--usertag',
#    dest='tag',
#    default='',
#    type='string',
#    help='user tag')

#parser.add_option('',
#    '--gdb-url',
#    default=False,
#    type='string',
#    help='url of GraceDB e.g. ')


(opts, args) = parser.parse_args()

opts.skip_gracedb_upload = (opts.gracedb_id==None) or opts.skip_gracedb_upload

#=================================================
# read relevant stuff from config file
#=================================================
config = SafeConfigParser()
config.read( opts.config )

ifo = config.get('general','ifo')
tag = config.get('general','usertag')

realtimedir = config.get('general','realtimedir')
gdbdir = config.get('gdb general','main_gdb_dir')

if not opts.skip_gracedb_upload:
    if config.has_option('gdb general', 'gdb_url'):
        gracedb = GraceDb(config.get('gdb general', 'gdb_url'))
    else:
        gracedb = GraceDb()

##########################################
### Find relevant files
##########################################

if opts.verbose:
    print 'Finding relevant *glitch*.xml files'
gchxml_filenames = sorted([filename for filename in 
                          idq.get_all_files_in_range(realtimedir, opts.start, opts.end, pad=0, suffix='.xml.gz') 
                          if opts.classifier == idq.extract_xml_name(filename) 
                          and 'glitch' in filename
                          and ifo in filename])
                        
if not gchxml_filenames:
    # no files found, print the message, and exit
    if not opts.skip_gracedb_upload:
        gracedb.writeLog(opts.gracedb_id, message="No iDQ glitch tables data from "+opts.classifier+" available for the candidate  at "+ifo)
    print "No glitch files found, exiting."
    sys.exit(0)

if opts.verbose:
    print "Found:"
    for filename in gchxml_filenames:
        print '\t' + filename

# ##############################################################
# ## Load and Merge xml files using in-memory sqlite database
# ##############################################################
# 
# load files into database
connection, cursor = idq_tables_dbutils.load_xml_files_into_database(\
    gchxml_filenames, verbose=opts.verbose)

# #########################################################
# ## remove redundant rows and any events outside of [opts.start, opts.end]
###########################################################

# remove redundant entries from the tables
idq_tables_dbutils.remove_redundant_entries(connection, cursor, verbose = opts.verbose)

# form two open segments using start and stop times
seglist = segments.segmentlist()
seglist.append(segments.segment([-segments.infinity(), lal.LIGOTimeGPS(opts.start)]))
seglist.append(segments.segment([lal.LIGOTimeGPS(opts.end), segments.infinity()]))

# delete glitch events that fall inside of these segments
idq_tables_dbutils.delete_glitch_events_in_segmentlist(connection, cursor, seglist)


###############################################################################
# ## save merged xmldoc
###############################################################################
#merged_xmldoc_filename = '%s/%s_idq_%s_glitch_%s%d-%d.xml' % (
#    opts.output_dir,
#    opts.ifo,
#    opts.classifier,
#    opts.tag,
#    int(opts.start),
#    int(opts.end - opts.start)
#    )
merged_xmldoc_filename = idq.gdb_xml(gdbdir, opts.classifier, ifo, tag, int(opts.start), int(opts.end-opts.start))

if opts.verbose:
    print 'saving ' + merged_xmldoc_filename
# exctract data base into xml file and write it to disk
ligolw_sqlite.extract(connection, merged_xmldoc_filename , verbose = opts.verbose) # exctract data base into xml file
connection.close()
if not opts.skip_gracedb_upload:
    #write log message to gracedb and upload file
    gracedb.writeLog(opts.gracedb_id, message="iDQ glitch tables " + ifo + ":", filename=merged_xmldoc_filename)
