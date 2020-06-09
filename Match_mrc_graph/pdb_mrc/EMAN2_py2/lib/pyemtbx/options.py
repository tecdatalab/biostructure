#!/home/luis98/Documents/git-proyects/biostructure/Match_mrc_graph/pdb_mrc/EMAN2_py2/extlib/bin/python
#
# Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
# Copyright (c) 2000-2006 Baylor College of Medicine
#
# This software is issued under a joint BSD/GNU license. You may use the
# source code in this file under either license. However, note that the
# complete EMAN2 and SPARX software packages have some GPL dependencies,
# so you are responsible for compliance with the licenses of these packages
# if you opt to use BSD licensing. The warranty disclaimer below holds
# in either instance.
#
# This complete copyright notice must be included in any revised version of the
# source code. Additional authorship citations may be added, but existing
# author citations must be preserved.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  2111-1307 USA
#
#

def parse_filter_params(filterparams):
    params = filterparams.split(":")
    filtername = params[0]

    if len(params) == 1:
        return (filtername, None)
    else:
        d = Dict()
        for param in params[1:]:
            key_values = param.split("=")
            d[key_values[0]] = EMObject(key_values[1])
        return (filtername, d)


def get_optionlist(argv):
    optionlist = []
    for arg1 in argv:
        if arg1[0] == "-":
            argname = arg1.split("=")
            optionlist.append(argname[0].lstrip("-"))
    return optionlist

def intvararg_callback(option, opt_str, value, parser):
#    print 'vararg_callback:'
#    print '\toption:', repr(option)
#    print '\topt_str:', opt_str
#    print '\tvalue:', value
#    print '\tparser:', parser
    
    v = [int(i) for i in value.split(',')]
    setattr(parser.values, option.dest, v)
    return

def floatvararg_callback(option, opt_str, value, parser):
#    print 'vararg_callback:'
#    print '\toption:', repr(option)
#    print '\topt_str:', opt_str
#    print '\tvalue:', value
#    print '\tparser:', parser
    
    v = [float(i) for i in value.split(',')]
    setattr(parser.values, option.dest, v)
    return
