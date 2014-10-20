#-----------------------------------------------------------------------bl-
#--------------------------------------------------------------------------
#
# QUESO - a library to support the Quantification of Uncertainty
# for Estimation, Simulation and Optimization
#
# Copyright (C) 2008,2009,2010,2011,2012,2013 The PECOS Development Team
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the Version 2.1 GNU Lesser General
# Public License as published by the Free Software Foundation.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc. 51 Franklin Street, Fifth Floor,
# Boston, MA  02110-1301  USA
#
#-----------------------------------------------------------------------el-

import re
import sys

def convertOne(mFileName, getChain=True, getLike=True, getTarget=True, chainIdent=r'[cC]hain\w*', logLikeIdent=r'[lL]ogLikelihood\w*', logTargetIdent=r'[lL]ogTarget\w*'):
    mFileNameBase = mFileName.rstrip('.m')
    dataEncloseTokens = ('[', ']')

    fIn = open(mFileName,'r')
    fileBuf = fIn.read()
    fIn.close()

    parseMap = {'chain': (getChain, chainIdent), 'logLikelihood': (getLike, logLikeIdent), 'logTarget': (getTarget, logTargetIdent)}
    for item in parseMap:
        if parseMap[item][0] == False:
            continue
        
        # Find tag for current item via search with regex: identString\s*=\s*(?=\[)
        parseTagRegEx = parseMap[item][1] + r'\s*=\s*' + '(?=' + re.escape(dataEncloseTokens[0]) + ')'
        match = re.search(parseTagRegEx, fileBuf)
        tag = ''
        if match:
            tag = match.group(0) + re.escape(dataEncloseTokens[0])
        else:
            print 'Error: Could not find identifier {0} in {1}. Regex string: {2}'.format(item, mFileName, parseTagRegEx)
            return -1

        # Extract data between '[ ]' for current item via search with regex: (?<=tag)[^\]]+
        parseDataRegEx = '(?<=' + tag + ')[^' + re.escape(dataEncloseTokens[1]) + ']+'
        match = re.search(parseDataRegEx, fileBuf)
        if match:
            outName = mFileNameBase + '_' + item + '.txt'
            fOut = open(outName, 'w')
            fOut.write(match.group(0))
            fOut.close()
        else:
            print 'Error: Could not find {0} data in {1}. Regex string: {2}'.format(item, mFileName, parseDataRegEx)
            return -1


if __name__ == "__main__":
    usage = 'Usage: {0} [-chain <yes|no>] [-like <yes|no>] [-target <yes|no>] [-chainIdent <string>] [-likeIdent <string>] [-targetIdent <string>] <file list>'.format(sys.argv[0])
    if len(sys.argv) < 2 or '--help' in sys.argv or '-h' in sys.argv:
        print usage
        sys.exit(1)

    # parse options first
    optMapGets = {'-chain': True, '-like': True, '-target': True}
    for opt in optMapGets:
        if opt in sys.argv:
            ind = sys.argv.index(opt)
            if ind+1 >= len(sys.argv) or sys.argv.count(opt) > 1:
                print usage
                sys.exit(1)
            val = sys.argv.pop(ind+1)
            if val.lower() == 'yes':
                optMapGets[opt] = True
            elif val.lower() == 'no':
                optMapGets[opt] = False
            else:
                print 'Unrecognized value for option {0}: {1}'.format(opt, val)
                sys.exit(1)
            sys.argv.pop(ind)

    optMapIdents = {'-chainIdent': r'[cC]hain\w*', '-likeIdent': r'[lL]ogLikelihood\w*', '-targetIdent': r'[lL]ogTarget\w*'}
    for opt in optMapIdents:
        if opt in sys.argv:
            ind = sys.argv.index(opt)
            if ind+1 >= len(sys.argv) or sys.argv.count(opt) > 1:
                print usage
                sys.exit(1)
            optMapIdents[opt] = sys.argv.pop(ind+1)
            sys.argv.pop(ind)

    # parse files one-by-one
    for mFileName in sys.argv[1:]:
        print 'Converting {0}...'.format(mFileName)
        convertOne(mFileName, getChain=optMapGets['-chain'], getLike=optMapGets['-like'], getTarget=optMapGets['-target'], chainIdent=optMapIdents['-chainIdent'], logLikeIdent=optMapIdents['-likeIdent'], logTargetIdent=optMapIdents['-targetIdent'])
