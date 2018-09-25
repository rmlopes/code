import os
import sys
notice = '''
#############################################################################
#  CoDe - A library for evolutionary computational devices                  #
#  Copyright (C) 2014 Rui L. Lopes                                          #
#                                                                           #
#  This program is free software: you can redistribute it and/or modify     #
#  it under the terms of the GNU General Public License as published by     #
#  the Free Software Foundation, either version 3 of the License, or        #
#  (at your option) any later version.                                      #
#                                                                           #
#  This program is distributed in the hope that it will be useful,          #
#  but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#  GNU General Public License for more details.                             #
#                                                                           #
#  A copy of the GNU General Public License is provided in the file         #
#  license.txt.  For more information, see <http://www.gnu.org/licenses/>.  #
#############################################################################
\n'''
targetdir = sys.argv[1]
dryrun=False
if len(sys.argv) > 2:
    dryrun = True


fileList = []
for root, subFolders, files in os.walk(targetdir):
    for f in files:
        tokens = f.split(".")
        if len(tokens) > 1 and tokens[-1] =='py':
            fileList.append(os.path.join(root,f))

if dryrun:
    for f in fileList: print f

if not dryrun:
    for f in fileList:
        fh = open(f,"r")
        lines = fh.readlines()
        fh.close()
        fh = open(f,"w")
        fh.write(notice)
        for l  in lines:
            fh.write(l)
        fh.close()
