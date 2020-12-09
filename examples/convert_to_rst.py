#!/usr/bin/env python

# COnvert all notebooks to .rst format and save them
# in the ../doc folder

from __future__ import print_function
import glob
import subprocess
import shutil
import os

notebooks = glob.glob("*.ipynb")

for nb in notebooks:
    
    root = nb.split(".")[0]
    
    cmd_line = 'ipython nbconvert --to rst %s' % (nb)
    
    print(cmd_line)
    
    subprocess.check_call(cmd_line,shell=True)
    
    # Now move the .rst file and the directory with the data
    # under ../doc
    
    try:
        
        os.remove("../doc/%s.rst" % root)
    
    except:
        
        pass
    
    files_dir = "%s_files" % root
    
    try:
        
        shutil.rmtree("../doc/%s" % files_dir)
    
    except:
        
        pass
    
    shutil.move("%s.rst" % root,"../doc")
    
    if os.path.exists(files_dir):
    
        shutil.move(files_dir, "../doc")
