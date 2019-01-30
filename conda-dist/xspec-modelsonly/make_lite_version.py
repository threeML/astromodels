import sys
import subprocess
import json
import os
import platform

archive = sys.argv[1]

# Untar archive
system = platform.system()
workdir = "__workdir_%s" % system
os.makedirs(system)
os.makedirs(workdir)

os.chdir(workdir)

subprocess.check_call("tar xvf %s" % archive, shell=True)

# Find model files and their sizes

file_names = [os.path.join(path, name) for path, _, filenames in os.walk("lib") for name in filenames]

file_sizes = [(name, os.path.getsize(name) / 1000.0 / 1000.0) for name in file_names]

# Select all model files above 1 Mb
large_files = filter(lambda x:x[1] > 1.0 and os.path.splitext(x[0])[1]=='.fits', file_sizes)

# Delete them
deleted = []
for (filename, size) in large_files:
    
    deleted.append(filename)
    
    os.remove(filename)

# Now update the files list
files = list(map(lambda x:x.replace("\n",""), open("info/files").readlines()))

for deleted_file in deleted:
    
    files.pop(files.index(deleted_file))

with open("info/files", "w+") as f:
    
    for filename in files:
        
        f.write("%s\n" % filename)

# Update also the paths.json file
paths = json.load(open("info/paths.json"))

new_paths = []

for filepath in paths['paths']:
    
    if filepath['_path'] not in deleted:
        
        new_paths.append(filepath)

paths['paths'] = new_paths

json.dump(paths, open("info/paths.json", "w+"))

# Now update the index
index = json.load(open("info/index.json"))
index["name"] = "xspec-modelsonly-lite"
json.dump(index, open("info/index.json", "w+"))

subprocess.check_call("tar c lib info | bzip2 > ../%s/%s-%s-0.tar.bz2" % (system, index["name"], index["version"]), shell=True)

os.chdir("..")

subprocess.check_call("rm -rf %s" % workdir, shell=True)
