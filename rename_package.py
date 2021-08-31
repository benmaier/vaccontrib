import sys
import pathlib

include_fileendings = ["py", "txt", "yml", "cfg", "rst", "md", "ini", "in"]
include_files = ["Makefile"]
exclude_files = ["rename_package.py"]

args = sys.argv[1:]
new_name = args[0]

this_dir = pathlib.Path('.').absolute()

for drc in this_dir.glob("**/*"):
    if drc.is_dir():
        new_dir = str(drc.absolute()).replace("bfmdummypackage",new_name)
        drc.rename(new_dir)

for file in this_dir.glob("**/*"):
    if file.is_file():
        new_file = str(file.absolute()).replace("bfmdummypackage",new_name)
        file.rename(new_file)

for file in this_dir.glob("**/*"):
    if file.is_file() and \
       ((file.suffix[1:] in include_fileendings) or (file.name in include_files)) and \
       file.name not in exclude_files:
        with open(file,'r') as f:
            text = f.read()
            new_text = text.replace("bfmdummypackage",new_name)
        with open(file,'w') as f:
            f.write(new_text)

