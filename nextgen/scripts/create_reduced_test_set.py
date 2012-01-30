"""
Creates a minimal raw sequencing data test set from a flowcell directory by copying a
subset of the tiles and just the files needed for downstream analysis, as well as 
modifying the relevant config files to reflect this.
    
Usage:
    create_reduced_test_set.py <flow cell dir> [options]
    
Options:
    
    -s, --subset=<FLOAT | INT>             Specify the size of the subset of tiles to use. 
                                           0 < FLOAT < 1 or 1 <= INT <= NUM_TILES. If < 1, 
                                           the specified fraction of tiles (rounded down) 
                                           will be used. If >= 1, the specified number of 
                                           tiles will be used. 
                                  
    -t, --target_dir=<target directory>    The destination directory where the minimal test
                                           set will be written to a sub-directory with the
                                           same name as the source directory. There must not
                                           already exist a folder for the flowcell within the 
                                           target directory.
                    

"""

import sys
import glob
import os
import re
import math
import xml.etree.ElementTree as xml
from random import randint
from shutil import (copy2,copytree,rmtree)

def main(run_dir,subset_size,target_dir):
    
    # Get the name of the run folder
    _,run_name = os.path.split(os.path.normpath(run_dir))

    # Change to the target directory if specified
    cwd = os.getcwd()
    if target_dir is not None:
        target_dir = os.path.abspath(target_dir)
        assert os.path.exists(target_dir)
        os.chdir(target_dir)
        
    # Create the folder hierarchy to the intensity files
    assert not os.path.exists(run_name)
    os.mkdir(run_name,0770)
    os.chdir(run_name)

    # Copy the data files
    copy_data_files(run_dir)
    # Copy the meta files
    print "Copying meta files"
    copy_meta_files(run_dir)
    
    # Change back to the original cwd
    os.chdir(cwd)
    
def copy_data_files(run_dir,subset_size):

    # Copy and parse the config.xml file to get the available tiles common to all lanes
    config_file = os.path.join(run_dir,'Data','Intensities','config.xml')
    all = _parse_config(config_file)
    
    # Randomly pick a subset of tiles to copy
    subset = []
    if subset_size < 1.0:
        desired = max(1,int(math.floor(subset_size*len(all))))
    else:
        desired = min(math.floor(subset_size),len(all))
    
    while len(subset) < desired:
        subset.append(all.pop(randint(1,len(all))-1))

    print "%s tiles available, will copy a subset of %s tiles" % (len(all),len(subset))
    
    to_copy = []
    size = 0
    idir = os.path.join(run_dir,'Data','Intensities')
    bcdir = os.path.join(run_dir,'Data','Intensities','BaseCalls')
    for tile in subset:
        globs = [os.path.join(idir,"*%s*" % tile),
                 os.path.join(idir,"L*","*%s*" % tile),
                 os.path.join(idir,"L*","C*","*%s*" % tile),
                 os.path.join(idir,"Offsets","*%s*" % tile),
                 os.path.join(bcdir,"*%s*" % tile),
                 os.path.join(bcdir,"L*","*%s*" % tile),
                 os.path.join(bcdir,"L*","C*","*%s*" % tile),
                 os.path.join(bcdir,"Matrix", "*%s*" % tile),
                 os.path.join(bcdir,"Phasing", "*%s*" % tile)]
        for g in globs:
            for item in glob.glob(g):
                to_copy.append([item,os.path.relpath(os.path.dirname(item),run_dir)])
                size += os.path.getsize(to_copy[-1][0])
             
    print "Will copy %s files, totalling %s Mbytes" % (len(to_copy),size/(1024*1024))
    
    # Copy the files, creating subdirectories as necessary
    for src_file, dst_dir in to_copy:
        if not os.path.exists(dst_dir):
            print "Creating directory %s" % dst_dir
            os.makedirs(dst_dir,0770)
        print "Copying file %s to %s" % (src_file,dst_dir)
        copy2(src_file,dst_dir)
    
    config_file_cpy = os.path.join('Data','Intensities',os.path.basename(config_file))
    print "Copying and updating config file %s" % config_file
    _update_config(config_file,config_file_cpy,subset)

    config_file = os.path.join(run_dir,'Data','Intensities','BaseCalls','config.xml')
    config_file_cpy = os.path.join('Data','Intensities','BaseCalls',os.path.basename(config_file))
    print "Copying and updating config file %s" % config_file
    _update_config(config_file,config_file_cpy,subset)

def copy_meta_files(run_dir):
    
    bcdir = os.path.join('Data','Intensities','BaseCalls')
    globs = ['*.xml',
             '*.csv',
             '*.txt',
             os.path.join(bcdir,'*.xsl'),
             os.path.join(bcdir,'*.htm'),
             os.path.join(bcdir,'Plots'),
             os.path.join('Data','reports'),
             os.path.join('Data','Status_Files'),
             os.path.join('Data','Status.htm'),
             'InterOp']
    
    for g in globs:
        path = os.path.join(run_dir,g)
        print "Copying %s..." % path
        for item in glob.glob(path):
            relpath = os.path.relpath(item,run_dir)
            if os.path.isdir(item):
                # Remove the directory if it already exists
                if os.path.exists(relpath):
                    rmtree(relpath)
                copytree(item,relpath)
            elif os.path.isfile(item):
                copy2(item,relpath)
    
def _file_in_subset(file,subset):
    return _file_tile(file) in subset

def _file_tile(file):
    m = re.search(r"_(\d{4})",file)
    if m is not None and len(m.groups()):
        return int(m.group(1))
    return -1

def _parse_config(config_file):
    all_tiles = []
    cfg = xml.parse(config_file)
    root = cfg.getroot()
    tiles = root.find('Run/TileSelection')
    if tiles is None:
        return all_tiles
    for lane in tiles.findall('Lane'):
        lane_tiles = []
        for tile in lane.findall('Tile'):
            lane_tiles.append(int(tile.text))
        if len(all_tiles) == 0:
            all_tiles = lane_tiles
        else:
            all_tiles = list(set(all_tiles).intersection(lane_tiles))
    return all_tiles

def _update_config(config_file_src,config_file_dst,subset):
    cfg = xml.parse(config_file_src)
    root = cfg.getroot()
    tiles = root.find('Run/TileSelection')
    if tiles is None:
        return
    for lane in tiles.findall('Lane'):
        for tile in lane.findall('Tile'):
            if int(tile.text) not in subset:
                lane.remove(tile)
    cfg.write(config_file_dst,encoding='UTF8')

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-s", "--subset", dest="subset", default=0.05)
    parser.add_option("-t", "--target_dir", dest="target_dir", default=None)
    parser.add_option("-v", "--verbose", dest="verbose", default=False, \
                                                        action="store_true")
    options, args = parser.parse_args()
    if len(args) == 1:
        run_dir, = args
    else:
        print __doc__
        sys.exit()
    main(os.path.abspath(run_dir),float(options.subset),options.target_dir)
    