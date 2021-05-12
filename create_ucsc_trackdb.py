#! /usr/bin/env python

str_main_track = '''track %s
superTrack on show
shortLabel %s
longLabel %s
'''

str_bw_track = '''
track %s
type bigWig
description %s
bigDataUrl %s
visibility 2
color 40,60,80
smoothingWindow OFF
maxHeightPixels 100:20:8
autoScale on
priority %s
graphTypeDefault bar
'''

def main():
    import os, sys

    outfile = sys.argv[1]
    ucsc_folder = sys.argv[2]
    bigwigs = sys.argv[3:]

    # translate passed folder for storing publicly available bigwigs into url
    big_data_url = ucsc_folder.replace('/project/MDL_ChIPseq', 'http://owww.molgen.mpg.de/~MDL_ChIPseq').replace('/project/ngsuploads/data', 'http://ngs.molgen.mpg.de')
    
    # sort file list by species and stage
    with open(outfile, 'w') as f:
        track_name = os.path.basename(os.path.dirname(outfile))
        f.write(str_main_track %(track_name, track_name, track_name))
        for i, bw in enumerate(bigwigs):
            bw_url = big_data_url + os.path.basename(bw)
            sample_name = os.path.basename(bw).split('.')[0]
            f.write(str_bw_track %(sample_name, sample_name, bw_url, i))
    
if __name__ == '__main__':
    main()
