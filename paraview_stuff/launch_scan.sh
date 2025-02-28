#!/bin/bash

pvpython_p2 --mesa scan_slice_cgpt.py 

#to make the video:
#ffmpeg -framerate 15 -i slice_0%03d.png -c:v libx264 -r 15 -pix_fmt yuv420p output.mp4
