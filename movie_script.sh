#!/bin/bash
# Parameters
r=0.00063         # Coords of zoom point
i=0.06156
total_zoom=0.0001 # Zoom by fractor of 1000
frames=100        # No. of frames in sequence

zoom_per_frame=`echo "scale=40; e(l($total_zoom)/$frames)" | bc -l`
frame_no=0

while true
do
    rad=`echo "scale=20; e(l($zoom_per_frame)*$frame_no)" | bc -l`
    echo ============================== frame = $frame_no ==============================
    echo running ./main $r $i $rad
    ./main $r $i $rad
    cp viewer/public/output/fractal.pbm $frame_no.ppm
    if test $frame_no == $frames
    then
        break
    fi
    frame_no=$(( frame_no + 1 ))
done

convert -delay 15 ?.ppm ??.ppm motion.gif
