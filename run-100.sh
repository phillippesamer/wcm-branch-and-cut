#!/bin/bash

instances_gnp_100_only50=(
"./input/alternate_g_np_100/100-1-0.gcc"  
"./input/alternate_g_np_100/100-1-1.gcc"  
"./input/alternate_g_np_100/100-1-2.gcc"  
"./input/alternate_g_np_100/100-1-3.gcc"  
"./input/alternate_g_np_100/100-1-4.gcc"  
"./input/alternate_g_np_100/100-2-0.gcc"  
"./input/alternate_g_np_100/100-2-1.gcc"  
"./input/alternate_g_np_100/100-2-2.gcc"  
"./input/alternate_g_np_100/100-2-3.gcc"  
"./input/alternate_g_np_100/100-2-4.gcc"  
"./input/alternate_g_np_100/100-3-0.gcc"  
"./input/alternate_g_np_100/100-3-1.gcc"  
"./input/alternate_g_np_100/100-3-2.gcc"  
"./input/alternate_g_np_100/100-3-3.gcc"  
"./input/alternate_g_np_100/100-3-4.gcc"  
"./input/alternate_g_np_100/100-4-0.gcc"  
"./input/alternate_g_np_100/100-4-1.gcc"  
"./input/alternate_g_np_100/100-4-2.gcc"  
"./input/alternate_g_np_100/100-4-3.gcc"  
"./input/alternate_g_np_100/100-4-4.gcc"  
"./input/alternate_g_np_100/100-5-0.gcc"  
"./input/alternate_g_np_100/100-5-1.gcc"  
"./input/alternate_g_np_100/100-5-2.gcc"  
"./input/alternate_g_np_100/100-5-3.gcc"  
"./input/alternate_g_np_100/100-5-4.gcc"  
"./input/alternate_g_np_100/100-6-0.gcc"  
"./input/alternate_g_np_100/100-6-1.gcc"  
"./input/alternate_g_np_100/100-6-2.gcc"  
"./input/alternate_g_np_100/100-6-3.gcc"  
"./input/alternate_g_np_100/100-6-4.gcc"  
"./input/alternate_g_np_100/100-7-0.gcc"  
"./input/alternate_g_np_100/100-7-1.gcc"  
"./input/alternate_g_np_100/100-7-2.gcc"  
"./input/alternate_g_np_100/100-7-3.gcc"  
"./input/alternate_g_np_100/100-7-4.gcc"  
"./input/alternate_g_np_100/100-8-0.gcc"  
"./input/alternate_g_np_100/100-8-1.gcc"  
"./input/alternate_g_np_100/100-8-2.gcc"  
"./input/alternate_g_np_100/100-8-3.gcc"  
"./input/alternate_g_np_100/100-8-4.gcc"  
"./input/alternate_g_np_100/100-9-0.gcc"  
"./input/alternate_g_np_100/100-9-1.gcc"  
"./input/alternate_g_np_100/100-9-2.gcc"  
"./input/alternate_g_np_100/100-9-3.gcc"  
"./input/alternate_g_np_100/100-9-4.gcc"  
"./input/alternate_g_np_100/100-10-0.gcc"  
"./input/alternate_g_np_100/100-10-1.gcc"  
"./input/alternate_g_np_100/100-10-2.gcc"  
"./input/alternate_g_np_100/100-10-3.gcc"  
"./input/alternate_g_np_100/100-10-4.gcc"  
)

instances_gnp_100=( 
"./input/alternate_g_np_100/100-1-0.gcc"  
"./input/alternate_g_np_100/100-1-1.gcc"  
"./input/alternate_g_np_100/100-1-2.gcc"  
"./input/alternate_g_np_100/100-1-3.gcc"  
"./input/alternate_g_np_100/100-1-4.gcc"  
"./input/alternate_g_np_100/100-2-0.gcc"  
"./input/alternate_g_np_100/100-2-1.gcc"  
"./input/alternate_g_np_100/100-2-2.gcc"  
"./input/alternate_g_np_100/100-2-3.gcc"  
"./input/alternate_g_np_100/100-2-4.gcc"  
"./input/alternate_g_np_100/100-3-0.gcc"  
"./input/alternate_g_np_100/100-3-1.gcc"  
"./input/alternate_g_np_100/100-3-2.gcc"  
"./input/alternate_g_np_100/100-3-3.gcc"  
"./input/alternate_g_np_100/100-3-4.gcc"  
"./input/alternate_g_np_100/100-4-0.gcc"  
"./input/alternate_g_np_100/100-4-1.gcc"  
"./input/alternate_g_np_100/100-4-2.gcc"  
"./input/alternate_g_np_100/100-4-3.gcc"  
"./input/alternate_g_np_100/100-4-4.gcc"  
"./input/alternate_g_np_100/100-5-0.gcc"  
"./input/alternate_g_np_100/100-5-1.gcc"  
"./input/alternate_g_np_100/100-5-2.gcc"  
"./input/alternate_g_np_100/100-5-3.gcc"  
"./input/alternate_g_np_100/100-5-4.gcc"  
"./input/alternate_g_np_100/100-6-0.gcc"  
"./input/alternate_g_np_100/100-6-1.gcc"  
"./input/alternate_g_np_100/100-6-2.gcc"  
"./input/alternate_g_np_100/100-6-3.gcc"  
"./input/alternate_g_np_100/100-6-4.gcc"  
"./input/alternate_g_np_100/100-7-0.gcc"  
"./input/alternate_g_np_100/100-7-1.gcc"  
"./input/alternate_g_np_100/100-7-2.gcc"  
"./input/alternate_g_np_100/100-7-3.gcc"  
"./input/alternate_g_np_100/100-7-4.gcc"  
"./input/alternate_g_np_100/100-8-0.gcc"  
"./input/alternate_g_np_100/100-8-1.gcc"  
"./input/alternate_g_np_100/100-8-2.gcc"  
"./input/alternate_g_np_100/100-8-3.gcc"  
"./input/alternate_g_np_100/100-8-4.gcc"  
"./input/alternate_g_np_100/100-9-0.gcc"  
"./input/alternate_g_np_100/100-9-1.gcc"  
"./input/alternate_g_np_100/100-9-2.gcc"  
"./input/alternate_g_np_100/100-9-3.gcc"  
"./input/alternate_g_np_100/100-9-4.gcc"  
"./input/alternate_g_np_100/100-10-0.gcc"  
"./input/alternate_g_np_100/100-10-1.gcc"  
"./input/alternate_g_np_100/100-10-2.gcc"  
"./input/alternate_g_np_100/100-10-3.gcc"  
"./input/alternate_g_np_100/100-10-4.gcc"  
"./input/alternate_g_np_100/100-11-0.gcc"  
"./input/alternate_g_np_100/100-11-1.gcc"  
"./input/alternate_g_np_100/100-11-2.gcc"  
"./input/alternate_g_np_100/100-11-3.gcc"  
"./input/alternate_g_np_100/100-11-4.gcc"  
"./input/alternate_g_np_100/100-12-0.gcc"  
"./input/alternate_g_np_100/100-12-1.gcc"  
"./input/alternate_g_np_100/100-12-2.gcc"  
"./input/alternate_g_np_100/100-12-3.gcc"  
"./input/alternate_g_np_100/100-12-4.gcc"  
"./input/alternate_g_np_100/100-13-0.gcc"  
"./input/alternate_g_np_100/100-13-1.gcc"  
"./input/alternate_g_np_100/100-13-2.gcc"  
"./input/alternate_g_np_100/100-13-3.gcc"  
"./input/alternate_g_np_100/100-13-4.gcc"  
"./input/alternate_g_np_100/100-14-0.gcc"  
"./input/alternate_g_np_100/100-14-1.gcc"  
"./input/alternate_g_np_100/100-14-2.gcc"  
"./input/alternate_g_np_100/100-14-3.gcc"  
"./input/alternate_g_np_100/100-14-4.gcc"  
"./input/alternate_g_np_100/100-15-0.gcc"  
"./input/alternate_g_np_100/100-15-1.gcc"  
"./input/alternate_g_np_100/100-15-2.gcc"  
"./input/alternate_g_np_100/100-15-3.gcc"  
"./input/alternate_g_np_100/100-15-4.gcc"  
"./input/alternate_g_np_100/100-16-0.gcc"  
"./input/alternate_g_np_100/100-16-1.gcc"  
"./input/alternate_g_np_100/100-16-2.gcc"  
"./input/alternate_g_np_100/100-16-3.gcc"  
"./input/alternate_g_np_100/100-16-4.gcc"  
"./input/alternate_g_np_100/100-17-0.gcc"  
"./input/alternate_g_np_100/100-17-1.gcc"  
"./input/alternate_g_np_100/100-17-2.gcc"  
"./input/alternate_g_np_100/100-17-3.gcc"  
"./input/alternate_g_np_100/100-17-4.gcc"  
"./input/alternate_g_np_100/100-18-0.gcc"  
"./input/alternate_g_np_100/100-18-1.gcc"  
"./input/alternate_g_np_100/100-18-2.gcc"  
"./input/alternate_g_np_100/100-18-3.gcc"  
"./input/alternate_g_np_100/100-18-4.gcc"  
"./input/alternate_g_np_100/100-19-0.gcc"  
"./input/alternate_g_np_100/100-19-1.gcc"  
"./input/alternate_g_np_100/100-19-2.gcc"  
"./input/alternate_g_np_100/100-19-3.gcc"  
"./input/alternate_g_np_100/100-19-4.gcc"  
"./input/alternate_g_np_100/100-20-0.gcc"  
"./input/alternate_g_np_100/100-20-1.gcc"  
"./input/alternate_g_np_100/100-20-2.gcc"  
"./input/alternate_g_np_100/100-20-3.gcc"  
"./input/alternate_g_np_100/100-20-4.gcc"  
"./input/alternate_g_np_100/100-21-0.gcc"  
"./input/alternate_g_np_100/100-21-1.gcc"  
"./input/alternate_g_np_100/100-21-2.gcc"  
"./input/alternate_g_np_100/100-21-3.gcc"  
"./input/alternate_g_np_100/100-21-4.gcc"  
"./input/alternate_g_np_100/100-22-0.gcc"  
"./input/alternate_g_np_100/100-22-1.gcc"  
"./input/alternate_g_np_100/100-22-2.gcc"  
"./input/alternate_g_np_100/100-22-3.gcc"  
"./input/alternate_g_np_100/100-22-4.gcc"  
"./input/alternate_g_np_100/100-23-0.gcc"  
"./input/alternate_g_np_100/100-23-1.gcc"  
"./input/alternate_g_np_100/100-23-2.gcc"  
"./input/alternate_g_np_100/100-23-3.gcc"  
"./input/alternate_g_np_100/100-23-4.gcc"  
"./input/alternate_g_np_100/100-24-0.gcc"  
"./input/alternate_g_np_100/100-24-1.gcc"  
"./input/alternate_g_np_100/100-24-2.gcc"  
"./input/alternate_g_np_100/100-24-3.gcc"  
"./input/alternate_g_np_100/100-24-4.gcc"  
"./input/alternate_g_np_100/100-25-0.gcc"  
"./input/alternate_g_np_100/100-25-1.gcc"  
"./input/alternate_g_np_100/100-25-2.gcc"  
"./input/alternate_g_np_100/100-25-3.gcc"  
"./input/alternate_g_np_100/100-25-4.gcc"  
)

instances_gnp_100_num_colours=( 
"7"  
"7"  
"7"  
"7"  
"7"  
"2"  
"2"  
"2"  
"2"  
"2"  
"10" 
"10" 
"10" 
"10" 
"10" 
"9"  
"9"  
"9"  
"9"  
"9"  
"6"  
"6"  
"6"  
"6"  
"6"  
"10" 
"10" 
"10" 
"10" 
"10" 
"7"  
"7"  
"7"  
"7"  
"7"  
"6"  
"6"  
"6"  
"6"  
"6"  
"6"  
"6"  
"6"  
"6"  
"6"  
"6"  
"6"  
"6"  
"6"  
"6"  
"4"  
"4"  
"4"  
"4"  
"4"  
"6"  
"6"  
"6"  
"6"  
"6"  
"4"  
"4"  
"4"  
"4"  
"4"  
"10" 
"10" 
"10" 
"10" 
"10" 
"5"  
"5"  
"5"  
"5"  
"5"  
"2"  
"2"  
"2"  
"2"  
"2"  
"11" 
"11" 
"11" 
"11" 
"11" 
"4"  
"4"  
"4"  
"4"  
"4"  
"5"  
"5"  
"5"  
"5"  
"5"  
"9"  
"9"  
"9"  
"9"  
"9"  
"3"  
"3"  
"3"  
"3"  
"3"  
"7"  
"7"  
"7"  
"7"  
"7"  
"8"  
"8"  
"8"  
"8"  
"8"  
"8"  
"8"  
"8"  
"8"  
"8"  
"11" 
"11" 
"11" 
"11" 
"11" 
)



# new instances, now on 100 vertices

idx=0
for entry in "${instances_gnp_100_only50[@]}";
do
    timestamp=$(date)
    $(echo "[$timestamp]  $entry" >> "$entry""_xp8a.out")
    echo "[$timestamp]  $entry"

    output=$(./cks_xp8_full  $entry  ${instances_gnp_100_num_colours[$idx]} >> "$entry""_xp8a.out" 2>&1)
    echo "$output"

    timestamp=$(date)
    $(echo -e "[$timestamp]  done\n" >> "$entry""_xp8a.out")

    ((++idx))
done


idx=0
for entry in "${instances_gnp_100_only50[@]}";
do
    timestamp=$(date)
    $(echo "[$timestamp]  $entry" >> "$entry""_xp8b.out")
    echo "[$timestamp]  $entry"

    output=$(./cks_xp8_msi  $entry  ${instances_gnp_100_num_colours[$idx]} >> "$entry""_xp8b.out" 2>&1)
    echo "$output"

    timestamp=$(date)
    $(echo -e "[$timestamp]  done\n" >> "$entry""_xp8b.out")

    ((++idx))
done


idx=0
for entry in "${instances_gnp_100_only50[@]}";
do
    timestamp=$(date)
    $(echo "[$timestamp]  $entry" >> "$entry""_xp9a.out")
    echo "[$timestamp]  $entry"

    output=$(./cks_xp9_full  $entry  ${instances_gnp_100_num_colours[$idx]} >> "$entry""_xp9a.out" 2>&1)
    echo "$output"

    timestamp=$(date)
    $(echo -e "[$timestamp]  done\n" >> "$entry""_xp9a.out")

    ((++idx))
done


idx=0
for entry in "${instances_gnp_100_only50[@]}";
do
    timestamp=$(date)
    $(echo "[$timestamp]  $entry" >> "$entry""_xp9b.out")
    echo "[$timestamp]  $entry"

    output=$(./cks_xp9_msi  $entry  ${instances_gnp_100_num_colours[$idx]} >> "$entry""_xp9b.out" 2>&1)
    echo "$output"

    timestamp=$(date)
    $(echo -e "[$timestamp]  done\n" >> "$entry""_xp9b.out")

    ((++idx))
done
