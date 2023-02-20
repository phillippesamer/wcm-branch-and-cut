#!/bin/bash

cr_selection=( 
"./input/itor2020-instances/10/010_010_001.gcc   ./input/itor2020-bounds/10/010_010_001.gcc  "   
"./input/itor2020-instances/10/010_020_001.gcc   ./input/itor2020-bounds/10/010_020_001.gcc  "   
"./input/itor2020-instances/10/010_030_001.gcc   ./input/itor2020-bounds/10/010_030_001.gcc  "   
"./input/itor2020-instances/10/010_040_001.gcc   ./input/itor2020-bounds/10/010_040_001.gcc  "   
"./input/itor2020-instances/10/010_050_001.gcc   ./input/itor2020-bounds/10/010_050_001.gcc  "   
"./input/itor2020-instances/20/020_010_001.gcc   ./input/itor2020-bounds/20/020_010_001.gcc  "   
"./input/itor2020-instances/20/020_020_001.gcc   ./input/itor2020-bounds/20/020_020_001.gcc  "   
"./input/itor2020-instances/20/020_030_001.gcc   ./input/itor2020-bounds/20/020_030_001.gcc  "   
"./input/itor2020-instances/20/020_040_001.gcc   ./input/itor2020-bounds/20/020_040_001.gcc  "   
"./input/itor2020-instances/20/020_050_001.gcc   ./input/itor2020-bounds/20/020_050_001.gcc  "   
"./input/itor2020-instances/30/030_010_001.gcc   ./input/itor2020-bounds/30/030_010_001.gcc  "   
"./input/itor2020-instances/30/030_020_001.gcc   ./input/itor2020-bounds/30/030_020_001.gcc  "   
"./input/itor2020-instances/30/030_030_001.gcc   ./input/itor2020-bounds/30/030_030_001.gcc  "   
"./input/itor2020-instances/30/030_040_001.gcc   ./input/itor2020-bounds/30/030_040_001.gcc  "   
"./input/itor2020-instances/30/030_050_001.gcc   ./input/itor2020-bounds/30/030_050_001.gcc  "   
"./input/itor2020-instances/40/040_010_001.gcc   ./input/itor2020-bounds/40/040_010_001.gcc  "   
"./input/itor2020-instances/40/040_020_001.gcc   ./input/itor2020-bounds/40/040_020_001.gcc  "   
"./input/itor2020-instances/40/040_030_001.gcc   ./input/itor2020-bounds/40/040_030_001.gcc  "   
"./input/itor2020-instances/40/040_040_001.gcc   ./input/itor2020-bounds/40/040_040_001.gcc  "   
"./input/itor2020-instances/40/040_050_001.gcc   ./input/itor2020-bounds/40/040_050_001.gcc  "   
"./input/itor2020-instances/50/050_010_001.gcc   ./input/itor2020-bounds/50/050_010_001.gcc  "   
"./input/itor2020-instances/50/050_020_001.gcc   ./input/itor2020-bounds/50/050_020_001.gcc  "   
"./input/itor2020-instances/50/050_030_001.gcc   ./input/itor2020-bounds/50/050_030_001.gcc  "   
"./input/itor2020-instances/50/050_040_001.gcc   ./input/itor2020-bounds/50/050_040_001.gcc  "   
"./input/itor2020-instances/50/050_050_001.gcc   ./input/itor2020-bounds/50/050_050_001.gcc  "   
"./input/itor2020-instances/60/060_010_001.gcc   ./input/itor2020-bounds/60/060_010_001.gcc  "   
"./input/itor2020-instances/60/060_020_001.gcc   ./input/itor2020-bounds/60/060_020_001.gcc  "   
"./input/itor2020-instances/60/060_030_001.gcc   ./input/itor2020-bounds/60/060_030_001.gcc  "   
"./input/itor2020-instances/60/060_040_001.gcc   ./input/itor2020-bounds/60/060_040_001.gcc  "   
"./input/itor2020-instances/60/060_050_001.gcc   ./input/itor2020-bounds/60/060_050_001.gcc  "   
"./input/itor2020-instances/70/070_010_001.gcc   ./input/itor2020-bounds/70/070_010_001.gcc  "   
"./input/itor2020-instances/70/070_020_001.gcc   ./input/itor2020-bounds/70/070_020_001.gcc  "   
"./input/itor2020-instances/70/070_030_001.gcc   ./input/itor2020-bounds/70/070_030_001.gcc  "   
"./input/itor2020-instances/70/070_040_001.gcc   ./input/itor2020-bounds/70/070_040_001.gcc  "   
"./input/itor2020-instances/70/070_050_001.gcc   ./input/itor2020-bounds/70/070_050_001.gcc  "   
"./input/itor2020-instances/80/080_010_001.gcc   ./input/itor2020-bounds/80/080_010_001.gcc  "   
"./input/itor2020-instances/80/080_020_001.gcc   ./input/itor2020-bounds/80/080_020_001.gcc  "   
"./input/itor2020-instances/80/080_030_001.gcc   ./input/itor2020-bounds/80/080_030_001.gcc  "   
"./input/itor2020-instances/80/080_040_001.gcc   ./input/itor2020-bounds/80/080_040_001.gcc  "   
"./input/itor2020-instances/80/080_050_001.gcc   ./input/itor2020-bounds/80/080_050_001.gcc  "   
"./input/itor2020-instances/90/090_010_001.gcc   ./input/itor2020-bounds/90/090_010_001.gcc  "   
"./input/itor2020-instances/90/090_020_001.gcc   ./input/itor2020-bounds/90/090_020_001.gcc  "   
"./input/itor2020-instances/90/090_030_001.gcc   ./input/itor2020-bounds/90/090_030_001.gcc  "   
"./input/itor2020-instances/90/090_040_001.gcc   ./input/itor2020-bounds/90/090_040_001.gcc  "   
"./input/itor2020-instances/90/090_050_001.gcc   ./input/itor2020-bounds/90/090_050_001.gcc  "   
"./input/itor2020-instances/100/100_010_001.gcc  ./input/itor2020-bounds/100/100_010_001.gcc "   
"./input/itor2020-instances/100/100_020_001.gcc  ./input/itor2020-bounds/100/100_020_001.gcc "   
"./input/itor2020-instances/100/100_030_001.gcc  ./input/itor2020-bounds/100/100_030_001.gcc "   
"./input/itor2020-instances/100/100_040_001.gcc  ./input/itor2020-bounds/100/100_040_001.gcc "   
"./input/itor2020-instances/100/100_050_001.gcc  ./input/itor2020-bounds/100/100_050_001.gcc "   
)

for entry in "${cr_selection[@]}";
do
    timestamp=$(date)
    $(echo "[$timestamp]  $entry" >> "$entry""_xp1.out")
    echo "[$timestamp]  $entry"

    output=$(./cks $entry >> "$entry""_xp1.out" 2>&1)
    echo "$output"

    timestamp=$(date)
    $(echo -e "[$timestamp]  done\n" >> "$entry""_xp1.out")
done
