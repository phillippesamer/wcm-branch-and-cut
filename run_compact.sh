#!/bin/bash

instances_selection_of_21=(
"./input/dimacs-GMWCS-GAM/1255693fce68.stp  -e"
"./input/dimacs-GMWCS-GAM/1fe159187d94.stp  -e"
"./input/dimacs-GMWCS-GAM/723a287b9872.stp  -e"
"./input/dimacs-GMWCS-GAM/723a2be354e1.stp  -e"
"./input/dimacs-GMWCS-GAM/72ae7cfa16cd.stp  -e"
"./input/dimacs-GMWCS-GAM/79536db666e0.stp  -e"
"./input/dimacs-MWCS-GAM/3a0d4427fe32.stp"
"./input/dimacs-MWCS-GAM/3a0dff0eb70.stp"
"./input/dimacs-MWCS-GAM/3a0d724ffec9.stp"
"./input/dimacs-MWCS-GAM/48e76a6886bc.stp"
"./input/dimacs-MWCS-GAM/3a0d2255a681.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-1000-a-0.6-d-0.5-e-0.5.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-1500-a-0.6-d-0.5-e-0.5.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-1000-a-1-d-0.25-e-0.25.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-1000-a-1-d-0.75-e-0.75.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-1500-a-1-d-0.25-e-0.75.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-1500-a-1-d-0.75-e-0.25.stp"
"./input/dimacs-MWCS-ACTMOD/drosophila001.stp"
"./input/dimacs-MWCS-ACTMOD/HCMV.stp"
"./input/dimacs-MWCS-ACTMOD/metabol_expr_mice_1.stp"
"./input/dimacs-MWCS-ACTMOD/metabol_expr_mice_3.stp"
)

instances_dimacs_mwcs_actmod=(
"./input/dimacs-MWCS-ACTMOD/drosophila001.stp"
"./input/dimacs-MWCS-ACTMOD/drosophila005.stp"
"./input/dimacs-MWCS-ACTMOD/drosophila0075.stp"
"./input/dimacs-MWCS-ACTMOD/HCMV.stp"
"./input/dimacs-MWCS-ACTMOD/lymphoma.stp"
"./input/dimacs-MWCS-ACTMOD/metabol_expr_mice_1.stp"
"./input/dimacs-MWCS-ACTMOD/metabol_expr_mice_2.stp"
"./input/dimacs-MWCS-ACTMOD/metabol_expr_mice_3.stp"
)

instances_dimacs_mwcs_gam=(
"./input/dimacs-MWCS-GAM/25e814a792c4.stp"
"./input/dimacs-MWCS-GAM/25e81700dead.stp"
"./input/dimacs-MWCS-GAM/25e83661bc4.stp"
"./input/dimacs-MWCS-GAM/25e83d7dbeea.stp"
"./input/dimacs-MWCS-GAM/25e857e14393.stp"
"./input/dimacs-MWCS-GAM/3a0d1335fe78.stp"
"./input/dimacs-MWCS-GAM/3a0d151a8ee0.stp"
"./input/dimacs-MWCS-GAM/3a0d17a83362.stp"
"./input/dimacs-MWCS-GAM/3a0d1a1e31cf.stp"
"./input/dimacs-MWCS-GAM/3a0d2255a681.stp"
"./input/dimacs-MWCS-GAM/3a0d226a0a5c.stp"
"./input/dimacs-MWCS-GAM/3a0d25c9a738.stp"
"./input/dimacs-MWCS-GAM/3a0d25f9bda3.stp"
"./input/dimacs-MWCS-GAM/3a0d2875c8cf.stp"
"./input/dimacs-MWCS-GAM/3a0d325af5cc.stp"
"./input/dimacs-MWCS-GAM/3a0d32b18854.stp"
"./input/dimacs-MWCS-GAM/3a0d33d2aa32.stp"
"./input/dimacs-MWCS-GAM/3a0d390c537e.stp"
"./input/dimacs-MWCS-GAM/3a0d435ee480.stp"
"./input/dimacs-MWCS-GAM/3a0d4427fe32.stp"
"./input/dimacs-MWCS-GAM/3a0d4ccc9b37.stp"
"./input/dimacs-MWCS-GAM/3a0d4dac5319.stp"
"./input/dimacs-MWCS-GAM/3a0d52ee8185.stp"
"./input/dimacs-MWCS-GAM/3a0d55ddd0a5.stp"
"./input/dimacs-MWCS-GAM/3a0d568fbd87.stp"
"./input/dimacs-MWCS-GAM/3a0d5dc4a759.stp"
"./input/dimacs-MWCS-GAM/3a0d5e4aac27.stp"
"./input/dimacs-MWCS-GAM/3a0d5e4aac27x.stp"
"./input/dimacs-MWCS-GAM/3a0d610beb4c.stp"
"./input/dimacs-MWCS-GAM/3a0d6505353b.stp"
"./input/dimacs-MWCS-GAM/3a0d6a21bbd5.stp"
"./input/dimacs-MWCS-GAM/3a0d6e97602a.stp"
"./input/dimacs-MWCS-GAM/3a0d724ffec9.stp"
"./input/dimacs-MWCS-GAM/3a0d73143aeb.stp"
"./input/dimacs-MWCS-GAM/3a0dff0eb70.stp"
"./input/dimacs-MWCS-GAM/48e7452da6ba.stp"
"./input/dimacs-MWCS-GAM/48e7526364af.stp"
"./input/dimacs-MWCS-GAM/48e76a6886bc.stp"
"./input/dimacs-MWCS-GAM/795313fd138b.stp"
)

instances_dimacs_mwcs_jmpalmk=(
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-1000-a-0.6-d-0.25-e-0.25.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-1000-a-0.6-d-0.25-e-0.5.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-1000-a-0.6-d-0.25-e-0.75.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-1000-a-0.6-d-0.5-e-0.25.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-1000-a-0.6-d-0.5-e-0.5.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-1000-a-0.6-d-0.5-e-0.75.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-1000-a-0.6-d-0.75-e-0.25.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-1000-a-0.6-d-0.75-e-0.5.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-1000-a-0.6-d-0.75-e-0.75.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-1000-a-1-d-0.25-e-0.25.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-1000-a-1-d-0.25-e-0.5.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-1000-a-1-d-0.25-e-0.75.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-1000-a-1-d-0.5-e-0.25.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-1000-a-1-d-0.5-e-0.5.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-1000-a-1-d-0.5-e-0.75.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-1000-a-1-d-0.75-e-0.25.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-1000-a-1-d-0.75-e-0.5.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-1000-a-1-d-0.75-e-0.75.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-1500-a-0.6-d-0.25-e-0.25.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-1500-a-0.6-d-0.25-e-0.5.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-1500-a-0.6-d-0.25-e-0.75.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-1500-a-0.6-d-0.5-e-0.25.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-1500-a-0.6-d-0.5-e-0.5.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-1500-a-0.6-d-0.5-e-0.75.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-1500-a-0.6-d-0.75-e-0.25.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-1500-a-0.6-d-0.75-e-0.5.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-1500-a-0.6-d-0.75-e-0.75.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-1500-a-1-d-0.25-e-0.25.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-1500-a-1-d-0.25-e-0.5.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-1500-a-1-d-0.25-e-0.75.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-1500-a-1-d-0.5-e-0.25.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-1500-a-1-d-0.5-e-0.5.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-1500-a-1-d-0.5-e-0.75.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-1500-a-1-d-0.75-e-0.25.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-1500-a-1-d-0.75-e-0.5.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-1500-a-1-d-0.75-e-0.75.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-500-a-0.62-d-0.25-e-0.25.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-500-a-0.62-d-0.25-e-0.5.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-500-a-0.62-d-0.25-e-0.75.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-500-a-0.62-d-0.5-e-0.25.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-500-a-0.62-d-0.5-e-0.5.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-500-a-0.62-d-0.5-e-0.75.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-500-a-0.62-d-0.75-e-0.25.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-500-a-0.62-d-0.75-e-0.5.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-500-a-0.62-d-0.75-e-0.75.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-500-a-1-d-0.25-e-0.25.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-500-a-1-d-0.25-e-0.5.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-500-a-1-d-0.25-e-0.75.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-500-a-1-d-0.5-e-0.25.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-500-a-1-d-0.5-e-0.5.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-500-a-1-d-0.5-e-0.75.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-500-a-1-d-0.75-e-0.25.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-500-a-1-d-0.75-e-0.5.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-500-a-1-d-0.75-e-0.75.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-750-a-0.647-d-0.25-e-0.25.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-750-a-0.647-d-0.25-e-0.5.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-750-a-0.647-d-0.25-e-0.75.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-750-a-0.647-d-0.5-e-0.25.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-750-a-0.647-d-0.5-e-0.5.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-750-a-0.647-d-0.5-e-0.75.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-750-a-0.647-d-0.75-e-0.25.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-750-a-0.647-d-0.75-e-0.5.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-750-a-0.647-d-0.75-e-0.75.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-750-a-1-d-0.25-e-0.25.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-750-a-1-d-0.25-e-0.5.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-750-a-1-d-0.25-e-0.75.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-750-a-1-d-0.5-e-0.25.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-750-a-1-d-0.5-e-0.5.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-750-a-1-d-0.5-e-0.75.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-750-a-1-d-0.75-e-0.25.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-750-a-1-d-0.75-e-0.5.stp"
"./input/dimacs-MWCS-JMPALMK/MWCS-I-D-n-750-a-1-d-0.75-e-0.75.stp"
)

instances_dimacs_gmwcs_gam=(
"./input/dimacs-GMWCS-GAM/1255693fce68.stp  -e"
"./input/dimacs-GMWCS-GAM/1fe159187d94.stp  -e"
"./input/dimacs-GMWCS-GAM/25e811bd74c7.stp  -e"
"./input/dimacs-GMWCS-GAM/25e811f2a7c9.stp  -e"
"./input/dimacs-GMWCS-GAM/25e8203fb7de.stp  -e"
"./input/dimacs-GMWCS-GAM/25e8226bf1cf.stp  -e"
"./input/dimacs-GMWCS-GAM/25e828ae41d9.stp  -e"
"./input/dimacs-GMWCS-GAM/25e8295af368.stp  -e"
"./input/dimacs-GMWCS-GAM/25e8344d42b8.stp  -e"
"./input/dimacs-GMWCS-GAM/25e83737e486.stp  -e"
"./input/dimacs-GMWCS-GAM/25e83aa98746.stp  -e"
"./input/dimacs-GMWCS-GAM/25e84772ec89.stp  -e"
"./input/dimacs-GMWCS-GAM/25e84f25f309.stp  -e"
"./input/dimacs-GMWCS-GAM/25e853c13426.stp  -e"
"./input/dimacs-GMWCS-GAM/25e8543881d1.stp  -e"
"./input/dimacs-GMWCS-GAM/25e85970c8bd.stp  -e"
"./input/dimacs-GMWCS-GAM/25e861a978f3.stp  -e"
"./input/dimacs-GMWCS-GAM/25e87515b1ff.stp  -e"
"./input/dimacs-GMWCS-GAM/25e875a3e06.stp   -e"
"./input/dimacs-GMWCS-GAM/25e87a601297.stp  -e"
"./input/dimacs-GMWCS-GAM/25e87e6a7023.stp  -e"
"./input/dimacs-GMWCS-GAM/25e87fbbac1f.stp  -e"
"./input/dimacs-GMWCS-GAM/25e8b025a25.stp   -e"
"./input/dimacs-GMWCS-GAM/723a125f899e.stp  -e"
"./input/dimacs-GMWCS-GAM/723a13c9c13.stp   -e"
"./input/dimacs-GMWCS-GAM/723a179e1563.stp  -e"
"./input/dimacs-GMWCS-GAM/723a1e88a66c.stp  -e"
"./input/dimacs-GMWCS-GAM/723a1ed5056a.stp  -e"
"./input/dimacs-GMWCS-GAM/723a21a4ba8f.stp  -e"
"./input/dimacs-GMWCS-GAM/723a24072832.stp  -e"
"./input/dimacs-GMWCS-GAM/723a2411c60c.stp  -e"
"./input/dimacs-GMWCS-GAM/723a287b9872.stp  -e"
"./input/dimacs-GMWCS-GAM/723a2be354e1.stp  -e"
"./input/dimacs-GMWCS-GAM/723a2f0af854.stp  -e"
"./input/dimacs-GMWCS-GAM/723a2fa882ad.stp  -e"
"./input/dimacs-GMWCS-GAM/723a3547be29.stp  -e"
"./input/dimacs-GMWCS-GAM/723a36cbd1.stp    -e"
"./input/dimacs-GMWCS-GAM/723a38eeece6.stp  -e"
"./input/dimacs-GMWCS-GAM/723a3ac69817.stp  -e"
"./input/dimacs-GMWCS-GAM/723a3b793235.stp  -e"
"./input/dimacs-GMWCS-GAM/723a4368043d.stp  -e"
"./input/dimacs-GMWCS-GAM/723a4a511045.stp  -e"
"./input/dimacs-GMWCS-GAM/723a4b4ca686.stp  -e"
"./input/dimacs-GMWCS-GAM/723a4bd22e72.stp  -e"
"./input/dimacs-GMWCS-GAM/723a4c5bc01.stp   -e"
"./input/dimacs-GMWCS-GAM/723a4eebd776.stp  -e"
"./input/dimacs-GMWCS-GAM/723a5541f17b.stp  -e"
"./input/dimacs-GMWCS-GAM/723a575cdc90.stp  -e"
"./input/dimacs-GMWCS-GAM/723a58b5a879.stp  -e"
"./input/dimacs-GMWCS-GAM/723a5c31d346.stp  -e"
"./input/dimacs-GMWCS-GAM/723a5d6eec94.stp  -e"
"./input/dimacs-GMWCS-GAM/723a5e7c968c.stp  -e"
"./input/dimacs-GMWCS-GAM/723a5ee26a8e.stp  -e"
"./input/dimacs-GMWCS-GAM/723a65f1f1e0.stp  -e"
"./input/dimacs-GMWCS-GAM/723a69873c1a.stp  -e"
"./input/dimacs-GMWCS-GAM/723a79d1673b.stp  -e"
"./input/dimacs-GMWCS-GAM/723a7cc8a8c4.stp  -e"
"./input/dimacs-GMWCS-GAM/72ae15bb11d6.stp  -e"
"./input/dimacs-GMWCS-GAM/72ae2402493.stp   -e"
"./input/dimacs-GMWCS-GAM/72ae4a8792fb.stp  -e"
"./input/dimacs-GMWCS-GAM/72ae6c7de2d7.stp  -e"
"./input/dimacs-GMWCS-GAM/72ae7cfa16cd.stp  -e"
"./input/dimacs-GMWCS-GAM/79536db666e0.stp  -e"
)

# xp14:

idx=1
for entry in "${instances_dimacs_gmwcs_gam[@]}";
do
    timestamp=$(date)
    $(echo "[$timestamp]  $entry" >> "$entry""_xp14.out")
    echo "[$timestamp] $idx/${#instances_dimacs_gmwcs_gam[@]}:  $entry"

    output=$(./wcm  $entry  >> "$entry""_xp14.out" 2>&1)
    echo "$output"

    timestamp=$(date)
    $(echo -e "[$timestamp]  done\n" >> "$entry""_xp14.out")

    ((++idx))
done


idx=1
for entry in "${instances_dimacs_mwcs_jmpalmk[@]}";
do
    timestamp=$(date)
    $(echo "[$timestamp]  $entry" >> "$entry""_xp14.out")
    echo "[$timestamp] $idx/${#instances_dimacs_mwcs_jmpalmk[@]}:  $entry"

    output=$(./wcm  $entry  >> "$entry""_xp14.out" 2>&1)
    echo "$output"

    timestamp=$(date)
    $(echo -e "[$timestamp]  done\n" >> "$entry""_xp14.out")

    ((++idx))
done


idx=1
for entry in "${instances_dimacs_mwcs_actmod[@]}";
do
    timestamp=$(date)
    $(echo "[$timestamp]  $entry" >> "$entry""_xp14.out")
    echo "[$timestamp] $idx/${#instances_dimacs_mwcs_actmod[@]}:  $entry"

    output=$(./wcm  $entry  >> "$entry""_xp14.out" 2>&1)
    echo "$output"

    timestamp=$(date)
    $(echo -e "[$timestamp]  done\n" >> "$entry""_xp14.out")

    ((++idx))
done


idx=1
for entry in "${instances_dimacs_mwcs_gam[@]}";
do
    timestamp=$(date)
    $(echo "[$timestamp]  $entry" >> "$entry""_xp14.out")
    echo "[$timestamp] $idx/${#instances_dimacs_mwcs_gam[@]}:  $entry"

    output=$(./wcm  $entry  >> "$entry""_xp14.out" 2>&1)
    echo "$output"

    timestamp=$(date)
    $(echo -e "[$timestamp]  done\n" >> "$entry""_xp14.out")

    ((++idx))
done
