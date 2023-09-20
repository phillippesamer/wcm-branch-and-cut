Maximize
   x1 +  x2 +  x3 +  x4 +  x5 +  x6 +  x7 +  x8 

Subject To
DEG1: x1 + x2 + x3 + x4 <= 1
DEG2: x5 + x6 + x7 + x8 <= 1
DEG3: x1 + x5 <= 1
DEG4: x2 + x6 <= 1
DEG5: x3 + x7 <= 1
DEG6: x4 + x8 <= 1
(1,2)-SEP_#1: x1 + x2 + x3 + x4 +x5+x6+x7+x8 -x1 -x5 -x2 -x6 -x3 -x7 -x4 -x8 <= 1
(3,4)-SEP_#1: x1 + x5 +x2+x6 -x1 -x2 -x3 -x4 -x5 -x6 -x7 -x8 <= 1
(3,5)-SEP_#1: x1 + x5 +x3+x7 -x1 -x2 -x3 -x4 -x5 -x6 -x7 -x8 <= 1
(3,6)-SEP_#1: x1 + x5 +x4+x8 -x1 -x2 -x3 -x4 -x5 -x6 -x7 -x8 <= 1
(4,5)-SEP_#1: x2 + x6 +x3+x7 -x1 -x2 -x3 -x4 -x5 -x6 -x7 -x8 <= 1
(4,6)-SEP_#1: x2 + x6 +x4+x8 -x1 -x2 -x3 -x4 -x5 -x6 -x7 -x8 <= 1
(5,6)-SEP_#1: x3 + x7 +x4+x8 -x1 -x2 -x3 -x4 -x5 -x6 -x7 -x8 <= 1

Bounds

Binaries
 x1  x2  x3  x4  x5  x6  x7  x8 

End
