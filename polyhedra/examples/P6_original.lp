Maximize
   x1 +  x2 +  x3 +  x4 +  x5 

Subject To
DEG1: x1 <= 1
DEG2: x1 + x2 <= 1
DEG3: x2 + x3 <= 1
DEG4: x3 + x4 <= 1
DEG5: x4 + x5 <= 1
DEG6: x5 <= 1
(1,3)-SEP_#1: x1 +x2+x3 -x1 -x2 <= 1
(1,4)-SEP_#1: x1 +x3+x4 -x1 -x2 <= 1
(1,4)-SEP_#2: x1 +x3+x4 -x2 -x3 <= 1
(1,5)-SEP_#1: x1 +x4+x5 -x1 -x2 <= 1
(1,5)-SEP_#2: x1 +x4+x5 -x2 -x3 <= 1
(1,5)-SEP_#3: x1 +x4+x5 -x3 -x4 <= 1
(1,6)-SEP_#1: x1 +x5 -x1 -x2 <= 1
(1,6)-SEP_#2: x1 +x5 -x2 -x3 <= 1
(1,6)-SEP_#3: x1 +x5 -x3 -x4 <= 1
(1,6)-SEP_#4: x1 +x5 -x4 -x5 <= 1
(2,4)-SEP_#1: x1 + x2 +x3+x4 -x2 -x3 <= 1
(2,5)-SEP_#1: x1 + x2 +x4+x5 -x2 -x3 <= 1
(2,5)-SEP_#2: x1 + x2 +x4+x5 -x3 -x4 <= 1
(2,6)-SEP_#1: x1 + x2 +x5 -x2 -x3 <= 1
(2,6)-SEP_#2: x1 + x2 +x5 -x3 -x4 <= 1
(2,6)-SEP_#3: x1 + x2 +x5 -x4 -x5 <= 1
(3,5)-SEP_#1: x2 + x3 +x4+x5 -x3 -x4 <= 1
(3,6)-SEP_#1: x2 + x3 +x5 -x3 -x4 <= 1
(3,6)-SEP_#2: x2 + x3 +x5 -x4 -x5 <= 1
(4,6)-SEP_#1: x3 + x4 +x5 -x4 -x5 <= 1

Bounds

Binaries
 x1  x2  x3  x4  x5 

End
