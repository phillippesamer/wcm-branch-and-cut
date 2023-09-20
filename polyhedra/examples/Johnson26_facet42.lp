I wrote on the .lp file
(2,6)-SEP_#3: x1 + x5 + x6 + x7 +x11+x12+x13 -x1 -x2 -x3 -x4 -x2 -x8 -x9 -x3 -x6 -x11 -x7 -x10 -x12 -x14 <= 1
(2,8)-SEP_#2: x1 + x5 + x6 + x7 +x4+x9+x13+x14 -x1 -x2 -x3 -x4 -x2 -x8 -x9 -x3 -x6 -x11 -x7 -x10 -x12 -x14 <= 1
(4,6)-SEP_#2: x5 + x8 + x10 +x11+x12+x13 -x1 -x2 -x3 -x4 -x2 -x8 -x9 -x3 -x6 -x11 -x7 -x10 -x12 -x14 <= 1
(4,8)-SEP_#2: x5 + x8 + x10 +x4+x9+x13+x14 -x1 -x2 -x3 -x4 -x2 -x8 -x9 -x3 -x6 -x11 -x7 -x10 -x12 -x14 <= 1

We see that {1,3,5,7} is the separator in all four cases! That is, dropping the positive terms we get
(2,6)-SEP_#3:  -x1 -x2 -x3 -x4      -x2 -x8 -x9       -x3 -x6 -x11      -x7 -x10 -x12 -x14 <= 1
(2,8)-SEP_#2:  -x1 -x2 -x3 -x4      -x2 -x8 -x9       -x3 -x6 -x11      -x7 -x10 -x12 -x14 <= 1
(4,6)-SEP_#2:  -x1 -x2 -x3 -x4      -x2 -x8 -x9       -x3 -x6 -x11      -x7 -x10 -x12 -x14 <= 1
(4,8)-SEP_#2:  -x1 -x2 -x3 -x4      -x2 -x8 -x9       -x3 -x6 -x11      -x7 -x10 -x12 -x14 <= 1

Rearranging, the actual MSI indicated in the .lp file are
(2,6)-SEP_#3: x5+x13 -2x2 -2x3     -x4         -x8 -x9 -x10           -x14 <= 1
(2,8)-SEP_#2: x5+x13 -2x2 -2x3                 -x8     -x10 -x11 -x12      <= 1
(4,6)-SEP_#2: x5+x13 -2x2 -2x3 -x1 -x4 -x6 -x7     -x9                -x14 <= 1
(4,8)-SEP_#2: x5+x13 -2x2 -2x3 -x1     -x6 -x7              -x11 -x12      <= 1

In turn, polymake finds the facet-defining inequality
42: x2 + x3 - x5 - x13 >= -1
that is,
42: x5 + x13 -x2 -x3 <= 1
which dominates all four MSI.
