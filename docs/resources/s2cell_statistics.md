---
title: S2 Cell Statistics
---

These are the approximate ranges of areas of
[S2 cells](/about/overview.md#s2cell-hierarchy)
at each level. The average size is that returned by `S2Cell.AverageArea()`,
guaranteed to be within a factor of 1.5 of the high and low end of the range.
Additionally, two random locations were chosen to provide example cell edge
lengths.

| **level** | **min area** | **max area** | **average area** | **units** |         | **Random cell 1 (UK) min edge length** | **Random cell 1 (UK) max edge length** | | **Random cell 2 (US) min edge length** | **Random cell 2 (US) max edge length** | **Number of cells** |
| --- | ---: | ---: | ---: | ---: | --- | ---: | ---: | --- | ---: | ---: | ---: |
| 00 |  85011012.19 |  85011012.19 |  85011012.19 | km<sup>2</sup> |       |    7842 km |   7842 km | |    7842 km |   7842 km |  6 |
| 01 |  21252753.05 |  21252753.05 |  21252753.05 | km<sup>2</sup> |       |    3921 km |   5004 km | |    3921 km |   5004 km |  24 |
| 02 |  4919708.23 |  6026521.16 |  5313188.26 | km<sup>2</sup> |  |    1825 km |   2489 km | |    1825 km |   2489 km |  96 |
| 03 |  1055377.48 |  1646455.50 |  1328297.07 | km<sup>2</sup> |  |     840 km |   1167 km | |    1130 km |   1310 km |  384 |
| 04 |   231564.06 |   413918.15 |   332074.27 | km<sup>2</sup> |  |     432 km |    609 km | |     579 km |    636 km |  1536 |
| 05 |    53798.67 |   104297.91 |    83018.57 | km<sup>2</sup> |  |     210 km |    298 km | |     287 km |    315 km |  6K |
| 06 |    12948.81 |    26113.30 |    20754.64 | km<sup>2</sup> |  |     108 km |    151 km | |     143 km |    156 km |  24K |
| 07 |     3175.44 |     6529.09 |     5188.66 | km<sup>2</sup> |  |      54 km |     76 km | |      72 km |     78 km |  98K |
| 08 |      786.20 |     1632.45 |     1297.17 | km<sup>2</sup> |  |      27 km |     38 km | |      36 km |     39 km |  393K |
| 09 |      195.59 |      408.12 |      324.29 | km<sup>2</sup> |  |      14 km |     19 km | |      18 km |     20 km |  1573K |
| 10 |       48.78 |      102.03 |       81.07 | km<sup>2</sup> |  |       7 km |      9 km | |       9 km |     10 km |  6M |
| 11 |       12.18 |       25.51 |       20.27 | km<sup>2</sup> |  |       3 km |      5 km | |       4 km |      5 km |  25M |
| 12 |        3.04 |        6.38 |        5.07 | km<sup>2</sup> |  |    1699  m |      2 km | |       2 km |      2 km |  100M |
| 13 |        0.76 |        1.59 |        1.27 | km<sup>2</sup> |  |     850  m |    1185 m | |    1123  m |    1225 m |  402M |
| 14 |        0.19 |        0.40 |        0.32 | km<sup>2</sup> |  |     425  m |    593  m | |     562  m |    613  m |  1610M |
| 15 |    47520.30 |    99638.93 |    79172.67 | m<sup>2</sup> |   |     212  m |    296  m | |     281  m |    306  m |  6B |
| 16 |    11880.08 |    24909.73 |    19793.17 | m<sup>2</sup> |   |     106  m |    148  m | |     140  m |    153  m |  25B |
| 17 |     2970.02 |     6227.43 |     4948.29 | m<sup>2</sup> |   |      53  m |     74  m | |      70  m |     77  m |  103B |
| 18 |      742.50 |     1556.86 |     1237.07 | m<sup>2</sup> |   |      27  m |     37  m | |      35  m |     38  m |  412B |
| 19 |      185.63 |      389.21 |      309.27 | m<sup>2</sup> |   |      13  m |     19  m | |      18  m |     19  m |  1649B |
| 20 |       46.41 |       97.30 |       77.32 | m<sup>2</sup> |   |       7  m |      9  m | |       9  m |     10  m |  7T |
| 21 |       11.60 |       24.33 |       19.33 | m<sup>2</sup> |   |       3  m |      5  m | |       4  m |      5  m |  26T |
| 22 |        2.90 |        6.08 |        4.83 | m<sup>2</sup> |   |     166 cm |      2  m | |       2  m |      2  m |  105T |
| 23 |        0.73 |        1.52 |        1.21 | m<sup>2</sup> |   |      83 cm |    116 cm | |     110 cm |    120 cm |  422T |
| 24 |        0.18 |        0.38 |        0.30 | m<sup>2</sup> |   |      41 cm |     58 cm | |      55 cm |     60 cm |  1689T |
| 25 |      453.19 |      950.23 |      755.05 | cm<sup>2</sup> |  |      21 cm |     29 cm | |      27 cm |     30 cm |  7e15 |
| 26 |      113.30 |      237.56 |      188.76 | cm<sup>2</sup> |  |      10 cm |     14 cm | |      14 cm |     15 cm |  27e15 |
| 27 |       28.32 |       59.39 |       47.19 | cm<sup>2</sup> |  |       5 cm |      7 cm | |       7 cm |      7 cm |  108e15 |
| 28 |        7.08 |       14.85 |       11.80 | cm<sup>2</sup> |  |       2 cm |      4 cm | |       3 cm |      4 cm |  432e15 |
| 29 |        1.77 |        3.71 |        2.95 | cm<sup>2</sup> |  |      12 mm |     18 mm | |      17 mm |     18 mm |  1729e15 |
| 30 |        0.44 |        0.93 |        0.74 | cm<sup>2</sup> |  |       6 mm |      9 mm | |       8 mm |      9 mm |  7e18 |
