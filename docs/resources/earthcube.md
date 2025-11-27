---
title: Earth Cube
---

The S2Cell hierarchy projects the Earth onto six top-level "face cells",
which may be unwrapped and flattened into a stairstep arrangement. Below
is a map of the six faces:

![](/devguide/img/s2cell_global.jpg)

A single space-filling curve covers the entire globe, snaking from one
face to the next. Note the traversal order in the odd-numbered faces is
the mirror of the order in even numbered faces.

~~~~
                  ^---->  <----^
                  | 4  |    5  |
                  |    |       |
                  x    V  ----->

          ^---->  <----^
          | 2  |    3  |
          |    |       |
          x    V  ----->

  ^---->  <----^
  | 0  |    1  |
  |    |       |
  x    v  ----->
~~~~

You can print out the following six images on a good color printer, glue
them onto sturdy cardboard squares and tape them together to form a
cube, to make a model S2Cell-Hierarchy Globe, which will let you
visually locate any cell at level 5 or higher.

CAVEAT: While these images are conceptually useful, the satellite imagery has
not been transformed correctly, so don't rely on these images to determine
which geographic features are contained by a given `S2CellId`.

[![](img/face0_disp.jpg)](img/face0.jpg)

[![](img/face1_disp.jpg)](img/face1.jpg)

[![](img/face2_disp.jpg)](img/face2.jpg)

[![](img/face3_disp.jpg)](img/face3.jpg)

[![](img/face4_disp.jpg)](img/face4.jpg)

[![](img/face5_disp.jpg)](img/face5.jpg)
