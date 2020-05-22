This project contains copies of the following files from HEALPix
(http://healpix.sourceforge.net):

cextern/chealpix/chealpix.{h,c}
    Copied from https://sourceforge.net/p/healpix/code/HEAD/tree/trunk/src/C/subs/

These files are licensed under the terms of the following license.

---

| !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
| 
|   Copyright (C) 1997-2016  HEALPix collaboration
| 
| 
|   HEALPix is free software; you can redistribute it and/or modify
|   it under the terms of the GNU General Public License as published by
|   the Free Software Foundation; either version 2 of the License, or
|   (at your option) any later version.
| 
|   HEALPix is distributed in the hope that it will be useful,
|   but WITHOUT ANY WARRANTY; without even the implied warranty of
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
|   GNU General Public License for more details.
| 
|   You should have received a copy of the GNU General Public License
|   along with HEALPix; if not, write to the Free Software
|   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
| 
|   For more information about HEALPix see http://healpix.sourceforge.net
| 
| !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
| 
| HISTORY
| 
|  This is a software package for HEALPix, the Hierarchical, Equal Area, and 
|  iso-Latitude Pixelation of the sphere, which was  originally devised 
|  in early 1997 by Krzysztof M. Gorski at the Theoretical Astrophysics 
|  Center in Copenhagen. Early development of the HEALPix concept and initial 
|  implementation was made in the spring of 1997 in  collaborative work 
|  of KMG, currently at JPL with Eric Hivon, currently at IAP. Benjamin D. Wandelt,
|  currently at IAP, has contributed critically to the 
|  further development of HEALPix and related mathematical methods.
| 
|  The  HEALPix-based  methods  contained  in  the  present  package,  their 
|  implementation in Fortran90, IDL, C, C++ and java routines and codes, and documentation 
|  were developed jointly by: 
|  Anthony J. Banday (banday AT mpa-garching.mpg.de),
|  Matthias Bartelmann (mbartelmann AT ita.uni-heidelberg.de),
|  Hans Kristian Eriksen (h.k.k.eriksen AT astro.uio.no),
|  Krzysztof M. Gorski (krzysztof.m.gorski AT jpl.nasa.gov),
|  Frode K. Hansen (fhansen AT mpa-garching.mpg.de), 
|  Eric Hivon (hivon AT iap.fr),
|  William O'Mullane (womullan AT skysrv.pha.jhu.edu),
|  Martin Reinecke (martin AT mpa-garching.mpg.de), and
|  Benjamin D. Wandelt (bwandelt AT iap.fr).
| 
|  The HEALPix package is available at this address: 
|  http://healpix.sourceforge.net
| 
| COPYRIGHTS
| 
|  Copyright 1997 by Krzysztof M. Gorski
|  for the HEALPix design and logo.
|   All rights reserved.
| 
|  Copyright 1997 by Eric Hivon and Krzysztof M. Gorski
|  for the HEALPix package version 0.80
|   All rights reserved.
| 
|  Copyright 1998 by Benjamin Wandelt, Eric Hivon, and Krzysztof M. Gorski
|  for the HotSpot package
|   All rights reserved.
| 
|  Copyright 1998 by K.M. Gorski, Eric Hivon, and B.D. Wandelt
|  for the HEALPix package version 0.90
|   All rights reserved.
| 
|  Copyright 1998 by A.J. Banday and M.S. Bartelmann
|  for the healgif package
|   All rights reserved.
| 
|  Copyright 1999 by K.M. Gorski, B.D. Wandelt, E. Hivon, F. Hansen, A.J. Banday,
|  and M. Bartelmann
|  for the HEALPix package version 1.00
|   All rights reserved.
| 
|  Copyright 2000 by K.M. Gorski, B.D. Wandelt, E. Hivon, F. Hansen, A.J. Banday,
|  and M. Bartelmann
|  for the HEALPix package version 1.10
|   All rights reserved.
| 
|  Copyright 2003 by K.M. Gorski, E. Hivon, M. Reinecke, A.J. Banday,
|  B.D. Wandelt, and M. Bartelmann
|  for the HEALPix package version 1.20
|   All rights reserved.
| 
|  Copyright 2003 by E. Hivon, Kenneth M. Ganga (ganga AT cdf.in2p3.fr), Reza Ansari (ansari AT lal.in2p3.fr)
|  for the HEALPix C package
|   All rights reserved.
| 
|  Copyright 2005 by E. Hivon, M. Reinecke, W. O'Mullane, H.K. Eriksen, K.M. Gorski, A.J. Banday
|  for the HEALPix package version 2.00
|   All rights reserved.
| 
|  Copyright 2005 by  H.K. Eriksen, Snorre Boasson (snorre.boasson AT ntnu.no), M. Reinecke, E. Hivon
|  for the MPI implementation of HEALPix F90 alm routines
|   All rights reserved.
| 
|  Copyright 2008 by E. Hivon, M. Reinecke, W. O'Mullane, H.K. Eriksen, K.M. Gorski, A.J. Banday, E. Joliet
|  for the HEALPix package version 2.10
|   All rights reserved.
| 
|  Copyright 2008 by G. Rocha
|  for the ngsims package
|   All rights reserved.
| 
|  Copyright 2008 by J. P. Leahy
|  for the ximview package
|   All rights reserved.
| 
|  Copyright 2008 by D. Larson
|  for the alice package
|   All rights reserved.
| 
|  Copyright 2010 by E. Hivon, M. Reinecke, W. O'Mullane, E. Joliet, K.M. Gorski, A.J. Banday
|  for the HEALPix package version 2.20
|   All rights reserved.
| 
|  Copyright 2012 by E. Hivon, M. Reinecke, W. O'Mullane, E. Joliet, K.M. Gorski, A.J. Banday
|  for the HEALPix package version 3.00
|   All rights reserved.
| 
|  Copyright 2013 by E. Hivon, M. Reinecke, W. O'Mullane, E. Joliet, K.M. Gorski, A.J. Banday
|  for the HEALPix package version 3.10
|   All rights reserved.
| 
|  Copyright 2013 by A. Zonca, C. Rosset
|  for the healpy package version 1.5.0
|   All rights reserved.
| 
|  Copyright 2014 by E. Hivon, M. Reinecke, W. O'Mullane, E. Joliet, K.M. Gorski, A.J. Banday
|  for the HEALPix package version 3.20
|   All rights reserved.
| 
|  Copyright 2014 by A. Zonca, L. Singer, C. Rosset
|  for the healpy package version 1.8.1
|   All rights reserved.
| 
|  Copyright 2015 by E. Hivon, M. Reinecke, W. O'Mullane, E. Joliet, K.M. Gorski, A.J. Banday
|  for the HEALPix package version 3.30
|   All rights reserved.
| 
|  Copyright 2015 by A. Zonca, L. Singer, C. Rosset
|  for the healpy package version 1.9.0
|   All rights reserved.
| 
|  Copyright 2016 by E. Hivon, M. Reinecke, W. O'Mullane, E. Joliet, K.M. Gorski, A.J. Banday
|  for the HEALPix package version 3.31
|   All rights reserved.
| 
|  Copyright 2015 by A. Zonca, L. Singer, C. Rosset
|  for the healpy package version 1.9.1
|   All rights reserved.
| 
| !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
| 
| ACKNOWLEDGMENTS
| 
| See http://healpix.sourceforge.net/credits.php
