// Copyright 2015 Google Inc. All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS-IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//

#ifndef S2_GEOMETRY_S2TEXTFORMAT_H_
#define S2_GEOMETRY_S2TEXTFORMAT_H_

// s2textformat contains a collection of functions for converting
// geometry to and from a human-readable format.  It is mainly
// intended for testing and debugging.  Be aware that the
// human-readable format is *not* designed to preserve the full
// precision of the original object, so it should not be used
// for data storage.

#include <string>
#include <vector>

#include "s2.h"
#include "s2latlngrect.h"

class S2LatLng;
class S2Loop;
class S2Polygon;
class S2Polyline;

namespace s2textformat {

// Given a latitude-longitude coordinate in degrees,
// return a newly allocated point.  Example of the input format:
//     "-20:150"
S2Point MakePoint(string const& str);

// Given a string of one or more latitude-longitude coordinates in degrees,
// return the minimal bounding S2LatLngRect that contains the coordinates.
// Example of the input format:
//     "-20:150"                     // one point
//     "-20:150, -20:151, -19:150"   // three points
S2LatLngRect MakeLatLngRect(string const& str);

// Given a string of latitude-longitude coordinates in degrees,
// return a newly allocated loop.  Example of the input format:
//     "-20:150, 10:-120, 0.123:-170.652"
// The strings "empty" or "full" create an empty or full loop respectively.
S2Loop* MakeLoop(string const& str);

// Similar to MakeLoop(), but returns an S2Polyline rather than an S2Loop.
S2Polyline* MakePolyline(string const& str);

// Given a sequence of loops separated by semicolons, return a newly
// allocated polygon.  Loops are automatically normalized by inverting them
// if necessary so that they enclose at most half of the unit sphere.
// (Historically this was once a requirement of polygon loops.  It also
// hides the problem that if the user thinks of the coordinates as X:Y
// rather than LAT:LNG, it yields a loop with the opposite orientation.)
//
// Examples of the input format:
//     "10:20, 90:0, 20:30"                                  // one loop
//     "10:20, 90:0, 20:30; 5.5:6.5, -90:-180, -15.2:20.3"   // two loops
//     ""       // the empty polygon (consisting of no loops)
//     "full"   // the full polygon (consisting of one full loop)
//     "empty"  // **INVALID** (a polygon consisting of one empty loop)
S2Polygon* MakePolygon(string const& str);

// Like MakePolygon(), except that it does not normalize loops (i.e., it
// gives you exactly what you asked for).
S2Polygon* MakeVerbatimPolygon(string const& str);

// Parse a string in the same format as MakeLatLngRect, and return the
// corresponding vector of S2LatLng points.
void ParseLatLngs(string const& str, std::vector<S2LatLng>* latlngs);

// Parse a string in the same format as MakeLatLngRect, and return the
// corresponding vector of S2Point values.
void ParsePoints(string const& str, std::vector<S2Point>* vertices);

// Convert a point, lat-lng rect, loop, polyline, or polygon to the string
// format above.
string ToString(S2Point const& point);
string ToString(S2LatLngRect const& rect);
string ToString(S2Loop const* loop);
string ToString(S2Polyline const* polyline);
string ToString(S2Polygon const* polygon);
string ToString(std::vector<S2Point> const& points);
string ToString(std::vector<S2LatLng> const& points);

}  // namespace s2textformat

#endif  // S2_GEOMETRY_S2TEXTFORMAT_H_
