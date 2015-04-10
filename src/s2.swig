// Copyright 2006 Google Inc. All Rights Reserved.

%include base/swig/google.swig

%{
#include <strstream>

#include "s2cellid.h"
#include "s2region.h"
#include "s2cap.h"
#include "s2latlng.h"
#include "s2latlngrect.h"
#include "s2loop.h"
#include "s2polygon.h"
#include "s2polyline.h"
#include "s2regioncoverer.h"
#include "s2cell.h"
#include "s2cellunion.h"
%}

// The PACKED macro makes SWIG think that we're declaring a variable of type
// S2CellId named PACKED.  We don't need it so we clobber it with an empty body.
#define PACKED

// The DECLARE_POD macro makes SWIG think the specified type's defined twice.
// We don't need it, so we can just remove it.
#define DECLARE_POD(TypeName)

// Warning 510 is "friend function 'operator +' ignored.". We can't do anything
// about that.  Warning 389 is 'operator[] ignored' (these are alternate
// accessor methods on R1Interval and S1Interval), which is fine.
#pragma SWIG nowarn=510,389

// If we don't ignore this, the wrapper ends up assigning to None
%ignore S2CellId::None;

#ifdef SWIGPYTHON

%inline %{
  static PyObject *FromS2CellId(const S2CellId &cell_id) {
    return SWIG_NewPointerObj(new S2CellId(cell_id), SWIGTYPE_p_S2CellId, 1);
  }
%}

%typemap(in, numinputs=0)
std::vector<S2CellId> *OUTPUT(std::vector<S2CellId> temp) {
  $1 = &temp;
}

%typemap(argout, fragment="t_output_helper")
std::vector<S2CellId> *OUTPUT {
  $result = t_output_helper($result, vector_output_helper($1, &FromS2CellId));
}

%apply std::vector<S2CellId> *OUTPUT {std::vector<S2CellId> *covering};
%apply std::vector<S2CellId> *OUTPUT {std::vector<S2CellId> *output};

%{
template<class T, class TR>
static PyObject* array_output_helper(const T* arr, int size,
                                     PyObject* (*converter)(const TR x)) {
  if (arr == NULL) {
    // Return a nice out-of-band value if the pointer is NULL.
    Py_INCREF(Py_None);
    return Py_None;
  }

  PyObject* const lst = PyList_New(size);
  if (lst == NULL)
    return NULL;

  for (int i = 0; i < size; i++) {
    PyObject* const obj = converter(arr[i]);
    if (!obj) {
      Py_DECREF(lst);
      return NULL;
    }
    PyList_SET_ITEM(lst, i, obj);
  }

  return lst;
}
%}

%typemap(in, numinputs=0)
S2CellId *OUTPUT_ARRAY_4(S2CellId temp[4]) {
  $1 = temp;
}

%typemap(argout, fragment="t_output_helper")
S2CellId *OUTPUT_ARRAY_4 {
  $result = t_output_helper($result, array_output_helper($1, 4, &FromS2CellId));
}

%apply S2CellId *OUTPUT_ARRAY_4 {S2CellId neighbors[4]};

// This overload shadows the one the takes vector<uint64>&, and it
// does not work anyway.
%ignore S2CellUnion::Init(vector<S2CellId> const& cell_ids);

// The SWIG code which picks between overloaded methods doesn't work
// when given a list parameter.  SWIG_Python_ConvertPtrAndOwn calls
// SWIG_Python_GetSwigThis, doesn't find the 'this' attribute and gives up.
// To avoid this problem rename the Polyline::Init methods so they aren't
// overloaded.
%rename(InitFromS2LatLngs) S2Polyline::Init(vector<S2LatLng> const& vertices);
%rename(InitFromS2Points) S2Polyline::Init(vector<S2Point> const& vertices);

%apply int *OUTPUT {int *next_vertex};

// Avoid leaking memory when holding S2Points in python.
class S2Point {
 public:
  ~S2Point();
};

// The extensions below exist because of the difficulty swigging S2Point.

// This alternate method of S2Loop::vertex() returns a S2LatLng instead.
%extend S2Loop {
 public:
  S2LatLng GetS2LatLngVertex(int i) {
    return S2LatLng(self->vertex(i));
  }
};

// This alternate method of S2Cell::GetVertex() returns a S2LatLng instead.
%extend S2Cell {
 public:
  S2LatLng GetS2LatLngVertex(int k) {
    return S2LatLng(self->GetVertex(k));
  }
};

// This alternate method of S2Cell::GetEdge() returns a S2LatLng instead.
%extend S2Cell {
 public:
  S2LatLng GetS2LatLngEdge(int k) {
    return S2LatLng(self->GetEdge(k));
  }
};

%template() std::vector<S2Point>;
%template() std::vector<S2LatLng>;

#endif

%include "util/geometry/r1interval.h"
%include "util/geometry/s1angle.h"
%include "util/geometry/s1interval.h"
%include "util/geometry/s2cellid.h"
%include "util/geometry/s2region.h"
%include "util/geometry/s2cap.h"
%include "util/geometry/s2latlng.h"
%include "util/geometry/s2latlngrect.h"
%include "util/geometry/s2loop.h"
%include "util/geometry/s2polygon.h"
%include "util/geometry/s2polyline.h"
%include "util/geometry/s2regioncoverer.h"
%include "util/geometry/s2cell.h"
%include "util/geometry/s2cellunion.h"

%define USE_STREAM_INSERTOR_FOR_STR(type)
  %extend type {
    string __str__() {
      std::ostrstream output;
      output << *self << std::ends;
      return output.str();
    }
  }
%enddef

%define USE_EQUALS_FOR_EQ_AND_NE(type)
  %extend type {
    bool __eq__(const type& other) {
      return *self == other;
    }

    bool __ne__(const type& other) {
      return *self != other;
    }
  }
%enddef

%define USE_COMPARISON_FOR_LT_AND_GT(type)
  %extend type {
    bool __lt__(const type& other) {
      return *self < other;
    }

    bool __gt__(const type& other) {
      return *self > other;
    }
  }
%enddef

%define USE_STD_HASH_FOR_HASH(type)
  %extend type {
    size_t __hash__() {
      return HASH_NAMESPACE::hash<type>()(*self);
    }
  }
%enddef

USE_STREAM_INSERTOR_FOR_STR(S1Angle)
USE_STREAM_INSERTOR_FOR_STR(S1Interval)
USE_STREAM_INSERTOR_FOR_STR(S2CellId)
USE_STREAM_INSERTOR_FOR_STR(S2Cap)
USE_STREAM_INSERTOR_FOR_STR(S2LatLng)
USE_STREAM_INSERTOR_FOR_STR(S2LatLngRect)

USE_EQUALS_FOR_EQ_AND_NE(S2CellId)
USE_COMPARISON_FOR_LT_AND_GT(S2CellId)
USE_STD_HASH_FOR_HASH(S2CellId)
