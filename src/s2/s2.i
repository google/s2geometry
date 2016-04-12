%include "std_vector.i"
%include "std_string.i"

%template() std::vector<unsigned long long>;
%template() std::vector<S2CellId>;
%template() std::vector<S2Point>;
%template() std::vector<S2LatLng>;

%apply int {int32};
%apply unsigned long long {uint64};
%apply std::string {string};
%apply std::vector<unsigned long long> const & {std::vector<uint64> const &};

// Standard Google convention is to ignore all functions and methods, and
// selectively add back those for which wrapping is both required and
// functional.
%define %ignoreall %ignore ""; %enddef
%define %unignore %rename("%s") %enddef
%define %unignoreall %rename("%s") ""; %enddef

%include "s2_common.i"
