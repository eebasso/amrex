#ifndef AMREX_FIELD_LOCATION_H_
#define AMREX_FIELD_LOCATION_H_
#include <AMReX_Config.H>

#include <AMReX_Dim3.H>
#include <AMReX_IndexType.H>


namespace amrex {

enum class FieldLocation { NODE, EDGE, FACE, CELL };

bool
is_field_location_ok (const FieldLocation& loc, const Dim3& ix)
{
    constexpr int n = static_cast<int>(IndexType::NODE);
    constexpr int c = static_cast<int>(IndexType::CELL);

    if        (loc == FieldLocation::NODE) {
        return (ix.x == n && ix.y == n && ix.z == n);
    } else if (loc == FieldLocation::EDGE) {
        return (ix.x == c && ix.y == n && ix.z == n) ||
               (ix.x == n && ix.y == c && ix.z == n) ||
               (ix.x == n && ix.y == n && ix.z == c);
    } else if (loc == FieldLocation::FACE) {
        return (ix.x == n && ix.y == c && ix.z == c) ||
               (ix.x == c && ix.y == n && ix.z == c) ||
               (ix.x == c && ix.y == c && ix.z == n);
    } else if (loc == FieldLocation::CELL) {
        return (ix.x == c && ix.y == c && ix.z == c);
    }
    else {
        amrex::Abort("is_field_location_ok: Unknown FieldLocation type");
        return false;
    }
}

bool
is_field_location_ok (const FieldLocation& loc, const IndexType& ix) {
    constexpr int n = static_cast<int>(IndexType::NODE);
    constexpr int c = static_cast<int>(IndexType::CELL);


}

bool
does_indextype_match_dim3



}



#endif /* AMREX_FIELD_LOCATION_H_ */