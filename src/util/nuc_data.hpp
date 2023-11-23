#ifndef NUC_DATA_HEADER
#define NUC_DATA_HEADER

#include "util/constants.hpp"
#include "util/types.hpp"

namespace osiris {

struct Isotope {
  int A;
  int Z;
  // TODO load a mass table
  real mass;
};

/// @brief incident particle type
enum class Proj : bool {
  proton,
  neutron,
};

template <Proj p> constexpr real mass() {
  if constexpr (p == Proj::neutron)
    return constants::n_mass_amu;
  if constexpr (p == Proj::proton)
    return constants::p_mass_amu;
}

template <Proj p> constexpr int charge() {
  if constexpr (p == Proj::neutron)
    return 0;
  if constexpr (p == Proj::proton)
    return 1;
}

template <Proj p>
/// @brief 2j+1, where j is the intrinsic spin of the projectile
constexpr int spin() {
  if constexpr (p == Proj::neutron)
    return 2;
  if constexpr (p == Proj::proton)
    return 2;
}

} // namespace osiris

#endif
