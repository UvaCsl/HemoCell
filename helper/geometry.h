#include "palabos3D.h"
#include "palabos3D.hh"
#include <cmath>
#include <limits>

#ifndef GEOMETRY_H
#define GEOMETRY_H

namespace geom {

/// \brief A domain functional indicating true inside a rectangular box.
/// The BoxDomain functional evaluates to true inside its box, where inclusive
/// edges are considered using plb::contained().
class BoxDomain : public plb::DomainFunctional3D {
public:
  explicit BoxDomain(plb::Box3D box) : box{box} {}

  /// Returns true inside the rectangular box (inclusive edges).
  bool operator()(plb::plint x, plb::plint y, plb::plint z) const override {
    return plb::contained(x, y, z, box);
  }

  /// The relevant bounding box of the BoxDomain functional is known.
  const plb::Box3D& bounding_box() const { return box; }

  BoxDomain *clone() const override { return new BoxDomain(*this); }
private:
  plb::Box3D box;
};

/// \brief An ellipsoidal domain functional returning true inside the ellipsoid.
/// The EllipseDomain represents an ellipsoidal domain which evaluates to true
/// inside the ellipsoid. The ellipsoid can be given a desired center location
/// and radius for each dimension. However, the EllipseDomain does current not
/// support rotation.
class EllipseDomain : public plb::DomainFunctional3D {
public:
  /// The ellipsoid can be constructed with only x,y coordinates and x,y radii.
  /// Doing so initialises and ellipsoid at z=0 with infinite radius in the z
  /// direction. This mimics a two-dimensional projection and ensures the
  /// ellipsoidal domain is constant in z-direction.
  EllipseDomain(double x, double y, double rx, double ry)
      : EllipseDomain(x, y, 0, rx, ry, std::numeric_limits<double>::max())
  {}

  /// Default constructor specifying the x, y, z coordinates and radii.
  EllipseDomain(double cx, double cy, double cz, double rx, double ry, double rz)
      : cx{ cx },
        cy{ cy },
        cz{ cz },
        rx{ rx },
        ry{ ry },
        rz{ rz }
  {}

  /// Returns true if a point (x, y, z) is contained within the ellipsoid. This
  /// evaluates the parametric equation: x^2/a^2 + y^2/b^2 + z^2/c^2 = 1, where
  /// the equation is translated to the center of the ellipsoid.
  bool operator()(plb::plint x, plb::plint y, plb::plint z) const override {
    auto dx = x - cx;
    auto dy = y - cy;
    auto dz = z - cz;
    return (dx*dx)/(rx*rx) + (dy*dy)/(ry*ry) + (dz*dz)/(rz*rz) < 1;
  }

  EllipseDomain *clone() const override { return new EllipseDomain(*this); }

private:
  double cx = 0;
  double cy = 0;
  double cz = 0;
  double rx = 0;
  double ry = 0;
  double rz = 0;
};

/// \brief The union of multiple domain functionals.
/// This domain functional builds the union of any number of domain functionals
/// (derived from plb::DomainFunctional3D). The union of these domains considers
/// the region where any of the included domains returns true.
class Union : public plb::DomainFunctional3D {
public:
  Union(std::initializer_list<std::reference_wrapper<plb::DomainFunctional3D>>
            domains)
      : domains{domains}
  {}

  /// Return true if any of the included domains is true.
  bool operator()(plb::plint x, plb::plint y, plb::plint z) const override {
    auto predicate = [x, y, z](const plb::DomainFunctional3D &domain) {
      return domain(x, y, z);
    };
    return std::any_of(domains.begin(), domains.end(), predicate);
  }

  Union *clone() const override { return new Union(*this); };

private:
  std::vector<std::reference_wrapper<plb::DomainFunctional3D>> domains;
};

/// \brief Create the boolean difference of the first and subsequent domains.
/// This creates the boolean difference between the first domain with respect to
/// any of the subsequent domains. The difference is considered true if the
/// first domain is true and _none of_ the other domains are true.
///
/// NOTE: when only a single domain is provided, the difference returns false.
class Difference : public plb::DomainFunctional3D {
public:
  Difference(
      std::initializer_list<std::reference_wrapper<plb::DomainFunctional3D>> domains)
      : domains{ domains } {}

  /// Return true if the first domain is true and none of the others.
  bool operator()(plb::plint x, plb::plint y, plb::plint z) const override {
    if (domains.empty())
      return false;

    if (!domains.front()(x, y, z))
      // No need to subtract subsequent domains if the first domain is empty.
      return false;

    auto predicate = [x, y, z](const plb::DomainFunctional3D& domain) {
      return domain(x, y, z);
    };
    return std::none_of(std::next(domains.begin()), domains.end(), predicate);
  }

  Difference *clone() const override { return new Difference(*this); }

private:
  std::vector<std::reference_wrapper<plb::DomainFunctional3D>> domains;
};

/// \brief The boolean intersection of a number of domain functionals.
/// This creates the boolean intersection between the provided domains. The
/// intersection returns true when _all of_ the domains are true.
///
/// NOTE: when no domains are provided, the intersection always returns true.
class Intersection : public plb::DomainFunctional3D {
public:
  Intersection(
      std::initializer_list<std::reference_wrapper<plb::DomainFunctional3D>>
          domains)
      : domains{domains} {}

  /// Returns true only if all domains are true.
  bool operator()(plb::plint x, plb::plint y, plb::plint z) const override {
    if (domains.empty())
      // No domains have an empty intersection, whereas std::all_of returns
      // true for emtpy ranges.
      return false;

    auto predicate = [x, y, z](const plb::DomainFunctional3D &domain) {
      return domain(x, y, z);
    };
    return std::all_of(domains.begin(), domains.end(), predicate);
  }

  Intersection *clone() const override { return new Intersection(*this); }

private:
  std::vector<std::reference_wrapper<plb::DomainFunctional3D>> domains;
};

} // namespace geom

#endif /* GEOMETRY_H */
