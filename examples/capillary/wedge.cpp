#include "wedge.h"

///
/// \brief Geometry definition for a constricted capillary channel.
///
/// Defines a rectangular channel with a wedge-shaped constriction in the
/// channel. The channel is parametric in its length, height, wedge length, and
/// gap size between both opposing wedges. The gap is symmetric, where the
/// height of the top and bottom wedges are equal. Each wedge reduces half its
/// height along the length of the channel, as given by its slope parameter.
///
/// ---------------------------------------------------------------------------
///                         |                  |
///                         |        |----------
///                         |--------
///
///                         |--------
///                         |        |----------
///                         |                  |
/// ---------------------------------------------------------------------------
///
/// \tparam T
template <typename T>
class TriangleShapeDomain3D : public plb::DomainFunctional3D {
public:
  TriangleShapeDomain3D(unsigned channel_length, unsigned channel_height,
                        unsigned wedge_length, unsigned gap_size) {
    bottom_height = (channel_height - gap_size) / 2;
    top_height = channel_height - 1 - bottom_height;
    start_wedge = (channel_length - wedge_length) / 2;
    end_wedge = start_wedge + wedge_length;
    // half triangle reduction along triangle length
    slope = T(top_height) / 2.0 / T(wedge_length);
  }

  bool operator()(plint x, plint y, plint) const override {
    auto dy = slope * (x - start_wedge);
    return (x > start_wedge && x <= end_wedge) &&
           (y <= bottom_height - dy || y >= top_height + dy);
  }

  TriangleShapeDomain3D<T> *clone() const override {
    return new TriangleShapeDomain3D<T>(*this);
  }

private:
  plint start_wedge;
  plint end_wedge;
  plint top_height;
  plint bottom_height;
  T slope;
};

/// \brief Length `refDirN` with square cross section of 36x36 cells.
std::tuple<unsigned, unsigned, unsigned> Wedge::domain_size(unsigned resolution) {
  return std::make_tuple(resolution, 36, 36);
}

int Wedge::geometry(plb::MultiBlockLattice3D<double, DESCRIPTOR> *&lattice, unsigned resolution) {
  auto nx = 0;
  auto ny = 0;
  auto nz = 0;
  std::tie(nx, ny, nz) = Wedge::domain_size(resolution);

  unsigned wedge_length = 50;
  unsigned gap_size = 12;

  defineDynamics(*lattice, lattice->getBoundingBox(),
                 new TriangleShapeDomain3D<T>(nx, ny, wedge_length, gap_size),
                 new plb::BounceBack<T, DESCRIPTOR>);

  return 0;
}

plb::Array<double, 3> Wedge::driving_force(double force) {
  auto scaling = hemo::param::dx * hemo::param::dx * hemo::param::dt *
                 hemo::param::dt / hemo::param::dm;
  return {force * scaling, 0, 0};
}
