#include "bifurcation.h"

std::tuple<unsigned, unsigned, unsigned>
Bifurcation::domain_size(unsigned resolution) {
  return std::make_tuple(8 * resolution, resolution, resolution);
}

/// \brief Geometry definition for the bifurcating capillary channel.
/// The geometry is defined as a single channel with periodicity along the
/// x-direction. The channel splits into two branches, where an inner dividing
/// solid region is defined.
int Bifurcation::geometry(plb::MultiBlockLattice3D<double, DESCRIPTOR> *&lattice,
                    unsigned resolution, double capillary_diameter) {
  auto nx = 0;
  auto ny = 0;
  auto nz = 0;
  std::tie(nx, ny, nz) = Bifurcation::domain_size(resolution);

  unsigned wall = 2;

  // The dimensions of the inner ellipse is derived from the outer ellipse.
  auto outer_rx = ny - 2 * wall;
  auto outer_ry = 0.5 * outer_rx;
  auto inner_ry = outer_ry - capillary_diameter;
  auto inner_rx = outer_rx * inner_ry / outer_ry;

  if (inner_ry < 1 || inner_rx < 1) {
    plb::pcout << "Error generating geometry Bifurcation capillary: "
               << "zero radii encountered for inner ellipse."
               << "Found major, minor axis: " << inner_rx << ", " << inner_ry
               << std::endl;
    return -1;
  }

  // The center coordinates of the left-most ellipse.
  auto cx = 0.1875 * nx;
  auto cy = ny * 0.5 - 1;

  // Find the intersection between the x-coordinate bound of the ellipses and
  // the outer ellipse to determine the size of the inlet channel. This solves
  // the equation given by the quadratic ellipse equation and a line:
  // a*x**2 + b*x + c = 0, where a = 1.0 (and is therefore left out).
  auto ellipse_start = cx - outer_rx + capillary_diameter;
  auto b = -2.0 * cy;
  auto c = cy * cy - outer_ry * outer_ry *
                         (1 - (ellipse_start - cx) * (ellipse_start - cx) /
                                  (outer_rx * outer_rx));
  auto d = b * b - 4.0 * c;

  if (d < 0.0) {
    plb::pcout
        << "(Capillary) Negative  discriminant: Check the geometry setup!"
        << std::endl;
    return -1;
  }

  unsigned y_top = 0;
  unsigned y_bot = 0;

  if (d == 0) {
    y_top = -b / 2.0;
    y_bot = y_top;
  } else {
    y_top = std::ceil((-b + std::sqrt(d)) / (2.0));
    y_bot = std::floor((-b - std::sqrt(d)) / (2.0));
  }

  // The outer region: This can be expressed by the difference of the full
  // simulation domain and the inner regions, i.e. the center divider and the
  // outer walls near the inlet/outlet and center regions.
  auto bounding_box = geom::BoxDomain(lattice->getBoundingBox());
  auto inlet = geom::BoxDomain(plb::Box3D(0, nx, y_bot - 1, y_top - 1, 0, nz));
  auto left_outer_ellipse = geom::EllipseDomain(cx, cy, outer_rx, outer_ry);
  auto right_outer_ellipse = geom::EllipseDomain(nx - cx, cy, outer_rx, outer_ry);
  auto center = geom::BoxDomain(plb::Box3D(cx, nx - cx, wall, ny - wall - 2, 0, nz));

  defineDynamics(*lattice, lattice->getBoundingBox(),
                 new geom::Difference({bounding_box, inlet, left_outer_ellipse,
                                       right_outer_ellipse, center}),
                 new plb::BounceBack<T, DESCRIPTOR>());

  // This includes the left and right inner ellipsoids together with a
  // rectangular block to connect both ellipsoids.
  auto middle = geom::BoxDomain(plb::Box3D(cx, nx - cx, 
                                           wall + capillary_diameter - 1,
                                           ny - wall - capillary_diameter - 1, 0, nz));
  auto left_ellipse = geom::EllipseDomain(cx, cy, inner_rx, inner_ry);
  auto right_ellipse = geom::EllipseDomain(nx - cx, cy, inner_rx, inner_ry);

  defineDynamics(*lattice, lattice->getBoundingBox(),
                 new geom::Union({left_ellipse, middle, right_ellipse}),
                 new plb::BounceBack<T, DESCRIPTOR>());
  return 0;
}

plb::Array<double, 3> Bifurcation::driving_force(double force) {
  auto scaling = hemo::param::dx * hemo::param::dx * hemo::param::dt *
                 hemo::param::dt / hemo::param::dm;
  return {force * scaling, 0, 0};
}
