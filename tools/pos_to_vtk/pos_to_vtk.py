from cell import CellType
from cells import Cells
import click
import pyvista as pv

# This defines the main entry point for the command-line utility. This sets up
# all arguments and options following Click's specifications.


@click.command()
@click.argument('RBC', nargs=-1, type=click.Path(exists=True, allow_dash=True))
@click.option('--plt', 'plt', multiple=True, type=click.Path(exists=True),
              help="""Path for platelet position PLT files, typically
              'PLT.pos'. This argument can be supplied multiple times if
              multiple PLT files are to be considered.""")
@click.option('-v', '--verbose', 'verbose', is_flag=True,
              help="""Increase verbosity. This prints the file path to stdout
              when '-o' or '--output' is used.""")
@click.option('-o', '--output', 'output', type=click.Path(writable=True),
              help="""Write data to the output file. This combines all
              'RBC.pos' and 'PLT.pos' files into a single file. The output
              format is inferred from the output file's extension as long as it
              is supported through the `meshio' package. For large cell counts,
              using the XDMF format is suggested.""")
@click.option('--show-axes/--hide-axes', 'axes', is_flag=True,
              default=True, help="""Show xyz axes orientation with labels.""")
@click.option('--bounding-box', 'bounding_box', nargs=3, type=float,
              help="""Specify a bounding box by its `x, y, z` upper bounds,
              such that the drawn box is defined with bounds=(0 x 0 y 0 z).""")
@click.option('--stl', 'stl',
              help="""Draws the wireframe of the provided surface mesh,
              commonly given in the surface mesh although other surface meshes
              can be provided if `PyVista` provides a corresponding reader
              through its `get_reader` method. For STL files the unit of
              millimetre (mm) is assumed.""")
@click.option('--clip', 'clip', is_flag=True,
              help="""Clips all cells that fall outside the given bounding box
              (`--bounding-box`) or STL surface mesh (`--stl`). The
              visualisation and output will only considers cells that
              are fully enclosed by the bounding box or surface mesh.""")
def cli(rbc, plt, verbose, output, axes, bounding_box, clip, stl):
    """Convert and visualise cell position files used in HemoCell [0].

    The script reads from any number of red blood cell position files (RBC.pos)
    and allows for direct visualisation with PyVista [1] or conversion to XDMF
    for inspection with Paraview through `-o, --output`.

    The RBC positions can be provided through `stdin` by specifying a dash
    (`-`) as argument and piping the input, e.g. `cat RBC.pos | pos_to_vtk -`.

    Optionally, any number of platelets position files (PLT.pos) can be
    supplied by the `--plt` argument, which might be repeated any number of
    times to specify multiple PLT.pos files.

    [0] https://hemocell.eu

    [1] https://dev.pyvista.org/index.html
    """
    if len(rbc) == 0 and len(plt) == 0:
        message = """No input files provided. When passing an RBC through
`stdin`, please provide a dash (`-`) as the argument."""
        raise click.UsageError(message)

    if (bounding_box is not None) and (stl is not None):
        message = """Options `--bounding-box` and `--stl` are exclusive and
cannot be combined. Please provide only one."""
        raise click.UsageError(message)

    rbcs = Cells(flatten([Cells.from_pos(rbc) for rbc in rbc]))
    plts = Cells(flatten([Cells.from_pos(plt) for plt in plt]))

    if len(rbcs) == 0:
        # The script terminates early when no RBCs are present. This
        # most-likely corresponds with users error, e.g. by providing the wrong
        # file path or an empty file.
        raise click.UsageError("No RBC cells to display.""")

    if bounding_box:
        bounds = flatten(zip((0, 0, 0), bounding_box))
        wireframe = pv.Box(bounds)

    if stl:
        wireframe = pv.get_reader(stl).read()
        wireframe.scale((1e3, 1e3, 1e3))

    if clip:
        rbcs = rbcs.clip(wireframe, bounding_box)
        plts = plts.clip(wireframe, bounding_box)

    if (plt and len(plt) == 0) and len(rbcs) == 0:
        raise click.UsageError("No cells remain after clipping with ...")

    if output:
        rbcs = create_glyphs(rbcs, CellType.RBC)
        merged_cells = rbcs.merge(plts) if plt else rbcs
        pv.save_meshio(output, merged_cells)

        if verbose:
            click.echo(f"Written to: '{click.format_filename(output)}'")
        return 0

    plotter = pv.Plotter()

    rbcs_glyph = create_glyphs(rbcs, CellType.RBC)
    rbcs_actor = plotter.add_mesh(rbcs_glyph, color="red")

    if plt:
        plts_glyph = create_glyphs(plts, CellType.PLT)
        plts_actor = plotter.add_mesh(plts_glyph, color="yellow")

    if bounding_box or stl:
        plotter.add_mesh(wireframe, style='wireframe')

    message = f'RBC: {len(rbcs)}\nPLT: {len(plts)}\n'

    if clip:
        # Only when the domain is clipped, a bounding domain is known either
        # through the given bounding box or by the given STL surface mesh.
        # These allow to determine an estimate of the domain's volume and
        # thereby estimate the expected hematocrit given the number of RBCs
        # left after clipping.
        volume = 1e-6 * wireframe.volume
        hc_percentage = 100 * rbcs.hematocrit(wireframe)
        message += f'Domain volume estimate (µm^3): {volume:.2f}\n'
        message += f'Hematocrit (vol%): {hc_percentage:.2f}\n'

    plotter.add_text(message, position='lower_right')

    # A simple widget illustrating the x-y-z axes orientation of the domain.
    if axes:
        plotter.add_axes(line_width=5, labels_off=False)

    # Two radio buttons are added that enable to show/hide the glyphs
    # corresponding to all RBCs or PLTs, where the lambda functions are given
    # to implement the callback function toggling the right set of glyphs.
    plotter.add_checkbox_button_widget(lambda x: rbcs_actor.SetVisibility(x),
                                       position=(5, 12),
                                       color_on="red",
                                       value=True)
    plotter.add_checkbox_button_widget(lambda x: plts_actor.SetVisibility(x),
                                       position=(5, 62),
                                       color_on="yellow",
                                       value=True)

    return plotter.show()


def flatten(list_of_lists):
    """Flatten a list of lists into a single consecutive list."""
    return [val for sublist in list_of_lists for val in sublist]


def create_point_cloud(cells, scaling=1):
    """Create a point cloud `pv.Polydata` of a set of pos files.

    The cells are defined at their centers with a directional unit vector given
    their orientation. The scaling is used to change the magnitude of the
    glyphs to visualise the difference between red blood cells and platelets.
    """
    points = pv.PolyData([cell.center() for cell in cells])
    points["vectors"] = [cell.direction() for cell in cells]
    points["scalars"] = [scaling for cell in cells]
    points.set_active_vectors("vectors")
    return points


def create_glyphs(cell_files, cell_type):
    """Create a glyph-based representation of the cell files."""
    point_cloud = create_point_cloud(cell_files, cell_type.scaling())

    # Only draw the progress bar when there is sufficient work to be done.
    draw_progress_bar = point_cloud.n_cells > 25_000

    glyphs = point_cloud.glyph(
        geom=cell_type.geometry(),
        orient="vectors",
        scale="scalars",
        progress_bar=draw_progress_bar)

    return glyphs
