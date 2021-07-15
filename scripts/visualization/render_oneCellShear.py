"""Example script to automate Blender rendering for oneCellShear example."""
import bpy
from math import radians
from mathutils import Matrix, Euler, Quaternion
import os

# Input parameters
t_0 = 0000
t_f = 0000
timestep = 2000

x3dpath = "../../examples/oneCellShear/tmp/x3d/"
renderpath = "../../examples/oneCellShear/tmp/renders/"
final_render_filename = "oneCellShear_render"
render_engine = "CYCLES"  # also try 'BLENDER_EEVEE'


# Helper functions ------------------------------------------------------------
def removelighting():
    for light in list(bpy.data.lights):
        bpy.data.lights.remove(light)


def removemeshes():
    for mesh in list(bpy.data.meshes):
        bpy.data.meshes.remove(mesh)


def removecameras():
    for camera in list(bpy.data.cameras):
        bpy.data.cameras.remove(camera)


def removematerials():
    for material in list(bpy.data.materials):
        bpy.data.materials.remove(material)


def remove_duplicate_objects():
    for obj in bpy.data.objects:
        if "." in obj.name:
            bpy.data.objects.remove(obj)


def clearscene():
    bpy.ops.object.select_all(action='SELECT')
    bpy.ops.object.delete(use_global=False, confirm=False)
    removecameras()
    removemeshes()
    removelighting()
    removematerials()
    #removescenes()


def importx3d(filename, objectname):

    #import the x3d scene
    bpy.ops.import_scene.x3d(filepath=filename,
                             filter_glob="*.x3d;*.wrl",
                             axis_forward='Z',
                             axis_up='Y')

    #rename mesh and object
    for i in range(len(bpy.data.meshes)):
        if "IndexedFaceSet" in bpy.data.meshes[i].name:
            bpy.data.meshes[i].name = objectname
    for i in range(len(bpy.data.objects)):
        if "IndexedFaceSet" in bpy.data.objects[i].name:
            bpy.data.objects[i].name = objectname

    #Remove any camera and lighting that may be present in imported x3d scene
    removelighting()
    removecameras()


def translate_all(x=0, y=0, z=0):
    for i in range(len(bpy.data.objects)):
        bpy.data.objects[i].location.x += x
        bpy.data.objects[i].location.y += y
        bpy.data.objects[i].location.z += z


def scale_all(factor=0.5):
    for i in range(len(bpy.data.objects)):
        bpy.data.objects[i].scale[0] *= factor
        bpy.data.objects[i].scale[1] *= factor
        bpy.data.objects[i].scale[2] *= factor


def rotate_all(x=0, y=0, z=0):
    for i in range(len(bpy.data.objects)):
        obj = bpy.data.objects[i]
        # define the rotation

        rot_mat = Matrix.Rotation(
            radians(x), 4, 'X'
        )  # you can also use as axis Y,Z or a custom vector like (x,y,z)

        # decompose world_matrix's components, and from them assemble 4x4 matrices
        orig_loc, orig_rot, orig_scale = obj.matrix_world.decompose()
        orig_loc_mat = Matrix.Translation(orig_loc)
        orig_rot_mat = orig_rot.to_matrix().to_4x4()
        orig_scale_mat = Matrix.Scale(orig_scale[0], 4,
                                      (.1, 0, 0)) @ Matrix.Scale(
                                          orig_scale[1], 4,
                                          (0, 1, 0)) @ Matrix.Scale(
                                              orig_scale[2], 4, (0, 0, 1))

        # assemble the new matrix
        obj.matrix_world = orig_loc_mat @ rot_mat @ orig_rot_mat @ orig_scale_mat

        rot_mat = Matrix.Rotation(
            radians(y), 4, 'Y'
        )  # you can also use as axis Y,Z or a custom vector like (x,y,z)

        # decompose world_matrix's components, and from them assemble 4x4 matrices
        orig_loc, orig_rot, orig_scale = obj.matrix_world.decompose()
        orig_loc_mat = Matrix.Translation(orig_loc)
        orig_rot_mat = orig_rot.to_matrix().to_4x4()
        orig_scale_mat = Matrix.Scale(orig_scale[0], 4,
                                      (.1, 0, 0)) @ Matrix.Scale(
                                          orig_scale[1], 4,
                                          (0, 1, 0)) @ Matrix.Scale(
                                              orig_scale[2], 4, (0, 0, 1))

        # assemble the new matrix
        obj.matrix_world = orig_loc_mat @ rot_mat @ orig_rot_mat @ orig_scale_mat

        rot_mat = Matrix.Rotation(
            radians(z), 4, 'Z'
        )  # you can also use as axis Y,Z or a custom vector like (x,y,z)

        # decompose world_matrix's components, and from them assemble 4x4 matrices
        orig_loc, orig_rot, orig_scale = obj.matrix_world.decompose()
        orig_loc_mat = Matrix.Translation(orig_loc)
        orig_rot_mat = orig_rot.to_matrix().to_4x4()
        orig_scale_mat = Matrix.Scale(orig_scale[0], 4,
                                      (.1, 0, 0)) @ Matrix.Scale(
                                          orig_scale[1], 4,
                                          (0, 1, 0)) @ Matrix.Scale(
                                              orig_scale[2], 4, (0, 0, 1))

        # assemble the new matrix
        obj.matrix_world = orig_loc_mat @ rot_mat @ orig_rot_mat @ orig_scale_mat


def rotate_one(obj_name, x=0, y=0, z=0):

    obj = bpy.data.objects[obj_name]
    # define the rotation

    rot_mat = Matrix.Rotation(
        radians(x), 4,
        'X')  # you can also use as axis Y,Z or a custom vector like (x,y,z)

    # decompose world_matrix's components, and from them assemble 4x4 matrices
    orig_loc, orig_rot, orig_scale = obj.matrix_world.decompose()
    orig_loc_mat = Matrix.Translation(orig_loc)
    orig_rot_mat = orig_rot.to_matrix().to_4x4()
    orig_scale_mat = Matrix.Scale(orig_scale[0], 4, (.1, 0, 0)) @ Matrix.Scale(
        orig_scale[1], 4, (0, 1, 0)) @ Matrix.Scale(orig_scale[2], 4,
                                                    (0, 0, 1))

    # assemble the new matrix
    obj.matrix_world = orig_loc_mat @ rot_mat @ orig_rot_mat @ orig_scale_mat

    rot_mat = Matrix.Rotation(
        radians(y), 4,
        'Y')  # you can also use as axis Y,Z or a custom vector like (x,y,z)

    # decompose world_matrix's components, and from them assemble 4x4 matrices
    orig_loc, orig_rot, orig_scale = obj.matrix_world.decompose()
    orig_loc_mat = Matrix.Translation(orig_loc)
    orig_rot_mat = orig_rot.to_matrix().to_4x4()
    orig_scale_mat = Matrix.Scale(orig_scale[0], 4, (.1, 0, 0)) @ Matrix.Scale(
        orig_scale[1], 4, (0, 1, 0)) @ Matrix.Scale(orig_scale[2], 4,
                                                    (0, 0, 1))
    # assemble the new matrix
    obj.matrix_world = orig_loc_mat @ rot_mat @ orig_rot_mat @ orig_scale_mat
    rot_mat = Matrix.Rotation(
        radians(z), 4,
        'Z')  # you can also use as axis Y,Z or a custom vector like (x,y,z)
    # decompose world_matrix's components, and from them assemble 4x4 matrices
    orig_loc, orig_rot, orig_scale = obj.matrix_world.decompose()
    orig_loc_mat = Matrix.Translation(orig_loc)
    orig_rot_mat = orig_rot.to_matrix().to_4x4()
    orig_scale_mat = Matrix.Scale(orig_scale[0], 4, (.1, 0, 0)) @ Matrix.Scale(
        orig_scale[1], 4, (0, 1, 0)) @ Matrix.Scale(orig_scale[2], 4,
                                                    (0, 0, 1))
    # assemble the new matrix
    obj.matrix_world = orig_loc_mat @ rot_mat @ orig_rot_mat @ orig_scale_mat


def translate_one(obj_name, x=0, y=0, z=0):
    bpy.data.objects[obj_name].location.x += x
    bpy.data.objects[obj_name].location.y += y
    bpy.data.objects[obj_name].location.z += z


def scale_one(obj_name, xfactor=0.5, yfactor=0.5, zfactor=0.5):
    bpy.data.objects[obj_name].scale[0] *= xfactor
    bpy.data.objects[obj_name].scale[1] *= yfactor
    bpy.data.objects[obj_name].scale[2] *= zfactor


def smooth_mesh(mesh_name):
    mesh = bpy.data.meshes[mesh_name]
    for f in mesh.polygons:
        f.use_smooth = True


# Main render loop ------------------------------------------------------------
if not os.path.isdir(renderpath):
    os.mkdir(renderpath)

for time in range(t_0, t_f + timestep, timestep):
    # First clear the scene of everything
    clearscene()

    # set the render engine
    bpy.context.scene.render.engine = render_engine

    # Uncomment to set the number of cores devoted to rendering
    #for scene in bpy.data.scenes:
    #    scene.render.threads_mode = 'FIXED'
    #    scene.render.threads = 18

    # Import the geometries at the correct time step
    importx3d(filename=x3dpath + 'RBC.' + str(time).zfill(12) + ".x3d",
              objectname="RBC")

    # Rotate and scale geometries
    # NOTE: this will strongly depend on the view that you want. Performing any
    # scaling, framing, and object (cells, lights, cameras) placement, are
    # always a bit of a hack through the API compared to performing these
    # operations manually in the GUI.
    rotate_all(y=0, x=-90, z=0)
    scale_all(factor=1.0)
    translate_all(x=30, z=40, y=0)

    # Also play around to get RBC in the correct view
    scale_one("RBC", xfactor=100000, yfactor=100000, zfactor=100000)

    # Remove duplicate cells since HemoCell can save duplicates when cells are
    # being copied from one core to another.
    remove_duplicate_objects()

    # Add a sun
    bpy.ops.object.light_add(type='SUN', align='WORLD', location=(0, 100, -10))
    bpy.data.objects["Sun"].data.energy = 10  #3 #Give it its energy
    # Rotate if you want to give a nice shade
    rotate_one("Sun", x=-65, z=10)

    # Add a Camera and make sure its pointing at the cell
    bpy.ops.object.camera_add(enter_editmode=False,
                              align='WORLD',
                              location=(28.928, 1.0406, 42.655),
                              rotation=(0.036486, -0.0511929, 180.636))
    bpy.context.scene.camera = bpy.data.objects["Camera"]

    # Add material for RBC and smoothing---------------------------------------
    # NOTE: This is all done using NODES. Much easier to understand in GUI
    smooth_mesh("RBC")
    bpy.data.materials.new("RBC")
    bpy.data.meshes["RBC"].materials.append(bpy.data.materials["RBC"])
    bpy.data.materials["RBC"].use_nodes = True
    # Play around with the color
    bpy.data.materials["RBC"].node_tree.nodes["Principled BSDF"].inputs[
        0].default_value = (1, 0.9, 0.847222, 1)
    bpy.data.materials["RBC"].node_tree.nodes["Principled BSDF"].inputs[
        1].default_value = 0.88
    bpy.data.materials["RBC"].node_tree.nodes["Principled BSDF"].inputs[
        3].default_value = (1, 0, 0, 1)
    bpy.data.materials["RBC"].node_tree.nodes["Principled BSDF"].inputs[
        5].default_value = 0.16
    bpy.data.materials["RBC"].node_tree.nodes["Principled BSDF"].inputs[
        7].default_value = 0.67
    # Add some texture to the cell
    bpy.data.materials["RBC"].node_tree.nodes.new(type="ShaderNodeTexNoise")
    bpy.data.materials["RBC"].node_tree.nodes["Noise Texture"].inputs[
        3].default_value = 6
    bpy.data.materials["RBC"].node_tree.nodes["Noise Texture"].inputs[
        2].default_value = 5
    # Link RGB curves to vertex colors
    bpy.data.materials["RBC"].node_tree.links.new(
        bpy.data.materials["RBC"].node_tree.nodes["Noise Texture"].outputs[1],
        bpy.data.materials["RBC"].node_tree.nodes["Material Output"].inputs[2])

    # Render Properties--------------------------------------------------------
    bpy.context.scene.render.film_transparent = True
    bpy.context.scene.render.image_settings.color_mode = 'RGBA'
    bpy.context.scene.render.image_settings.file_format = 'PNG'

    # Uncomment the following to use a white background
    # NOTE: This is quite complicated in Blender
    #switch on nodes and get reference
    #bpy.context.scene.use_nodes = True
    #tree = bpy.context.scene.node_tree
    #clear default nodes
    #for node in tree.nodes:
    #       tree.nodes.remove(node)
    #create input image node
    #alphaover = tree.nodes.new(type="CompositorNodeAlphaOver")
    #composite = tree.nodes.new(type="CompositorNodeComposite")
    #Rlayers = tree.nodes.new(type="CompositorNodeRLayers")
    #bpy.data.scenes["Scene"].node_tree.nodes["Alpha Over"].premul = 1
    #bpy.data.scenes["Scene"].node_tree.nodes["Alpha Over"].inputs[1].default_value = (1, 1, 1, 1)
    #turn off filmic
    #bpy.context.scene.view_settings.view_transform = 'Standard'
    #link nodes
    #links = tree.links
    #link1 = links.new(bpy.data.scenes["Scene"].node_tree.nodes["Render Layers"].outputs[0],
    #    bpy.data.scenes["Scene"].node_tree.nodes["Alpha Over"].inputs[2])
    #linl2 = links.new(bpy.data.scenes["Scene"].node_tree.nodes["Alpha Over"].outputs[0],
    #    bpy.data.scenes["Scene"].node_tree.nodes["Composite"].inputs[0])

    # Set the filename for the render output PNG
    bpy.context.scene.render.filepath = renderpath + final_render_filename + '_' + str(
        time).zfill(12) + '.png'

    # Turn on denoiser if using the CYCLES engine
    if bpy.context.scene.render.engine == 'CYCLES':
        bpy.context.scene.cycles.samples = 150
        bpy.context.scene.cycles.use_denoising = True
        bpy.context.scene.cycles.denoiser = 'NLM'

    # Preform rendering
    bpy.ops.render.render(use_viewport=True, write_still=True)
