from setuptools import setup

setup(
    name="pos_to_vtk",
    version="0.1.0",
    py_modules=["pos"],
    install_requires=[
            "Click",
            "pyvista",
            "numpy",
            "scipy",
            "tqdm",
    ],
    entry_points={
        'console_scripts': [
            'pos_to_vtk = pos_to_vtk:cli',
        ],
    },
)
