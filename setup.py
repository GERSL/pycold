# from os.path import exists
# from os.path import join
# from os.path import dirname
"""

# Check contents of whel
python setup.py bdist_wheel
unzip -l dist/pycold-0.1.0-cp38-cp38-linux_x86_64.whl

pip install dist/pycold-0.1.0-cp38-cp38-linux_x86_64.whl
python -c "import pycold"
"""
from setuptools import find_packages
from skbuild import setup


def parse_version(fpath):
    """
    Statically parse the version number from a python file
    """
    import ast
    from os.path import exists
    if not exists(fpath):
        raise ValueError('fpath={!r} does not exist'.format(fpath))
    with open(fpath, 'r') as file_:
        sourcecode = file_.read()
    pt = ast.parse(sourcecode)
    class VersionVisitor(ast.NodeVisitor):
        def visit_Assign(self, node):
            for target in node.targets:
                if getattr(target, 'id', None) == '__version__':
                    self.version = node.value.s
    visitor = VersionVisitor()
    visitor.visit(pt)
    return visitor.version


VERSION = parse_version('tool/python/pycold/__init__.py')  # needs to be a global var for git tags

if __name__ == '__main__':
    setup(
        package_dir={'': 'tool/python/'},
        name="pycold",
        version=VERSION,
        description="python implementation of COntinuous monitoring of Land disturbances algorithm",
        install_requires=[
            # See requirements.txt for a more complete list
            'numpy >= 1.19.2',
        ],
        # ext_modules = cythonize([sccd_extension]),
        # include_dirs=[numpy.get_include()],
        author="Su Ye",
        author_email="remotesensingsuy@gmail.com",
        packages=find_packages(where='tool/python', include='pycold.*'),
        # cmake_install_dir='tool/python/pycold',
        # packages=['pycold'],
    )
