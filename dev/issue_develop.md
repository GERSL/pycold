There is an issue installing a package in development mode when the base python
module folder is not in the root of the repo. I.e. the case where there is a
source directory:


I will use the module "pycold" and its submodule "pycold.demo" as the examples
for this issue.

```
$REPO/setup.py
$REPO/src/python/pycold/__init__.py
$REPO/src/python/pycold/demo/__init__.py
```


It seems to happen correctly when one "package", but fails with multiple packages.

The 'egg-base' seems to be set incorrectly.


I added debug print statements to setuptools/command/develop.py, setuptools/command/egg_info.py and skbuild/command/egg_info.py

Given:

```
    package_dir={
        '': 'src/python/',
        'pycold': 'src/python/pycold',
        'pycold.demo': 'src/python/pycold/demo',
        'pycold.demo.resources': 'src/python/pycold/demo/resources',
    },
```


When I specify 

```
packages=['pycold'],
```

I get:

```
    running develop
    !!! <FINALIZE EGG INFO> !!!
    self.egg_base = None
    !!! <--FINALIZE EGG INFO> !!!
    <setuptools-egg-info>
    self.distribution = <skbuild.setuptools_wrap.setup.<locals>.BinaryDistribution object at 0x7fe5a08860a0>
    self.distribution.package_dir = {'pycold.demo.resources': 'src/python/pycold/demo/resources', 'pycold.demo': 'src/python/pycold/demo', 'pycold': '_skbuild/linux-x86_64-3.8/cmake-install/src/python/pycold'}
    self.egg_info = 'pycold.egg-info'
    self.egg_base = '.'
    </setuptools-egg-info>
    self.egg_base = '.'
    !!! <//FINALIZE EGG INFO> !!!
    self.egg_base = '.'
    self.install_dir = '/home/joncrall/.pyenv/versions/3.8.6/envs/pyenv3.8.6/lib/python3.8/site-packages'
    self.egg_path = '/home/joncrall/code/pycold'
    self.setup_path = '.'
    running egg_info
    writing pycold.egg-info/PKG-INFO
```


But in the case where

```
packages=['pycold', 'pycold.demo', 'pycold.demo.resources'],
```

I get:

    running develop
    !!! <FINALIZE EGG INFO> !!!
    self.egg_base = 'src/python'
    self.egg_base = 'src/python'
    !!! <--FINALIZE EGG INFO> !!!
    <setuptools-egg-info>
    self.distribution = <skbuild.setuptools_wrap.setup.<locals>.BinaryDistribution object at 0x7f779a094a00>
    self.distribution.package_dir = {'pycold': '_skbuild/linux-x86_64-3.8/cmake-install/src/python/pycold'}
    self.egg_info = 'src/python/pycold.egg-info'
    self.egg_base = 'src/python'
    </setuptools-egg-info>
    self.egg_base = 'src/python'
    !!! <//FINALIZE EGG INFO> !!!
    self.egg_base = 'src/python'
    self.install_dir = '/home/joncrall/.pyenv/versions/3.8.6/envs/pyenv3.8.6/lib/python3.8/site-packages'
    self.egg_path = '/home/joncrall/code/pycold/src/python'
    self.setup_path = '../../'
    running egg_info
    writing src/python/pycold.egg-info/PKG-INFO



This is a problem because without specifying these packages, pycold.demo does
not get packaged in the wheel (although it is available in development mode)



If I remove details in `package_dir`:

```
    package_dir={
        '': 'src/python/',
    },
    packages=['pycold'],
```

I get:

```
    running develop
    !!! <FINALIZE EGG INFO> !!!
    self.egg_base = 'src/python'
    self.egg_base = 'src/python'
    !!! <--FINALIZE EGG INFO> !!!
    <setuptools-egg-info>
    self.distribution = <skbuild.setuptools_wrap.setup.<locals>.BinaryDistribution object at 0x7fd7c70c3eb0>
    self.distribution.package_dir = {'pycold': '_skbuild/linux-x86_64-3.8/cmake-install/src/python/pycold'}
    self.egg_info = 'src/python/pycold.egg-info'
    self.egg_base = 'src/python'
    </setuptools-egg-info>
    self.egg_base = 'src/python'
    !!! <//FINALIZE EGG INFO> !!!
    self.egg_base = 'src/python'
    self.install_dir = '/home/joncrall/.pyenv/versions/3.8.6/envs/pyenv3.8.6/lib/python3.8/site-packages'
    self.egg_path = '/home/joncrall/code/pycold/src/python'
    self.setup_path = '../../'
    running egg_info
    writing src/python/pycold.egg-info/PKG-INFO
```

And then 

```
    package_dir={
        '': 'src/python/',
    },
    packages=['pycold', 'pycold.demo', 'pycold.demo.resources'],

```


```
    running develop
    !!! <FINALIZE EGG INFO> !!!
    self.egg_base = None
    !!! <--FINALIZE EGG INFO> !!!
    <setuptools-egg-info>
    self.distribution = <skbuild.setuptools_wrap.setup.<locals>.BinaryDistribution object at 0x7fdd61d7f820>
    self.distribution.package_dir = {'pycold.demo.resources': 'src/python/pycold/demo/resources', 'pycold.demo': 'src/python/pycold/demo', 'pycold': '_skbuild/linux-x86_64-3.8/cmake-install/src/python/pycold'}
    self.egg_info = 'pycold.egg-info'
    self.egg_base = '.'
    </setuptools-egg-info>
    self.egg_base = '.'
    !!! <//FINALIZE EGG INFO> !!!
    self.egg_base = '.'
    self.install_dir = '/home/joncrall/.pyenv/versions/3.8.6/envs/pyenv3.8.6/lib/python3.8/site-packages'
    self.egg_path = '/home/joncrall/code/pycold'
    self.setup_path = '.'
    running egg_info
    writing pycold.egg-info/PKG-INFO
```
