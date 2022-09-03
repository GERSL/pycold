# syntax=docker/dockerfile:1.3.0-labs

# This dockerfile uses new-ish buildkit syntax. 
# Details on how to run are on the bottom of the file.
# (docker devs: todo unconsequential heredocs)

FROM pyenv:310

ENV HOME=/root

RUN <<EOF
#!/bin/bash
apt update -q
DEBIAN_FRONTEND=noninteractive apt install -q -y --no-install-recommends \
        build-essential zlib1g-dev libgsl-dev gfortran \
        ffmpeg tmux jq tree
apt-get clean 
rm -rf /var/lib/apt/lists/*
EOF

SHELL ["/bin/bash", "--login", "-c"]

RUN echo $(pwd)

COPY requirements.txt /pycold/
COPY requirements /pycold/requirements

# Setup primary dependencies
RUN <<EOF
#!/bin/bash
source $HOME/activate
cd /pycold
# Always use the latest Python build tools
python -m pip install pip -U 
python -m pip install -r requirements.txt
python -m pip install -r requirements/build.txt
python -m pip install -r requirements/headless.txt
python -m pip install -r requirements/runtime.txt
python -m pip install -r requirements/optional.txt
python -m pip install -r requirements/tests.txt
python -m pip install -r requirements/gdal.txt
python -m pip install setuptools==63.2.0  
EOF


# Copy the pycold source (being careful to not copy host build artifacts)

COPY .git /pycold/.git

#### Repo root files
#COPY *.sh /pycold/
#COPY *.py /pycold/
#COPY *.toml /pycold/
#COPY *.txt /pycold/

#### Subdirectory: resources
#COPY requirements/*.txt /pycold/requirements/

#### Subdirectory: tool/notebook
#COPY tool/notebook/*.ipynb /pycold/tool/notebook/

#### Subdirectory: tool/matlab
#COPY tool/matlab/*.m /pycold/tool/matlab/

#### Subdirectory: tests
#COPY tests/*.py /pycold/tests/

#### Subdirectory: tests/resources/feature_maps
#COPY tests/resources/*.npy /pycold/tests/resources
#COPY tests/resources/*.csv /pycold/tests/resources
#COPY tests/resources/*.yaml /pycold/tests/resources
#COPY tests/resources/*.pickle /pycold/tests/resources
#COPY tests/resources/*.pickle /pycold/tests/resources

#### Subdirectory: tests/resources/feature_maps
#COPY tests/resources/feature_maps/*.npy   /pycold/tests/resources/feature_maps/
#COPY tests/resources/feature_maps/*.model /pycold/tests/resources/feature_maps/
#COPY tests/resources/feature_maps/*.npy   /pycold/tests/resources/feature_maps/

#### Subdirectory: tests/resources/cm_maps
#COPY tests/resources/cm_maps/*.npy        /pycold/tests/resources/cm_maps/

#### Subdirectory: src/cxx
#COPY src/cxx/*.c /pycold/src/cxx/
#COPY src/cxx/*.h /pycold/src/cxx/
#COPY src/cxx/*.txt /pycold/src/cxx/

#### Subdirectory: src/python/pycold
#COPY src/python/pycold/*.py /pycold/src/python/pycold/
#COPY src/python/pycold/*.yaml /pycold/src/python/pycold/
#COPY src/python/pycold/*.txt /pycold/src/python/pycold/
#COPY src/python/pycold/*.pyx /pycold/src/python/pycold/

#### Subdirectory: src/python/pycold/imagetool
#COPY src/python/pycold/imagetool/*.py /pycold/src/python/pycold/imagetool/
#COPY src/python/pycold/imagetool/*.yaml /pycold/src/python/pycold/imagetool/
#COPY src/python/pycold/imagetool/*.tif /pycold/src/python/pycold/imagetool/


WORKDIR /pycold

RUN <<EOF
#!/bin/bash
ls -al
tree
git reset --hard HEAD
ls -al
EOF


RUN <<EOF
#!/bin/bash
source $HOME/activate
pip install -r requirements/build.txt
python -m pip install setuptools==63.2.0  
chmod +x run_developer_setup.sh
./run_developer_setup.sh
EOF


## Run simple tests
RUN <<EOF
#!/bin/bash
source $HOME/activate
tree
echo "Start simple tests"
python -c "import pycold; print(pycold.__file__)"
python -c "import pycold; print(pycold.__version__)"
EOF


################
### __DOCS__ ###
################
RUN <<EOF


# https://www.docker.com/blog/introduction-to-heredocs-in-dockerfiles/
echo "
# docker login
# docker pull docker/dockerfile:1.3.0-labs
cd $HOME/code/pycold

DOCKER_BUILDKIT=1 docker build --progress=plain \
    -t "pyenv:310" \
    --build-arg PYTHON_VERSION=3.10.5 \
    -f ./dockerfiles/pyenv.Dockerfile .

DOCKER_BUILDKIT=1 docker build --progress=plain \
    -t "pycold:310" \
    -f ./dockerfiles/pycold.Dockerfile .

docker run -w /pycold -it pycold:310 bash
# docker buildx build -t "pyenv3.10" -f ./pyenv.Dockerfile --build-arg BUILD_STRICT=1 .
"


EOF

