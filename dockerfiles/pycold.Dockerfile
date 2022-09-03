# syntax=docker/dockerfile:1.3.0-labs

# This dockerfile uses new-ish buildkit syntax. 
# Details on how to run are on the bottom of the file.
# (docker devs: todo unconsequential heredocs)

FROM pyenv:310

ARG HOME=/root

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
python -m pip install -r requirements/gdal.txt
python -m pip install setuptools==63.2.0  
EOF


# Stage the pycold source
COPY run_developer_setup.sh /pycold/
COPY setup.py /pycold/
COPY pyproject.toml /pycold/
COPY CMakeLists.txt /pycold/
COPY requirements.txt /pycold/
COPY requirements /pycold/requirements
COPY src /pycold/src
COPY tool /pycold/tool
COPY tests /pycold/tests
#COPY . /pycold

#ARG BUILD_STRICT=0

RUN <<EOF
#!/bin/bash
source $HOME/activate
set -x
cd /pycold
ls 
pwd
pip install -r requirements/build.txt
chmod +x run_developer_setup.sh
source run_developer_setup.sh
EOF


# Run simple tests
RUN <<EOF
#!/bin/bash
source $HOME/activate
cd /pycold
ls 
tree
echo "Start simple tests"
python -c "import pycold; print(pycold.__file__)"
EAGER_IMPORT=1 python -c "import pycold; print(pycold.__version__)"
EAGER_IMPORT=1 python -m pycold --help
EOF

# Copy over the rest of the repo
COPY . /pycold


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

docker run --runtime=nvidia -it pycold:310 bash
# docker buildx build -t "pyenv3.10" -f ./pyenv.Dockerfile --build-arg BUILD_STRICT=1 .
"


EOF

