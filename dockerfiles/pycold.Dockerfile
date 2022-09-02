# syntax=docker/dockerfile:1.3.0-labs

# This dockerfile uses new-ish buildkit syntax. 
# Details on how to run are on the bottom of the file.
# (docker devs: todo unconsequential heredocs)

FROM nvidia/cuda:11.4.3-cudnn8-devel-ubuntu20.04

ARG HOME=/root
ARG PYENV_ROOT=/root/.pyenv
ARG CHOSEN_PYTHON_VERSION=3.10.5

## Install Prerequisites 
RUN <<EOF
#!/bin/bash
apt update -q
DEBIAN_FRONTEND=noninteractive apt install -q -y --no-install-recommends \
        bzip2 \
        ca-certificates \
        git \
        libgl1-mesa-glx \
        libglib2.0-0 \
        libsm6 \
        libxext6 \
        libxrender1 \
        mercurial \
        subversion \
        wget \
        make build-essential libssl-dev zlib1g-dev \
        libbz2-dev libreadline-dev libsqlite3-dev wget curl llvm libncurses5-dev \
        libncursesw5-dev xz-utils tk-dev libffi-dev liblzma-dev python-openssl \
        ffmpeg tmux jq
apt-get clean 
rm -rf /var/lib/apt/lists/*
EOF

## Install pyenv
RUN <<EOF
#!/bin/bash
git clone https://github.com/pyenv/pyenv.git -b v2.3.3 $PYENV_ROOT 
(cd $PYENV_ROOT && src/configure && make -C src)
EOF


## Use pyenv to compile an optimized Python
RUN <<EOF
#!/bin/bash
export PATH="$PYENV_ROOT/bin:$PATH"
eval "$($PYENV_ROOT/bin/pyenv init -)"

#PROFILE_TASK="-m test.regrtest --pgo test_array test_base64 test_binascii test_binhex test_binop test_c_locale_coercion test_csv test_json test_hashlib test_unicode test_codecs test_traceback test_decimal test_math test_compile test_threading test_time test_fstring test_re test_float test_class test_cmath test_complex test_iter test_struct test_slice test_set test_dict test_long test_bytes test_memoryview test_io test_pickle"

#PYTHON_CONFIGURE_OPTS="--enable-shared --enable-optimizations --with-computed-gotos --with-lto"

#PYTHON_CFLAGS="-march=native -O2 -pipe" 

#PROFILE_TASK=$PROFILE_TASK \
#PYTHON_CFLAGS="$PYTHON_CFLAGS" \
#PYTHON_CONFIGURE_OPTS="$PYTHON_CONFIGURE_OPTS" \
pyenv install $CHOSEN_PYTHON_VERSION 
EOF


## Create a default Python virtualenv
RUN <<EOF
#!/bin/bash
export PATH="$PYENV_ROOT/bin:$PATH"
eval "$($PYENV_ROOT/bin/pyenv init -)"
pyenv global $CHOSEN_PYTHON_VERSION

PYENV_PREFIX=$(pyenv prefix)
python -m venv $PYENV_PREFIX/envs/pyenv$CHOSEN_PYTHON_VERSION

BASHRC_CONTENTS='
# Add the pyenv command to our environment if it exists
export HOME="/root"
export PYENV_ROOT="$HOME/.pyenv"
if [ -d "$PYENV_ROOT" ]; then
    export PATH="$PYENV_ROOT/bin:$PATH"
    eval "$($PYENV_ROOT/bin/pyenv init -)"
    source $PYENV_ROOT/completions/pyenv.bash
    export PYENV_PREFIX=$(pyenv prefix)
fi

# Optionally auto-activate the chosen pyenv pyenv environment
CHOSEN_PYTHON_VERSION=3.10.5
if [ -d "$PYENV_PREFIX/envs/pyenv$CHOSEN_PYTHON_VERSION" ]; then
    source $PYENV_PREFIX/envs/pyenv$CHOSEN_PYTHON_VERSION/bin/activate
fi
'
echo "$BASHRC_CONTENTS" >> $HOME/.bashrc
# Write a secondary script for non-interactive usage
echo "$BASHRC_CONTENTS" >> $HOME/activate
EOF

#### UNCOMMENT FOR DEBUGGING
## precache big pip packages if we are debugging steps after this
#RUN <<EOF
##!/bin/bash
#source $HOME/activate
#pip install pip -U
#pip install torch==1.11.0
#EOF


# Stage the pycold source
COPY setup.py /pycold/
COPY pyproject.toml /pycold/
COPY CMakeLists.txt /pycold/CMakeLists.txt
COPY requirements.txt /pycold/requirements.txt
COPY run_developer_setup.sh /pycold/run_developer_setup.sh
COPY requirements /pycold/requirements
COPY src /pycold/src
COPY tests /pycold/tests
COPY tool /pycold/tool
COPY . /pycold

SHELL ["/bin/bash", "--login", "-c"]

RUN echo $(pwd)

ARG BUILD_STRICT=0

# Setup primary dependencies
RUN <<EOF
#!/bin/bash
source $HOME/activate

# Always use the latest Python build tools
python -m pip install pip setuptools wheel build -U

cd /pycold
bash run_developer_setup.py

#if [ "$BUILD_STRICT" -eq 1 ]; then
#    echo "BUILDING STRICT VARIANT"
#    #pip install -e /pycold[runtime-strict,optional-strict,headless-strict]
#    pip install -e /pycold[runtime-strict,optional-strict]
#else
#    echo "BUILDING LOOSE VARIANT"
#    pip install -e /pycold[optional]
#    # python -m pip install dvc[all]>=2.13.0
#    # pip install awscli
#fi

EOF


# Finalize more fickle dependencies
RUN <<EOF
#!/bin/bash
source $HOME/activate

cd /pycold
if [ "$BUILD_STRICT" -eq 1 ]; then
    echo "FINALIZE STRICT VARIANT DEPS"
    sed 's/>=/==/g' requirements/gdal.txt > requirements/gdal-strict.txt
    pip install -r requirements/gdal-strict.txt
else
    echo "FINALIZE LOOSE VARIANT DEPS"
    pip install -r requirements/gdal.txt
fi

EOF


# Install other useful tools
RUN <<EOF
#!/bin/bash
source $HOME/activate

# python -m pip install dvc[all]>=2.13.0
# pip install scikit-image>=0.18.1
pip install awscli

EOF


# Run simple tests
RUN <<EOF

#!/bin/bash
source $HOME/activate

echo "Start simple tests"
EAGER_IMPORT=1 python -c "import pycold; print(pycold.__version__)"
EAGER_IMPORT=1 python -m pycold --help
EOF

## Copy over the rest of the repo
#COPY . /pycold


################
### __DOCS__ ###
################
RUN <<EOF


# https://www.docker.com/blog/introduction-to-heredocs-in-dockerfiles/
echo "
# docker login
# docker pull docker/dockerfile:1.3.0-labs
cd $HOME/code/pycold
DOCKER_BUILDKIT=1 docker build --progress=plain -t "pycold:310" -f ./dockerfiles/pycold.Dockerfile  .

docker run --runtime=nvidia -it pycold:310 bash
# docker buildx build -t "pyenv3.10" -f ./pyenv.Dockerfile --build-arg BUILD_STRICT=1 .
"


EOF

