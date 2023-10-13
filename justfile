build_no_agc:
    #!/usr/bin/env bash
    #pushd WFA2-lib
    #make all
    #popd

    rustup default stable
    cargo build -p pgr-db --release --no-default-features
    cargo build -p pgr-bin --release --no-default-features
    cargo install --path pgr-bin --no-default-features

    pushd pgr-tk/
    maturin build --release --no-default-features
    maturin build --release --skip-auditwheel --no-default-features
    popd

build:
    #!/usr/bin/env bash
    rustup default stable

    ## if necessary, you can instal libclang / clang using Anaconda 
    ## and set LIBCLANG_PATH to point to the libclang for cbindgen dependence clang-sys
    # export LIBCLANG_PATH=$HOME/miniconda3/lib

    ## if necessary, install maturin with `cargo install --locked maturin`
    # cargo install --locked maturin

    cargo build -p pgr-db --release
    cargo build -p pgr-bin --release
    cargo install --path pgr-bin

    pushd pgr-tk/
    maturin build --release
    maturin build --release --skip-auditwheel
    popd
