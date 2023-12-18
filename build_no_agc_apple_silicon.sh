#pushd WFA2-lib
#make all
#popd

#rustup default stable
rustup default stable-aarch64-apple-darwin
cargo build -p pgr-db --release --no-default-features
cargo build -p pgr-bin --release --no-default-features
cargo install --path pgr-bin --no-default-features

pushd pgr-tk/
maturin build --release --no-default-features
maturin build --release --skip-auditwheel --no-default-features
popd
