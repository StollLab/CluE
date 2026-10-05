## CluE
CluE (Cluester Evolution) is a quantum spin dynamics program to simulate central spin decoherence,  CluE implements Yang and Liu's cluster correlation expansion [1,2].  For usage information please refer to the manual.

### Installation
Before installing, make sure you have the latest version of
[Rust](https://www.rust-lang.org/)[3]
installed.

Before compiling CluE, it is a good idea to check that everything is working properly by running the tests. 
```
cargo test
```
To build Clue, run the following
```
cargo build --release
```
This will build the binary in `clue_oxide/target/release`; to make clue globally available either add this directory to your system path or establish as alias.
In Bash, add
```
alias clue="path/to/clue/target/release/clue_oxide"
```
to your `.bash_aliases` file.

## pyCluE
Within the CluE source directory, navigate to the pyclue directory.
```
cd pyclue
```

The Python interface uses [maturin](https://github.com/PyO3/maturin)[4] to compile.
To install maturin in a Python virtual environment use the following.
```
python3 -m venv <path/to/virtual/environment>
source <path/to/virtual/environment/bin/activate >
pip install maturin
```
And the use the following to build the Python interface.
```
maturin build
``` 
One potential issue is that when maturin tries to compile CluE, it can fail to see the operating system unique
flags, and will try to use them all. To account for this, open CluE’s Cargo.toml file and comment out
everything under the unneeded operating system section. For example, to compile on Linux add comments
as shown below.
```
[target.'cfg(unix)'.dependencies]
ndarray -linalg = { version = "0.15", features = ["openblas -static"] }

#[target.'cfg(windows) '.dependencies.ndarray -linalg]
#version = '0.15.0 '
#features = ['intel -mkl ']
```
Once built, navigate to `target/wheels`, source the desired Python environment, and use pip to install the
wheel.
```
cd target/wheels
source <path/to/installation/environment/bin/activate >
pip install <pyclue.whl>
```

## References
<a id="1">[1]</a> 
Yang, W.; Liu, R. B. Decoherence of Coupled Electron Spins via Nuclear Spin Dynamics in Quantum Dots. Phys. Rev. B 2008, 77 (8), 085302. https://doi.org/10.1103/PhysRevB.77.085302.

<a id="2">[2]</a> 
Yang, W.; Liu, R.-B. Quantum Many-Body Theory of Qubit Decoherence in a Finite-Size Spin Bath. II. Ensemble Dynamics. Phys. Rev. B 2009, 79 (11), 115320. https://doi.org/10.1103/PhysRevB.79.115320.

<a id="3">[3]</a>
Foundation, R. Rust A language empowering everyone to build reliable and efficient software. https://www.rust-lang.org/.

<a id="3">[4]</a>
PyO3 PyO3/maturin https://github.com/PyO3/maturin.

