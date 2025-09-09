Liquid-entropy-analysis-tools (LEAT) is package for analyzing entropy in liquid systems, including tools for  calculating pair distribution function (PDF), 3-body distribution function (3BDF), two-body entropy (S2), three-body entropy (S3), etc. It requires input and output files in VASP: INCAR, OUTCAR, POTCAR and XDATCAR.

It requires GNU Scientific Library (GSL), which, in general, is installed by default in Linux operating systems. Or you can refer to [GSL website](https://www.gnu.org/software/gsl/).

To install, 

    git clone https://github.com/yangh618/Liquid-entropy-analysis-tools.git
    cd Liquid-entropy-analysis-tools
    ./install

After compilation, all utilities will be available in the bin/ directory. To use them conveniently, either add bin/ to your PATH, or  create symbolic links to the binaries from a directory already in your PATH.

You can try out the examples located in the examples/ directory. For detailed usage of each command, please refer to the README file in the corresponding directory in utils/.

If you are using 3BDF tools, please cite

> Huang, Y. and Widom, M., 2024. Entropy approximations for simple fluids. Physical Review E, 109(3), p.034130.

If you are using S2 tools, please cite

> Huang, Y., Widom, M. and Gao, M.C., 2022. Ab initio free energies of liquid metal alloys: Application to the phase diagrams of Li-Na and K-Na. Physical Review Materials, 6(1), p.013802.

