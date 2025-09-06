                              ____________

                                 README

                               Yang Huang
                              ____________


                            <2023-04-27 Thu>


Table of Contents
_________________

1. Purpose
2. Dependency
3. Compilation
4. Usage
.. 1. Input files
.. 2. Output files





1 Purpose
=========

  Calculate two-body entropy from pair correlation function.


2 Dependency
============

  GSL - GNU Scientific Library is required. In general, it is installed
  by default in Linux operating systems. Or you can refer to [GSL
  website].


[GSL website] <https://www.gnu.org/software/gsl/>


3 Compilation
=============

  Type "make" and it will generate a "pdf2s2_v2" binary.


4 Usage
=======

  "pdf2s2_v2 -x XDATCAR" for calculation or "pdf2s2_v2 -h" for help
  page.


4.1 Input files
~~~~~~~~~~~~~~~

  1) Mass This file contains three rows:
     - type of species.
     - atomic numbers.
     - atomic masses [a.u.] in atomic units.

  2) Trun This file contains in value which is the temperature of MD
     simulation.

  3) XDATCAR Standard output of VASP runs.

  4) pdf files Pair distribution functions. You can generate those using
     "pdfxdat" command.


4.2 Output files
~~~~~~~~~~~~~~~~

  1) pdf.??.s2 Two-body entropy as a function of integral distance of
     each kind of pair in a form of
      R [A]  S2(fluct) [kB]  S2(info) [kB]  S2(tot) [kB]

  2) tot.s2 Total two-body entropy as a function of integral distance a
     form of
      R [A]  S2(fluct) [kB]  S2(info) [kB]  S2(tot) [kB]
