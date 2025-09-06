
# Table of Contents

1.  [Purpose](#org3af003a)
2.  [Dependency](#orgc582e17)
3.  [Compilation](#org802b803)
4.  [Usage](#org666a043)
    1.  [Input files](#org62cf43d)
    2.  [Output files](#orgd1dd586)



<a id="org3af003a"></a>

# Purpose

Calculate two-body entropy from pair correlation function.


<a id="orgc582e17"></a>

# Dependency

GSL - GNU Scientific Library is required. In general, it is installed
by default in Linux operating systems. Or you can refer to [GSL
website](https://www.gnu.org/software/gsl/).


<a id="org802b803"></a>

# Compilation

Type "make" and it will generate a "pdf2s2\_v2" binary.


<a id="org666a043"></a>

# Usage

"pdf2s2\_v2 -x XDATCAR" for calculation or "pdf2s2\_v2 -h" for help page.


<a id="org62cf43d"></a>

## Input files

1.  Mass
    This file contains three rows:
    -   type of species.
    -   atomic numbers.
    -   atomic masses [a.u.] in atomic units.

2.  Trun
    This file contains in value which is the temperature of MD simulation.

3.  XDATCAR
    Standard output of VASP runs.

4.  pdf files
    Pair distribution functions. You can generate those using "pdfxdat"
    command.


<a id="orgd1dd586"></a>

## Output files

1.  pdf.??.s2
    Two-body entropy as a function of integral distance of each kind of
    pair in a form of
    
    <table border="2" cellspacing="0" cellpadding="6" rules="groups" frame="hsides">
    
    
    <colgroup>
    <col  class="org-left" />
    
    <col  class="org-left" />
    
    <col  class="org-left" />
    
    <col  class="org-left" />
    </colgroup>
    <tbody>
    <tr>
    <td class="org-left">R [A]</td>
    <td class="org-left">S2(fluct) [kB]</td>
    <td class="org-left">S2(info) [kB]</td>
    <td class="org-left">S2(tot) [kB]</td>
    </tr>
    </tbody>
    </table>

2.  tot.s2
    Total two-body entropy as a function of integral distance  a form of
    
    <table border="2" cellspacing="0" cellpadding="6" rules="groups" frame="hsides">
    
    
    <colgroup>
    <col  class="org-left" />
    
    <col  class="org-left" />
    
    <col  class="org-left" />
    
    <col  class="org-left" />
    </colgroup>
    <tbody>
    <tr>
    <td class="org-left">R [A]</td>
    <td class="org-left">S2(fluct) [kB]</td>
    <td class="org-left">S2(info) [kB]</td>
    <td class="org-left">S2(tot) [kB]</td>
    </tr>
    </tbody>
    </table>

