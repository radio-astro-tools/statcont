STATCONT - A statistical continuum level determination method for line-rich sources
-----------------------------------------------------------------------------------

STATCONT is a python-based tool designed to determine the continuum 
emission level in line-rich spectral data. The tool inspects the 
intensity distribution of a given spectrum and automatically determines 
the continuum level by using differeng statistical approaches. The 
different methods included in STATCONT have been tested against 
synthetic data. We conclude that the sigma-clipping algorithm provides 
the most accurate continuum level determination, together with 
information on the uncertainty in its determination. This uncertainty 
is used to correct the final continuum emission level, resulting in the 
here-called 'corrected sigma-clipping method' or cSCM. The cSCM has 
been tested against synthetic data cubes reproducing typical conditions 
found in astronomical line-rich sources. In general, we obtain 
accuracies of < 10 % in the continuum determination, and < 5 % in most 
cases. The main products of STATCONT are the continuum emission level, 
together with its uncertainty, and data cubes containing only spectral 
line emission, i.e. continuum-subtracted data cubes. STATCONT also 
includes the option to estimate the spectral index or variation of the 
continuum emission with frequency.

If you find STATCONT useful, please cite/refer to:
[Sanchez-Monge, Schilke, Ginsburg, Cesaroni and Schmiedeke 2018, A&A, 609, A101](https://ui.adsabs.harvard.edu/abs/2018A%26A...609A.101S/abstract)
Access to the article in [PDF](https://hera.ph1.uni-koeln.de/~sanchez/my_documents/publications/SanchezMonge_etal_2018.pdf) format.

More information about STATCONT can be found in the webpage [http://www.astro.uni-koeln.de/~sanchez/statcont](http://www.astro.uni-koeln.de/~sanchez/statcont)

This README document is intented to help you install and use STATCONT using some test-cases provided below. In the following you will find:
  - [**Installation instructions**](https://github.com/radio-astro-tools/statcont#installation-instructions)
  - [**Required Python packages**](https://github.com/radio-astro-tools/statcont#required-python-packages)
  - [**Examples and test cases**](https://github.com/radio-astro-tools/statcont#examples-and-test-cases)
    - [**Determining the continuum in single spectrum files (ASCII files)**](https://github.com/radio-astro-tools/statcont/blob/master/README.md#determining-the-continuum-in-single-spectrum-files-ascii-files)
    - [**Determining the continuum in a 3D cube file (FITS files)**](https://github.com/radio-astro-tools/statcont/blob/master/README.md#determining-the-continuum-in-a-3d-cube-file-fits-files)
  - [**Publications making use of STATCONT**](https://github.com/radio-astro-tools/statcont#publications-making-use-of-statcont)

------------------------------------
Installation instructions
------------------------------------

STATCONT uses the [ASTROPY](https://www.astropy.org/) package-template and is
fully compatible with the ASTROPY ecosystem. It is freely available for at
the [GitHub](https://github.com/) repository [Radio Astro Tools](http://radio-astro-tools.github.io/),
and in this webpage. The only required software to use STATCONT is Python and Astropy.

You can directly install STATCONT by typing in a terminal session in 
your computer (you may need sudo permissions depending on the 
configuration of your system):

```
pip install https://github.com/radio-astro-tools/statcont/archive/master.zip
```

You can also clone STATCONT in your computer from the GitHub repository. For
this, create a directory and move there. Then type:

```
git init
git clone https://github.com/radio-astro-tools/statcont
cd statcont
python setup.py install
```

Alternatively, STATCONT can also be downloaded locally as a zip file 
from: [master.zip](https://github.com/radio-astro-tools/statcont/archive/master.zip).
In order to install it, download the file to a directory in your computer:

```
unzip master.zip
cd statcont
python setup.py install
```

If you get a permission error, this means that you do not have the required
administrative access to install STATCONT to your Python installation. In
this case you can use the option --user when installing the setup.py script.
You can follow these instructions:

```
python setup.py install --user --install-scripts="~/bin/"
export PATH=~/bin/:$PATH
```

Following the installation, you have immediate access to STATCONT in 
your computer by typing "statcont" in a terminal session. For example, 
inspect the help by doing:

```
statcont --help
```

---------------------
Required Python packages
---------------------

STATCONT uses the following Python packages:

- [argparse](https://docs.python.org/3/library/argparse.html) : version >= 1.1
- [astropy](https://www.astropy.org/) : version >= 4.2.1
- [datetime](https://docs.python.org/3/library/datetime.html)
- [math](https://docs.python.org/3/library/math.html)
- [matplotlib](https://matplotlib.org/) : version >= 3.5.0
- [numpy](https://numpy.org/) : version >= 1.21.4
- [os](https://docs.python.org/3/library/os.html)
- [scipy](https://scipy.org/) : version >= 1.7.3

The versions indicated above have been recently tested, older versions
worked fine in previous tests and are expected to continue working.

---------------------
Examples and test cases
---------------------

In the following we explain how to execute the main tasks of STATCONT. 
A set of test cases is provided in this [test_cases.tar.gz](https://hera.ph1.uni-koeln.de/~sanchez/software/STATCONT/test_cases.tar.gz) file. 
Download the file to your computer and follow these instructions:

```
gunzip test_cases.tar.gz
tar -xvf test_cases.tar
```

This creates a directory called `statcont-tests`. Inside, you will
find a directory called data that contains three other subdirectories
`MAP_TESTS`, `SPEC_TESTS`, and `SPINDEX`

STATCONT requires of a directory data where the files to be processed 
are stored. By executing STATCONT, another directory called products 
will be generated. The files to be processed, can be directly saved in 
data or in subdirectories. In the examples provided here, we have a set 
of single-spectrum files (in ASCII format) saved in `SPEC_TESTS` and
`SPINDEX`, and FITS cubes in `MAP_TESTS`.


Determining the continuum in single spectrum files (ASCII files)
----------------------------------------------------------------

```
   statcont -p SPEC_TESTS -s my_emission -n 1
```

  - The option -p indicates the subdirectory in data that contains
    the file to be analyzed
  - The option -s indicates the name of the ASCII file to be
    analyzed, without the extension [.dat]
  - The option -n indicates the rms noise level (in the units of the
    data) of the data to be analyzed. In this case, it is 1 K

If you want to determine the continuum level:

```
   statcont -p SPEC_TESTS -s my_emission -n 1 --continuum
```

  - The option --continuum makes use of the 'corrected sigma-clipping 
    algorithm' described in [Sanchez-Monge et al (2018)](https://ui.adsabs.harvard.edu/abs/2018A%26A...609A.101S/abstract), to determine
    the continuum level, the error in the continuum level, and to
    produce a file that contains only the line emission, i.e. a
    continuum-subtracted file

Using different methods to determine the continuum level. STATCONT 
contains a set of different statistical methods that can be used by the 
user at his/her convenience. You can select all them like this:

```
   statcont -p SPEC_TESTS -s my_emission -n 1 --call
```

Or you can select individual methods like:

```
   statcont -p SPEC_TESTS -s my_emission -n 1 --cmax
   statcont -p SPEC_TESTS -s my_emission -n 1 --cmean
   statcont -p SPEC_TESTS -s my_emission -n 1 --cmedian
   statcont -p SPEC_TESTS -s my_emission -n 1 --cpercent
   statcont -p SPEC_TESTS -s my_emission -n 1 --cGaussian
   statcont -p SPEC_TESTS -s my_emission -n 1 --cKDEmax
   statcont -p SPEC_TESTS -s my_emission -n 1 --csigmaclip
```

You can call several methods at once:

```
   statcont -p SPEC_TESTS -s my_emission -n 1 --cmax --cGaussian --csigmaclip
```

The different methods are explained in [Sanchez-Monge et al (2018)](https://ui.adsabs.harvard.edu/abs/2018A%26A...609A.101S/abstract)

If you want to remove the continuum from the original spectrum, in 
order to produce a line-only data file, you can use:

```
   statcont -p SPEC_TESTS - s my_emission -n 1 --csigmaclip --cfree
```

  - The option --cfree uses the 'corrected sigma-clipping algorithm'
    to determine the continuum, and removes it from the original
    datafile.

You can use other example files, like for example:

```
   statcont -p SPEC_TESTS -s my_absorption -n 1 --continuum
```

And you can select multiple files simultaneously, as long as they are 
saved in the same subdirectory, and they are considered to have the 
same rms noise level (option -n ):

```
   statcont -p SPEC_TESTS -s my_emission my_absorption my_broad-lines -n 1 --continuum
```

The products can be found in products/SPEC_TESTS
You can produce plots of the spectrum analyzed with the continuum 
levels by using the option --plots. As an example:

```
   statcont -p SPEC_TESTS -s my_emission -n 1 --continuum --plots
```

In this case the plot is saved in products/SPEC_TESTS/plots/my_emission_1_1.png
You can use all the continuum methods and plot them all together, like:

```
   statcont -p SPEC_TESTS -s my_emission -n 1 --call --plots
```

Have a look now at the plot products/SPEC_TESTS/plots/my_emission_1_1.png (see image below)

![STATCONT all methods](test-cases/my_emission_1_1.png)
 

Determining the continuum in a 3D cube file (FITS files)
--------------------------------------------------------

```
   statcont -p MAP_TESTS -i SYNTHETIC_cube -n 1
```

  - The option -p indicates the subdirectory in data that contains
    the file to be analyzed
  - The option -i indicates the name of the FITS file to be analyzed,
    without the extension [.fits]
  - The option -n indicates the rms noise level (in the units of the
    data) of the data to be analyzed. In this case, it is 1 K

If you want to determine the continuum level:

```
   statcont -p MAP_TESTS -i SYNTHETIC_cube -n 1 --continuum
```

This process analyzes each individual pixel, determining the continuum 
level, and then combines all the pixels to produce a continuum FITS 
image with the label "_continuum". Simultaneously, the --csigmaclip 
method provides information on the error in the determination of the 
continuum that is saved as a FITS image with the label "_noise". 
Finally, a line-only FITS datacube is also produced. All these files 
are saved in the directory products/MAP_TESTS

All the other options applicable to the single-spectrum ASCII files are 
also available for the FITS images (e.g. different continuum methods, 
creation of plots).

If your original FITS file is too large and you just want to determine 
the continuum level of a small portion you can indicate it like this:

```
   statcont -p MAP_TESTS -i SYNTHETIC_cube -n 1 --continuum --cutout 25 25 6
```

  - The --cutout option allows to select a central pixel (in this case
    25, 25) and the number of pixels in each direction of the final
    image (in this case 6). With this option, the products are saved
    with the label "_cutout"

If you have multiple FITS files at different frequencies, you can use 
the option --spindex to determine, first, the continuum level of every 
single file, and then the spectral index, i.e. the variation of the 
continuum emission with frequency.


-----------------------------------
Publications making use of STATCONT
-----------------------------------

This is a list of publications using STATCONT:

 - *Radio outburst from a massive (proto)star. When accretion turns into ejection*<br /> by Cesaroni et al. 2018, A&A, 612, A103 ([link](https://ui.adsabs.harvard.edu/abs/2018A%26A...612A.103C/abstract))

 - *Distributed Star Formation throughout the Galactic Center Cloud Sgr B2*<br /> by Ginsburg et al. 2018, ApJ, 853, 171 ([link](https://ui.adsabs.harvard.edu/abs/2018ApJ...853..171G/abstract))

 - *The physical and chemical structure of Sagittarius B2 - II. Continuum millimeter emission of SgrB2(M) and SgrB2(N) with ALMA*<br /> by Sanchez-Monge et al. 2017, A&A, 604, A6 ([link](https://ui.adsabs.harvard.edu/abs/2017A%26A...604A...6S/abstract))

 - *Chasing disks around O-type (proto)stars: Evidence from ALMA observations*<br />  by Cesaroni et al. 2017, A&A, 602, A59 ([link](https://ui.adsabs.harvard.edu/abs/2017A%26A...602A..59C/abstract))
