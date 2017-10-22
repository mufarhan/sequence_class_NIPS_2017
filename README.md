# Efficient-Approximation-Algorithms-for-Strings-Kernel-Based-Sequence-Classification

References:
-----------
Farhan, M., Tariq, J., Zaman, A., Shabbir, M., Khan, I.  "Efficient Approximation Algorithms for Strings Kernel Based Sequence Classification." NIPS, 2017

This implements new approximation based kernel algorithms for strings (see references).

Data Sources:
-------------

Ding-Dubchak - http://ranger.uta.edu/~chqding/protein/ <br />
Music Genre - http://opihi.cs.uvic.ca/sound/genres/ <br />
ISMIR Contest - http://ismir2004.ismir.net/genre_contest/ <br />
Artist20 - https://labrosa.ee.columbia.edu/projects/artistid/ <br />

Installation:
-------------
Download Source Code folder to your desired directory, and make sure, environment variable for java is set
> cd Source Code <br />
> javac Main.java 

This will produce executable class files for string kernel

computations:
-------------
Main.class : computes mismatch kernel matrices

Usage:
------
This function takes a text file with sequences and output a text file with a kernel matrix.

E.g. to compute mismatch(8,2) kernel for sequences with alphabet size=1024: <br />
> java Main music.genre.txt 1000 8 2 1024 <br />

this will create Kernel-k8-m2.txt file with kernel matrix <br />

String kernels are called with the following parameters: <br />
java Main \<Sequence-file\> \<# of Sequences\> \<k\> \<m\> \<AlphabetSize\> <br />
<br />
where <br />
\<Sequence-file\> is the file with sequence data: <br />
one sequence per line (line should end with line feed), with sequence elements separated by space, <br />
all sequence elements are assumed to be in the range [0, \<AlphabetSize\> - 1]. <br />
See Datasets folder for an example of the sequence file format. <br />
\<k\>,\<m\>,\<b\>,\<sigma\> are corresponding kernel parameters. (see references and help for a particular function for details) <br />
\<# of Sequences\> is the total number of sequences. <br />
\<AlphabetSize\> is the size of the alphabet. <br />

Output kernel matrix is written into <br />
  Kernel-k\<k\>-m\<m\>.txt <br />
file. <br />

Authors:
--------
Muhammad Farhan
14030031@lums.edu.pk <br />
Imdad Ullah Khan
imdad.khan@lums.edu.pk
