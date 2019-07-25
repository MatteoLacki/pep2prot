# From Peptidome to Proteome via Mainzer Dome

Our software aims at infering protein levels from observed peptide levels in a mass spectrometry experiment.
The procedures here are fairly general, but for now we effectively can turn an IsoQuant peptide report into
a protein report.

To install our software, you will need python3.
To get it, follow instructions here https://www.python.org/downloads/

Then, if you use windows, add python and pip to the system path (google it out).

Then, open the terminal and run:
```{bash}
pip install git+https://github.com/MatteoLacki/pep2prot.git
```

or even simpler:
```{bash}
pip install pep2prot
```

On Linux, you're done: check out
```{bash}
pep2prot -h
```

On Windows, you might want to additionally use the sendTo explorer option.
To do this, run 

```{bash}
console:sendto
```
in the windows menu command line (not powershell, not cmd, the thing where you search for software after pressing the windows button).
This will open up the sendto folder.
To this folder, add a shortcut of the file **_pep2prot**, which can be found in the **bin** folder of our python module (you have to figure out where it was installed. You can do it, for instance, by openining python prompt and calling 

```{python}
import pep2prot
pep2prot.__file__
```
although you might want to google that one out too).
With the shortcut in the sendto folder, you can now use **pep2prot** without ever needing to touch the terminal again.
To do this, simply open the explorer, go to the folder where both the peptide report and the fasta file are, highlight them, right click, send to and choose **_pep2prot**.


Best Regards,
Matteo Lacki
