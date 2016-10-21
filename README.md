#IonNIPT
==============

#### ThermoFisher Ion Torrent plugin to detect fetal trisomies and estimate fetal fraction

## 1 Tools and dependencies
Tools used from other projects are shown below.
![Screenshot](https://raw.githubusercontent.com/AllanSSX/IonNIPT/master/IonNIPT.PNG)

- Python   
--pickle 
--joblib
--scipy
--sklearn   

## 2 Installation

Simply clone the repository:

`git clone --recursive https://github.com/AllanSSX/IonNIPT.git`

You need to change some variables in IonNIPT.py:

- l.33: path to the trained model (Sanefalcon)
- l.36: the scaling factor value for Defrag
- l.37: the percYonMales value also for Defrag

Moreover, you need to traine Sanefalcon and Defrag on your data and push the results into data/ (Sanefalcon: nucleosome track, train model; Wisecondor: gcccount, reftable, and gcc/pickle for males and females).
See Sanefalcon and Wisecondor manual for further details.