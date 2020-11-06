## FCS Manager
I first developed FCS Manager for my research position at McGill University in 2018 to convert FCS files into PANDAS DataFrames for easier exploratory analysis than was available using FlowJo at the time. 

I used FCS Manager as a proof-of-concept to test automation capabilities for a high-throughput antibody screening data pipeline. Since then, I've used this FCS Manager to produce analyses included in at least one high-level grant application. 

**Only tested and safe for use with FCS3.0 and FCS3.1 files.** If you intend to use FCS Manager in your research, I **strongly urge** you to cross-check your analyses with FlowJo or similar flow package before integrating this code into your project. </br>**Tested on FCS data from BD FACSAria, BD FACSCalibur, and guava easyCyte machines.**

## Usage
Where ```fcs``` is the path to an FCS file e.g. ```'./data/experiment.fcs'``` and ```sample_number``` is the number of samples represented in that FCS file, usually ```0``` for single-experiment files, though some machines record multiple experiment data in a wrapper FCS file, with a subfolder containing individual FCS files - in which case ```sample_number``` would be ```30``` if 30 experiments (individual FCS files) contained in the wrapper FCS file.

```
import fcs_manager
fcs_manager.convertDF(fcs, sample_number)
```



Feel free to copy, modify, use, and distribute as you see fit.

## Cite
If you find this work useful in your own research, please cite as follows:

```
@misc{Zenkov-FCS-Manager,
  author = {Zenkov, Ilia},
  title = {fcs-manager},
  year = {2020},
  publisher = {GitHub},
  journal = {GitHub repository},
  howpublished = {\url{https://github.com/IliaZenkov/fcs-manager}},
}
```
## Licence

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/IliaZenkov/fcs-manager/blob/master/LICENSE)
