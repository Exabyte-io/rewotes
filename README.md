# mlbands

A python package that implements automatic prediction of electronic band gaps for a set of materials based on training data.


## Installation

```ruby
pip install --upgrade mlbands
```

## Documentation and Usage on Google Colab (click below)

<a href="https://colab.research.google.com/drive/14GS_jUo_B6ojDip-Ak2VsGdrZL1Fgx5a?usp=sharing">
<img src="https://github.com/andrewrgarcia/powerxrd/blob/main/img/colab.png?raw=true" width="300" ></a>



# ML Band Gaps (Materials)

> Ideal candidate: skilled ML data scientist with solid knowledge of materials science.

# Overview

The aim of this task is to create a python package that implements automatic prediction of electronic band gaps for a set of materials based on training data.

# User story

As a user of this software I can predict the value of an electronic band gap after passing training data and structural information about the target material.

# Requirements

- suggest the bandgap values for a set of materials designated by their crystallographic and stoichiometric properties
- the code shall be written in a way that can facilitate easy addition of other characteristics extracted from simulations (forces, pressures, phonon frequencies etc)

# Expectations

- the code shall be able to suggest realistic values for slightly modified geometry sets - eg. trained on Si and Ge it should suggest the value of bandgap for Si49Ge51 to be between those of Si and Ge
- modular and object-oriented implementation
- commit early and often - at least once per 24 hours

# Timeline

We leave exact timing to the candidate. Must fit Within 5 days total.

# Notes

- use a designated github repository for version control
- suggested source of training data: materialsproject.org
