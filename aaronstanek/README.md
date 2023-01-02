# ML Band Gaps (Materials)

> Ideal candidate: skilled ML data scientist with solid knowledge of materials science.

# Overview

The aim of this task is to create a python package that implements automatic prediction of electronic band gaps for a set of materials based on training data.

# User story

- As a user of this software, I can download sample training data from Materials Project, using my API key.
- As a user of this software, I can easily input the structural information of a target material.
- As a user of this software, I can predict the value of an electronic band gap with known accuracy after passing training data and structural information about the target material.

# Requirements

- suggest the bandgap values for a set of materials designated by their crystallographic and stoichiometric properties
- provide an intuitive, object-oriented, and documented API for downloading training data from the internet, reading training data from a file, collecting structural information about a target material from a Python interface or from a file, cleaning data, training, testing, and predicting.
- the code shall be written in a way that can facilitate easy addition of other characteristics extracted from simulations (forces, pressures, phonon frequencies etc)
- the package shall focus on usability and rapid adoption

# Expectations

- the code shall be able to suggest realistic values for slightly modified geometry sets - eg. trained on Si and Ge it should suggest the value of bandgap for Si49Ge51 to be between those of Si and Ge
- modular and object-oriented implementation
- commit early and often - at least once per 24 hours

# Timeline

We leave exact timing to the candidate. Must fit Within 5 days total.

# Notes

- use a designated github repository for version control
- suggested source of training data: materialsproject.org
