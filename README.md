# Bub-PTV: Bubble based Particle Tracking Velocimetry
* This demo code is for [2022 ANS THD Best Paper Award](https://thd.ans.org/awards/bestpaper/) & Best Paper Award in NUTHOS-13 for conference paper [Measurement of velocity induced by steam condensation into a water pool by tracking the motion of bubbles](https://www.researchgate.net/publication/363405940_Measurement_of_velocity_induced_by_steam_condensation_into_a_water_pool_by_tracking_the_motion_of_bubbles)
* It also refers to two journal publications [**under review**]

## Introduction
We introduce an experimental approach to quantification of the velocity field using Bub-PTV in which the streamwise velocity is inferred by stereoscopic tracking of air bubbles entrained by the flow. This demo code is for bubble tracking of the tests using water injection into a water pool intended to verify the setup of the experiment (e.g. air generating system, stereo cameras) and provide databases for code development and validation.

## System requirements
1. The main functions are written in Matlab (version R2022b)

## Required packages
There is no third party function. Warning of `certain functions not available` can be solved by downloading the specific packages in Matlab.

# How to use!!!
- Run `main.m` for bubble tracking and matching
- Run `postProcessing.m` for post processing the tracked bubbles

## Acknowledgment 
The authors would like to thank the experimental team at Lappeenranta-Lahti University of Technology (Finland) for carrying out SEF experiments and for financial support from NKS (Nordic Nuclear Safety Research) NKS-THEOS project and from Swedish Radiation Safety Authority (SSM). The authors are grateful to Prof. Erik Fransén and Prof. Tony Lindeberg for their expertise in computer vision. 

## Questions
To get help on how to use the data or code, simply open an issue in the GitHub "Issues" section.
