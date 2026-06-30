<!-- 
Describe the dataset. See previous weeks for the general format of the
description. The description is the part of the readme.md file above "The Data";
everything else will be filled in from the other md files in this directory +
automatic scripts. We usually include brief introduction along the lines of
"This week we're exploring DATASET" or "The dataset this week comes from 
SOURCE".
-->

This data set gives morphometric data for 93 penguins from 18 species within 6 genera. It was inspired by the
now-classic "Palmer penguins data set". I attended a workshop where students were analyzing the Palmer-penguins
data set with a hierarchical model with individuals grouped by species. However, because there are only three species
represented in the Palmer data set (Adélie, Chinstrap, and Gentoo), this data set is not ideal for that purpose.
I found (with the assistance of a chatbot) the [AVONET dataset](https://opentraits.org/datasets/avonet.html) (Tobias
*et al.* 2022):

> The AVONET database contains comprehensive functional trait data for all birds, including six ecological variables, eleven continuous morphological traits, and information on range size and location. Raw morphological measurements are available from 90020 individuals of 11009 extant bird species sampled from 181 countries.

I selected the penguin data from the database. The data set has 10 different morphometric measurements of penguin beaks, wings, tails, etc. (although up to 12% of some measurements are missing). There is also crude phylogenetic data (in a Newick format data file, readable with the`ape` package (although if you're new to phylogenetic data, you might want to skip this unless you want a challenge/to learn something new); getting spatial/geographic species range data would take more work.

How do trait values covary within/across species and genera? Is there a good way to do ordination/visualization that handles the missingness of some of the traits nicely? Are there interesting ways to visualize these data in >2 dimensions?
