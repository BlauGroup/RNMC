# GMC - <span style="color: #0066CC"> Homogeneous Catalysis </span>

We apply our `GMC` module to simulating [cobalt-catalyzed hydroformylation](./https://pubs.rsc.org/en/content/articlehtml/2017/sc/c7sc03628k). This application is choosen as an example because it has widely been studied and is a relatively simple organometallic reaction.

The following diagram outlines the reaction network of cobalt-catalyzed hydroformylation. The circles indicate different species and the solid lines indicate reactions. The numerical values inside the circles correspond to the species IDs used in the Python code to generate the .sqlite files to run `GMC` (see <a href="{{ site.github.repository_url }}"> examples directory</a>).

<figure>
    <img src="catalysis.pdf"
         alt="homogeneous catalysis">
    <figcaption> Final reaction network for the hydroformylation reaction.  </figcaption>
</figure>