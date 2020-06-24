# mpiBWA

This software allows the alignment of high-throughput sequencing data from fastq files. These files can be pair or single ended, trimmed or not.
`mpiBWA` relies on [BWA MEM](https://github.com/lh3/bwa) and on the [Message Passing Interface (MPI)](https://en.wikipedia.org/wiki/Message_Passing_Interface) standard to perform the parallelisation of the alignment processing over multiple cores and nodes of high performance computing clusters. `mpiBWA` is fully reproducible with the original [BWA MEM](https://github.com/lh3/bwa) version. This means that, if you take the same number of threads in the serial version of BWA and in the MPI version, then you will obtain exactly the same results.


* [Installation](docs/INSTALL.md)
* [Documentation](docs/README.md)
* [Credits](#credits)
* [Acknowledgements](#acknowledgements)
* [Release notes](CHANGELOG.md)
* [Citation](#citation)

## Credits

This program has been developed by Frederic Jarlier from Institut Curie and Nicolas Joly from Institut Pasteur with the help of students from Paris Descartes University and supervised by Philippe Hupé from Institut Curie.

Contacts: [frederic.jarlier@curie.fr](mailto:frederic.jarlier@curie.fr]), [philippe.hupe@curie.fr](mailto:frederic.jarlier@curie.fr])

## Acknowledgements

* [TGCC France Génomique](https://www.france-genomique.org/plateformes-et-equipements/plateforme-tgcc-arpajon/)
* [Intel](https://www.intel.fr/content/www/fr/fr/homepage.html)
* [Université Paris Descartes](https://u-paris.fr/en/498-2)
* [Institut Pasteur](https://www.pasteur.fr)

## Citation

Jarlier F, Joly N, Fedy N et al. QUARTIC: QUick pArallel algoRithms for high-Throughput sequencIng data proCessing. F1000Research 2020, 9:240 ([https://doi.org/10.12688/f1000research.22954.2](https://doi.org/10.12688/f1000research.22954.2))

