# Multiple Sequence Alignment

The commande `msa` compute multiple sequence alignement of any partition of the pangenome. The command uses [mafft](https://mafft.cbrc.jp/alignment/software/) with default options to perform the alignment. Using multiple cpus with the `--cpu` argument is recommended as multiple alignment can be quite demanding in computational resources.

This command can be used as follow:

```bash
ppanggolin msa -p pangenome.h5
```

By default it will write the strict 'core' (genes that are present in absolutely all genomes) and remove any duplicated genes. Beware however that, if you have many genomes (over 1000), the core will likely be either very small or even empty if you have fragmented genomes.

It will write one MSA for each family. You can then provide the directory where the MSA are written to [IQ-TREE](https://github.com/Cibiv/IQ-TREE) for example, to do phylogenetic analysis.

### Modify the partition with `--partition`

You can change the partition which is written, by using the --partition option.

for example will compute MSA for all the persistent gene families.

```bash
ppanggolin msa -p pangenome.h5 --partition persistent
``` 

Supported partitions are `core`, `persistent`, `shell`, `cloud`, `softcore`, `accessory`. If you need specific filters, you can submit a request in the [issue tracker](https://github.com/labgem/PPanGGOLiN/issues) with your requirements. You can also directly implement the new filter and submit a Pull Request (instructions for contribution can be found [here](../dev/contribute.md)). Most filters should be quite straightforward to add.

### Chose to align dna or protein sequences with `--source`

You can specify whether to use `dna` or `protein` sequences for the MSA by using `--source`. It uses protein sequences by default.

```bash
ppanggolin msa -p pangenome.h5 --source dna
```

### Write a single whole MSA file with `--phylo` 

It is also possible to write a single whole genome MSA file, which many phylogenetic softwares accept as input, by using the `--phylo` option as such:

```bash
ppanggolin msa -p pangenome.h5 --phylo
```

This will contatenate all of the family MSA into a single MSA, with one sequence for each genome.