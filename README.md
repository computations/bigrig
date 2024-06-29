# Introduction

`bigrig` is a program to simulate historical biogeographical regions under the
DEC(+J) model. 

## A note on terminology

The terminology used around this topic can be confusing. Here I am going to
define some terms that I use, and relate them to the original Ree _et. al_
paper [1].


- **Region**: A segment of some area.
- **Range**: A specific assignment regions where a taxa can be found. The set of
  regions that a taxa occupies. Practically, for this program, a range is a
  binary string (e.g. 01010).
- **Singleton**: A range consisting of only one region (e.g. 01000).
- **Split**: A split is a cladogenesis event. This is to say, a fork in the
  binary tree. When a split occurs, a cladogenesis type is used to generate the
  daughter ranges.
- **Allopatry**: A type cladogenesis where the daughter ranges share no regions
  in common. Ree _et. al_ calls this "Scenario 2".
- **Sympatry**: A type of cladogenesis where the daughter ranges share _some_
  regions. Ree _et. al_ calls this "Scenario 3".
- **Copy**: A type of cladogenesis where the daughter ranges share _all_
  regions. Ree _et. al_ calls this "Scenario 1".
- **Jump**: A type of cladogenesis where the daughter ranges share no regions,
  and _one_ of the daughter range shares no regions with the parent range. Ree
  _et. al_ do not include this type in the original paper, and was introduced by
  Matzke [2].

For split types, it can be useful to write down their properties as a binary
operations. 

- **Allopatry**: `left ^ right == parent`.
- **Sympatry**: `(left | right == parent) && (left & right != parent)`.
- **Copy**: `left & right == parent`.
- **Jump**: `(left & parent == 0 || right & parent == 0) && (left == parent || right == parent)`.

Note that an additional requirement is that _at least_ one of the daughter
ranges must be a singleton. This is to say, speciation will only occur in one
region. Furthermore, the "copy" type split can only occur when the parent range
is a singleton.
 
[1]: https://doi.org/10.1111/j.0014-3820.2005.tb00940.x
[2]: https://doi.org/10.1093/sysbio/syu056

# Options

The required parameters are: 

- `--tree`: Path to the tree file used for the simulation
- `--root-range`: Range for the species at the root of the tree. Alternatively,
  the starting range. Only required if `range-count` is not specified.
- `--range-count`: Number of regions to simulate with. If `root-range` is not
  specified, then a random root range is generated.
- `-d/--dispersion`: Dispersion rate for the simulation.
- `-e/--extinction`: Extinction rate for the simulation.
- `-v/--allopatry`: Allopatry/vicariance rate for the simulation.
- `-s/--sympatry`:  Sympatry rate for the simulation.
- `-y/--copy`: Copy rate for the simulation.
- `-j/--jump`:  rate for the simulation.

Additionally, there are a number of optional parameters:

- `--config`: (Optional) Pass a YAML file containing the configuration for the
  program. The details for this file are detailed later.
- `--prefix`: (Optional) Prefix for the results file.

# Config file

The config file is an alternative way of specifying the program options. Options
written into a YAML file, which is then passed to bigrig using the `--config`
switch. The name of options differ slightly between the command line and the
config file. Here is the schema:

```.yaml
periods:
  - start: <FLOAT>
    rates:
      dispersion: <FLOAT>
      extinction: <FLOAT>
    cladogenesis:
      allopatry: <FLOAT>
      sympatry: <FLOAT>
      copy: <FLOAT>
      jump: <FLOAT>
  <...>
root-range: <ROOT-RANGE>
tree: <FILE>
redo: <BOOL>
debug-log: <BOOL>
output-format: [YAML|JSON]
prefix: <PATH>
mode: [FAST|SIM]
seed: <INT>
```

If both the a command line option and a config option are set, for example in
the command

```
bigrig --config config.yaml --redo
```

Then the value from the command line is used instead. The idea here is to have a
"base" config with all the normal values, and the other values can be played
with rapidly by changing them on the command line. When this is done, there is a
warning that is emitted, like this:

```
[WARN] The 'redo' option is specified in both the config file and the command line. Using the value from the command line
```

# Result files

A given simulation will always produce the following result files:

- `{prefix}.phy`: An alignment containing the tip ranges. This is to say, the
  ranges of the "extant" taxa.
- `{prefix}.all.phy`: An alignment containing _all_ ranges, including the inner
  nodes. If the inner nodes have no label, they are labeled by an internal id
  number, starting with 0. Which node is has which label is noted in the
  `{prefix}.annotated.nwk` file.
- `{prefix}.annotated.nwk`: This file contains:
  - Inner node labels
  - NHX encoded splits, with the keys:
    - `init-range`: Range of the species at the "top" of the node. The range to
      split.
    - `left-range`: The "left" result of the split. 
    - `right-range`: The "right" result of the split.
    - `split-type`: One of `allopatric`, `sympatric`, `singleton` or `jump`.

Additionally, if specified as a runtime parameter `bigrig` will produce a
`{prefix}.json` or `{prefix}.yaml`, which is will contain all the information in
the other files and information about dispersion and extinction events.

## An example run

Suppose we have the tree file `test.nwk`

```
((d:0.5058,(b:0.3126,c:0.5691):0.4320):0.5163,(e:0.8394,a:0.9062):0.1542);
```

and the config file `config.yaml`

```
rates:
  dispersion: 1.0
  extinction: 1.0
cladogenesis:
  allopatry: 1.0
  copy: 1.0
  jump: 1.0
  sympatry: 1.0
root-range: 11010
tree: test.nwk
```

then command `./bigrig --config config.yaml` might produce


```
[  0.00s] Running simulation with the following options:
[  0.00s]    Tree file: /home/user/wrk/test.nwk
[  0.00s]    Prefix: /home/user/wrk/test.nwk
[  0.00s]    Root range: 11010
[  0.00s]    Region count: 5
[  0.00s]    Rate parameters:
[  0.00s]        Dispersion(d): 1.00, Extinction(e): [1 00](1.00)
[  0.00s]    Cladogenesis parameters:
[  0.00s]        Allopatry(v): 1.00, Sympatry(s): 1.00, Copy(y): 1.00, Jump(j): 1.00
[  0.00s] Parsing tree
[  0.00s] Sampling from tree
[  0.00s] Writing results to files
[  0.00s] Done!
```

In addition to the log on `stdout`, there will be the 3 results files, as
described above. The first, `test.nwk.phy` might look like

```
5 5
a 10010
e 10101
c 00001
b 00100
d 11110
```

Notice that there are no extinct taxa. This is intentional, as all tips are
assumed to be modern and extant.

The other result file, `test.nwk.all.phy` might look like

```
9 5
0 11010
3 11010
a 10010
e 10101
1 10010
2 00010
c 00001
b 00100
d 11110
```

Here, the numbers refer to the internal nodes. To see which internal node gets
which label, we need to look at `test.nwk.annotated.nwk`.

```
((d[&&NHX:dist=11110],(b[&&NHX:dist=00100],c[&&NHX:dist=00001])2[&&NHX:parent-range=00010:left-range=00100:right-range=00010:split-type=jump])1[&&NHX:parent-range=10010:left-range=10010:right-range=01000:split-type=jump],(e[&&NHX:dist=10101],a[&&NHX:dist=10010])3[&&NHX:parent-range=11010:left-range=01000:right-range=10010:split-type=allopatric])0[&&NHX:parent-range=11010:left-range=00010:right-range=11010:split-type=sympatric]
```

This file contains the tree that was simulated, as well as the ranges for each
node in the tree, including tips. In addition, the full details of each split is
included.
