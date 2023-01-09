# Connected subpartition branch-and-cut
Branch-and-cut algorithm for finding an optimal connected k subpartition of a graph (disjoint vertex subsets inducing k connected subgraphs)


### Dependencies

- The LEMON library from the COIN-OR initiative
- The Gurobi Optimizer for our branch-and-cut framework

### To get LEMON

Download lemon-1.3.1 source code at https://lemon.cs.elte.hu/trac/lemon

Unpack the file (_e.g._ on /opt) and open a terminal on that directory

```
mkdir build
cd build
cmake ../
make
make check
[optional] sudo make install
```
If LEMON is not installed (skipping the last step above), we need to add -I flags to the makefile indicating where to find the corresponding headers.


### To get Gurobi

Download the [Gurobi package](https://www.gurobi.com/downloads/gurobi-software), first creating a login if you don't have one yet. Follow the [installation guide](https://www.gurobi.com/documentation/10.0/quickstart_linux/software_installation_guid.html). In a nutshell, here's what we usually do:

- Unpack the file, _e.g._ on /opt
- Edit the `.profile` and/or `.bashrc` files under your home folder so that the path environment variables correctly point to your installation, for example
```
export GUROBI_HOME="/opt/gurobi1000/linux64"
export PATH="${PATH}:${GUROBI_HOME}/bin"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${GUROBI_HOME}/lib"
```

- Register for a Gurobi license. Assuming you're in academia, you should visit the [license page](https://www.google.com/url?q=https%3A%2F%2Fwww.gurobi.com%2Fdownloads%2Fend-user-license-agreement-academic%2F&sa=D&sntz=1&usg=AOvVaw0YU98KLcE2IKVrvlVaHEjO) from your university network (or through their VPN), and run the grbgetkey tool to validate it. **N.B.** The specific command is indicated at the bottom of the License Detail page you get in the end of this step! Just copy and paste it on the terminal (after restarting your terminal so that the command is recognized).


### To compile and run our software:

```
git clone https://github.com/phillippesamer/connected-subpartition-branch-and-cut.git
cd connected-subpartition-branch-and-cut
```

Edit the Makefile initial definition of variable `GRB_PATH` to reflect the root folder where Gurobi is installed. Finally, compile and run the solver:

```
make
./cks [input file path] [number of subgraphs]
```
