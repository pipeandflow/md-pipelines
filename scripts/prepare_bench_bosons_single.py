import itertools
from utils import random_boson_positions

def flatten(list_of_lists):
    return list(itertools.chain.from_iterable(list_of_lists))

def stringify_elementwise(lst):
    return [str(x) for x in lst]


def bench_bosons_single(nbosons, boson_positions, socket_id, workdir):
    assert len(boson_positions) == nbosons

    boson_positions_config = str(flatten(boson_positions))

    boson_labels = str(["E"] * nbosons)
    boson_masses = str(["1.0"] * nbosons)
    bosons_list = str(list(range(nbosons)))

    config = ipi_config(boson_positions_config, boson_masses, boson_labels, bosons_list, socket_id, workdir)
    return config

def ipi_config(boson_positions, boson_masses, boson_labels, bosons_list, socket_id, workdir):
    return  f"""<!--REGTEST
COMMAND(4)    i-pi-driver -u -h REGTEST_SOCKET -m harm3d -o 1.21647924E-8
ENDREGTEST-->
<simulation threading='False' verbosity='low'>

<ffsocket mode='unix' name='driver'>
        <address> {socket_id}  </address>
</ffsocket>

<total_steps> 100 </total_steps>

<output prefix="{workdir}/data">
  <trajectory stride="100" filename="pos" cell_units="angstrom">positions{{angstrom}}</trajectory>
  <!--<trajectory stride="1" filename="xc" format="xyz">x_centroid{{angstrom}}</trajectory>-->
  <properties stride="100"> [ step, time{{femtosecond}}, conserved, temperature{{kelvin}}, kinetic_cv,
        potential, virial_fq ] </properties>
</output>

<prng>
  <seed> 18885 </seed>
</prng>

<system>

  <forces>
      <force forcefield="driver"></force>
  </forces>

  <initialize nbeads="32">
    <positions mode="manual" bead="0"> {boson_positions} </positions>
    <masses mode="manual"> {boson_masses} </masses>
    <labels mode="manual"> {boson_labels} </labels>
    <cell>
     [   2500, 0, 0, 0, 2500, 0, 0, 0, 2500 ]
    </cell>
    <velocities mode='thermal' units='kelvin'> 17.4 </velocities>
  </initialize>

  <normal_modes propagator='bab'>
          <nmts> 10 </nmts>
          <bosons> {bosons_list} </bosons>
  </normal_modes>

  <ensemble>
     <temperature units="kelvin"> 17.4 </temperature>
  </ensemble>

  <motion mode="dynamics">
    <fixcom> False </fixcom>
    <dynamics mode="nvt">
     <timestep units="femtosecond"> 1 </timestep>
      <thermostat mode='pile_l'>
            <tau units='femtosecond'>100</tau>
      </thermostat>

    </dynamics>
  </motion>

</system>

</simulation>
"""

if __name__ == '__main__':
    num_bosons = int(snakemake.params.num_bosons)
    boson_positions = random_boson_positions(num_bosons)
    config = bench_bosons_single(nbosons=num_bosons, boson_positions=boson_positions,
    socket_id=int(snakemake.params.socket_id), workdir=snakemake.params.workdir)
    with open(snakemake.output[0], 'w') as f:
        f.write(config)
