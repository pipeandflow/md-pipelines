from snakemake.utils import min_version
min_version("7.20.0")

configfile: "config.yaml"

print(config)
def get_socket_id(wildcards, *args, **argw):
    return 50000 + int(wildcards.num_bosons)
#----------------------------------------------------------------

rule all:
    input:
        "boson_scalabilty.png"


rule plot_scalabilty:
    input:
        expand("results/{num_bosons}_bosons_timing.dat", num_bosons=config['boson_numbers'])
    output:
        rules.all.input
    script:
        "scripts/plot_scalability.py"
    

rule time_boson_ipi_run:
    input:
        "tmp/ipi_input_{num_bosons}_bosons.xml"
    output:
        "results/{num_bosons}_bosons_timing.dat"
    params:
        socket_id=get_socket_id,
        workdir="tmp/ipi-run-{num_bosons}/"
    threads: 2
    script:
        "scripts/ipi_run_and_time.py"
        

rule prepare_input_for_ipi:
    output:
        rules.time_boson_ipi_run.input
    params:
        num_bosons="{num_bosons}",
        socket_id=get_socket_id,
        workdir=rules.time_boson_ipi_run.params.workdir
    script:
        "scripts/prepare_bench_bosons_single.py"
        