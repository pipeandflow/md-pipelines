import subprocess
import time
import csv
import os
import shutil
from utils import set_logger


def time_ipi(input_filename, socket_id):
    # Run i-pi

    # duplicated from ipi/ipi_tests/test_tools.py: Runner.run()

    call_ipi="i-pi " + input_filename
    clientcall = f"i-pi-driver -u -m harm3d -o 1.21647924E-8 -h {socket_id}"

    logger.debug("Start i-pi server")

    ipi = subprocess.Popen(
        call_ipi,
        # cwd=(self.tmp_dir),
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    # if len(clients) > 0:
    #     f_connected = False
    #     for client in clients:
    #         for i in range(50):
    #             if os.path.exists("/tmp/ipi_" + client[2]):
    #                 f_connected = True
    #                 break
    #             else:
    #                 time.sleep(0.5)
    #         if not f_connected:
    #             print("Could not find the i-PI UNIX socket.")
    #             return "Could not find the i-PI UNIX socket"

    # Run drivers by defining cmd2 which will be called, eventually
    driver = list()

    # for client in clients:
    #     if client[1] == "unix":
    #         clientcall = call_driver + " -m {} {} {} -u ".format(
    #             client[0], address_key, client[2]
    #         )
    #     elif client[1] == "inet":
    #         clientcall = call_driver + " -m {} {} {} -p {}".format(
    #             client[0], client[2], address_key, client[3]
    #         )

    #     else:
    #         raise ValueError("Driver mode has to be either unix or inet")

    cmd = clientcall

    # Add extra flags if necessary
    # if any("-" in str(s) for s in client):
    #     flag_indeces = [
    #         i for i, elem in enumerate(client) if "-" in str(elem)
    #     ]
    #     for i, ll in enumerate(flag_indeces):
    #         if i < len(flag_indeces) - 1:
    #             cmd += " {} {}".format(
    #                 client[ll],
    #                 ",".join(client[ll + 1 : flag_indeces[i + 1]][:]),
    #             )
    #         else:
    #             cmd += " {} {}".format(
    #                 client[ll], ",".join(client[ll + 1 :][:])
    #             )
    # print("cmd:", cmd)

    time.sleep(5)
    logger.debug("Start driver execution")


    start_time = time.time()

    driver.append(
        subprocess.Popen(cmd,
                         # cwd=(cwd),
                         shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    )

    logger.debug("Start wait for ipi process")

    ipi_error = ipi.communicate(timeout=snakemake.config['SINGLE_BENCH_TIMEOUT_SECONDS'])[1].decode("ascii")
    if ipi_error != "":
        logger.debug(ipi_error)
    # assert "" == ipi_error, "IPI ERROR OCCURRED: {}".format(ipi_error)

    end_time = time.time()

    logger.debug("i-pi run ended; time %s" % str(end_time - start_time))

    return end_time - start_time


def bench_ipi():
    # code duplicated from ipi/ipi_tests/test_tools.py: Runner.run()
    measured_time = time_ipi(snakemake.input[0], int(snakemake.params.socket_id))
    with open(snakemake.output[0], 'w') as csv_log:
        w = csv.DictWriter(csv_log, ["nbosons", "time"])
        w.writeheader()
        w.writerow({"nbosons": snakemake.wildcards.num_bosons, "time": measured_time})
        csv_log.flush()

#####################################################################3333

if __name__ == '__main__':
    logger = set_logger(snakemake.log[0])
    workdir = snakemake.params.workdir
    if os.path.isdir(workdir):
        shutil.rmtree(workdir)
    os.mkdir(workdir)
    bench_ipi()
