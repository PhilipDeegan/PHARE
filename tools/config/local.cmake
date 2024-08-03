cmake_minimum_required (VERSION 3.20.1)
project(configured_phare)
message("")
message("!!PHARE CONFIGURATED!!")
message("")

set (PHARE_MPIRUN_POSTFIX "${PHARE_MPIRUN_POSTFIX} --bind-to none")
                