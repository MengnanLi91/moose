#!/usr/bin/env python3
# This file is part of the MOOSE framework
# https://mooseframework.inl.gov
#
# All rights reserved, see COPYRIGHT for full restrictions
# https://github.com/idaholab/moose/blob/master/COPYRIGHT
#
# Licensed under LGPL 2.1, please see LICENSE for details
# https://www.gnu.org/licenses/lgpl-2.1.html

import contextlib
import shutil

import numpy as np
import pandas as pd
from fmpy import extract, instantiate_fmu, read_model_description
from fmpy.simulation import apply_start_values
from moosefmu import set_real

def moose_fmu_step_by_step(
    moose_filename: str,
    t0: float,
    t1: float,
    dt: float,
    flag: str,
    cmd: str,
    *,
    rtol: float = 1e-6,
    atol: float = 1e-9,
    time_tol: float | None = None,  # None -> auto dt/2
    step_csv: str = "run_fmu_step_by_step.csv",
):
    """Manual FMI 2.0 run + comparison with baseline CSV produced by simulate_moose_fmu()."""
    if time_tol is None:
        time_tol = 1e-12

    moose_model = extract(moose_filename)
    md = read_model_description(moose_model)
    fmu = instantiate_fmu(unzipdir=moose_model, model_description=md)

    try:
        vrs = {v.name: v.valueReference for v in md.modelVariables}

        # --- Initialization ---
        fmu.instantiate()
        fmu.setupExperiment(startTime=t0, stopTime=t1)
        fmu.enterInitializationMode()

        apply_start_values(
            fmu=fmu,
            model_description=md,
            start_values={
                "flag": flag,
                "moose_command": cmd,
                "server_name": "web_server",
                "max_retries": 10,
            },
        )

        fmu.exitInitializationMode()
        # ----------------------

        # --- Step loop ---
        rows = []
        t = t0
        while t < t1 - time_tol:
            step_size = min(dt, t1 - t)
            fmu.doStep(currentCommunicationPoint=t, communicationStepSize=step_size)

            moose_time = fmu.getReal([vrs["moose_time"]])[0]
################# This is an example how to set MOOSE input with MOOSE FMU dome rss model ######
            if t <= 1000:
                set_real(fmu, vrs, "mfr_in", 0.5)
                mfr_in_now = fmu.getReal([vrs["mfr_in"]])[0]
                print("mfr_in (FMU) =", mfr_in_now)
################################################################################################

################# This is an example how to get MOOSE output with MOOSE FMU dome rss model #####
            air_heatrate = fmu.getReal([vrs["air_heatrate"]])[0]
################################################################################################
            print(
                f"fmu_time={t:.3f} -> moose_time={moose_time:.6f} -> air_heatrate={air_heatrate:.6f}"
            )
            rows.append((t, moose_time, air_heatrate))

            t += step_size

        result = np.array(
            rows,
            dtype=[
                ("time", np.float64),
                ("moose_time", np.float64),
                ("air_heatrate", np.float64)
            ],
        )

        # Save our step-by-step results
        df_step = pd.DataFrame(result)
        df_step.to_csv(step_csv, index=False)

        return result

    finally:
        # Cleanup
        with contextlib.suppress(Exception):
            fmu.terminate()
        with contextlib.suppress(Exception):
            fmu.freeInstance()
        shutil.rmtree(moose_model, ignore_errors=True)




if __name__ == "__main__":

# DOME RSS Inputs:
# - mfr_in
# - mfr_out
# - T_in
# - T_out
# - T_air
# - reactor_power
#
# DOME RSS Outputs:
# - air_heatrate; positive value is heat loss from shield
#
# Water inlet/outlet heat rates are calculated as follows:
#   Q = mfr * cp * (T - T_ref)
# where T_ref is taken to be 0. Better would be:
#   Q = mfr * h(T)

    t0, t1, dt = 0, 6000, 2000
    moose_filename = "DomeTest.fmu"
    flag = "MULTIAPP_FIXED_POINT_END"
    cmd = "../../combined-opt -i dome_rss.i"
    result = moose_fmu_step_by_step(
        moose_filename, t0, t1, dt, flag, cmd)

    fmu_time = result["time"]
    moose_time = result["moose_time"]
    air_heatrate = result["air_heatrate"]

    for ti, di, air_heatrate in zip(fmu_time, moose_time, air_heatrate):
        print(
            f"fmu_time={ti:.1f} -> moose_time={di:.5f} -> air_heatrate={air_heatrate:.5f}"
        )
