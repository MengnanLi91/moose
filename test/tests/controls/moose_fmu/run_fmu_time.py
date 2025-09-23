from fmpy import extract, read_model_description, instantiate_fmu
from fmpy.simulation import apply_start_values
import numpy as np
import pandas as pd
import logging
import time
import shutil

# Configure root logger
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(name)s: %(message)s'
)
logger = logging.getLogger(__name__)

def moose_fmu_step_by_step(
    moose_filename: str,
    t0: float,
    t1: float,
    dt: float,
    *,
    time_tol: float | None = None,        # None -> auto dt/2
    step_csv: str = "run_fmu_step_by_step.csv"
):
    """
    Manual FMI 2.0 run + comparison with baseline CSV produced by simulate_moose_fmu().
    """
    if time_tol is None:
        time_tol = dt / 2.0

    unzipdir = extract(moose_filename)
    md = read_model_description(unzipdir)
    fmu = instantiate_fmu(unzipdir=unzipdir, model_description=md)

    try:
        vrs = {v.name: v.valueReference for v in md.modelVariables}

        # --- Initialization ---
        fmu.instantiate()
        fmu.setDebugLogging(True, [])
        fmu.setupExperiment(startTime=t0, stopTime=t1)
        fmu.enterInitializationMode()

        apply_start_values(
            fmu=fmu,
            model_description=md,
            start_values={
                "flag":             "TIMESTEP_BEGIN FINAL",
                "moose_executable": "../../../moose_test-opt",
                "moose_inputfile":  "fmu_time.i",
                "server_name":      "web_server",
                "max_retries":      10,
            },
        )

        fmu.exitInitializationMode()
        # ----------------------

        # --- Step loop ---
        rows = []
        t = t0
        step = 0
        while t < t1 - 1e-15:
            print(f"step : {step}")
            fmu.doStep(currentCommunicationPoint=t, communicationStepSize=dt)
            t = min(t + dt, t1)

            moose_time = fmu.getReal([vrs["moose_time"]])[0]

            print(f"fmu_time={t:.3f} → moose_time={moose_time:.6f} ")
            rows.append((t, moose_time))
            step = step + 1

        result = np.array(
            rows,
            dtype=[("time", np.float64), ("moose_time", np.float64)],
        )

        # Save our step-by-step results
        df_step = pd.DataFrame(result)
        df_step.to_csv(step_csv, index=False)

        return result

    finally:
        # Cleanup
        try:
            fmu.terminate()
        except Exception:
            pass
        try:
            fmu.freeInstance()
        except Exception:
            pass
        shutil.rmtree(unzipdir, ignore_errors=True)

def main():

    t0, t1, dt = 0, 2.0, 0.1
    moose_filename = 'MooseTime.fmu'
    result = moose_fmu_step_by_step(moose_filename, t0, t1, dt)

    fmu_time  = result["time"]
    dt        = result["moose_time"]

    for ti, di in zip(fmu_time, dt):
        print(f"fmu_time={ti:.1f} → moose_time={di:.5f} ")


if __name__ == "__main__":
    main()



