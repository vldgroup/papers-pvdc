import pandas as pd

TRAJ_FILES = [
    "/Users/litongwu/Desktop/PVDC/iter-training/pdn-md/3000-30/log/1_3528_npt_ramp_1200.log",
    "/Users/litongwu/Desktop/PVDC/iter-training/pdn-md/3000-30/log/2_3528_nvt_ramp_3000.log",
    "/Users/litongwu/Desktop/PVDC/iter-training/pdn-md/3000-30/log/3_3528_nvt_hold_3000.log",
    "/Users/litongwu/Desktop/PVDC/iter-training/pdn-md/3000-30/log/4_3528_npt_quench_300.log",
]

STEP_TO_PS = [0.0005, 0.00025, 0.00025, 0.001]
STAGE_NAMES = [
    "npt_ramp_1200",
    "nvt_ramp_3000",
    "nvt_hold_3000",
    "npt_quench_300",
]
TARGET_STAGE_NAME = "target"
OUTPUT_CSV = "md_scheme_with_temp.csv"

full_rows = []
total_global_offset = 0.0

for i, log_file in enumerate(TRAJ_FILES):
    times = []
    temps = []
    file_accumulator = 0.0
    stage_last_ps = 0.0
    last_raw_step = -1

    with open(log_file, "r") as handle:
        header_idx = -1
        step_idx = -1

        for line in handle:
            line = line.strip()
            if "Step" in line and "Temp" in line:
                headers = line.split()
                header_idx = headers.index("Temp")
                step_idx = headers.index("Step")
                continue

            if header_idx != -1:
                parts = line.split()
                if len(parts) > header_idx and parts[0].isdigit():
                    try:
                        raw_step = int(parts[step_idx])
                        temp = float(parts[header_idx])
                        if temp <= 1.0:
                            continue

                        if raw_step < last_raw_step:
                            file_accumulator += stage_last_ps

                        current_ps = raw_step * STEP_TO_PS[i]
                        times.append(current_ps + file_accumulator)
                        temps.append(temp)
                        last_raw_step = raw_step
                        stage_last_ps = current_ps
                    except ValueError:
                        continue

    if times:
        shifted_times = [(t - times[0]) + total_global_offset for t in times]
        print(f"{log_file}: Time {shifted_times[0]:.1f} to {shifted_times[-1]:.1f} ps")

        for t_ps, temp_k in zip(shifted_times, temps):
            full_rows.append(
                {
                    "stage": STAGE_NAMES[i],
                    "time_ps": t_ps,
                    "temp_k": temp_k,
                }
            )

        total_global_offset = shifted_times[-1]

for time_ps, temp_k in [(0, 300), (270, 3000), (370, 3000), (640, 300)]:
    full_rows.append(
        {
            "stage": TARGET_STAGE_NAME,
            "time_ps": time_ps,
            "temp_k": temp_k,
        }
    )

df = pd.DataFrame(full_rows)
df = df.sort_values(["stage", "time_ps"]).reset_index(drop=True)
df.to_csv(OUTPUT_CSV, index=False)
print(f"Saved {len(df)} rows to {OUTPUT_CSV}")
