import os
import pandas as pd
import altair as alt
import altairThemes
alt.themes.register("publishTheme", altairThemes.publishTheme)
alt.themes.enable("publishTheme")

# report memory and time usage by CCC tools on the human lymph node 
def convert_time(hours, minutes, seconds):
    fractional_minutes = minutes / 60
    fractional_seconds = seconds / 3600
    fractional_hours = hours + fractional_minutes + fractional_seconds

    return fractional_hours 

def plot_usage(usage_dict: dict):
    usage_data = {
        "Method": list(usage_dict.keys()),
        "Time": [metrics["Time"] for metrics in usage_dict.values()],
        "Memory": [metrics["Memory"] for metrics in usage_dict.values()],
    }
    df = pd.DataFrame(usage_data)

    time_order = df.sort_values('Time', ascending=True)['Method'].tolist()
    memory_order = df.sort_values('Memory', ascending=True)['Method'].tolist()

    if "TWCOM" in time_order:
        time_order.remove("TWCOM") # need to add an X in inkscape
        time_order.append("TWCOM")

    if "TWCOM" in memory_order:
        memory_order.remove("TWCOM")
        memory_order.append("TWCOM")


    time_chart = alt.Chart(df).mark_bar().encode(
        x = alt.X('Method', title = 'Method', axis = alt.Axis(labelAngle=-45), sort = time_order), # y: ascending, -y: descending 
        y = alt.Y('Time', title = 'Time (hours)'),
        tooltip = ['Method', 'Time']
    ).properties(
    )

    memory_chart = alt.Chart(df).mark_bar().encode(
        x = alt.X('Method', axis = alt.Axis(labelAngle=-45), title = 'Method', sort = memory_order), # y: ascending, -y: descending 
        y = alt.Y('Memory', title = 'Memory (GB)'),
        tooltip = ['Method', 'Memory']
    ).properties(
    )

    return time_chart, memory_chart 

def main():
    out_dir = "paper_figures/memory_time_usage"
    os.makedirs(out_dir, exist_ok = True)
    usage_dict = {
        "CellChat": {"Time": convert_time(2, 5, 57), "Memory": 4.52}, # memory: GB, job wall clock time
        "Giotto": {"Time": convert_time(6, 29, 8), "Memory": 29.39}, # GB CPU
        "NICHES": {"Time": convert_time(0, 16, 0), "Memory": 17.09}, # GB CPU
        "COMMOT": {"Time": convert_time(1, 16, 55), "Memory": 3.35}, # GB CPU
        "CytoSignal": {"Time": convert_time(0, 4, 0), "Memory": 9.82}, # GB CPU
        "CellNEST": {"Time": convert_time(13, 0, 0), "Memory": 2.11}, # CPU GB
        # "NicheCompass": {"Time": convert_time(0, 5, 0), "Memory": 0.00136}, # CPU GB  
        "TWCOM": {"Time": convert_time(0, 0, 0), "Memory": 0} # GPU GB 

    }
    time_chart, memory_chart = plot_usage(usage_dict)
    time_chart.save(os.path.join(out_dir, "time.html"))
    memory_chart.save(os.path.join(out_dir, "memory.html"))

if __name__ == "__main__":
    main()
