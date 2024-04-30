import argparse
import pandas as pd

def parse_duration(duration_str):
    """ Convert duration from 'XXh XXm XXs' to total seconds.
    
    Args:
        duration_str (str): A string representing the duration, which can include hours (h),
                            minutes (m), and seconds (s). Parts of the duration may be omitted.

    Returns:
        int: The total duration in seconds.
    """
    total_seconds = 0
    parts = duration_str.split()

    # Parse hours, minutes, and seconds
    for part in parts:
        if 'h' in part:
            hours = int(part.replace('h', ''))
            total_seconds += hours * 3600
        elif 'm' in part:
            minutes = int(part.replace('m', ''))
            total_seconds += minutes * 60
        elif 's' in part:
            seconds = int(part.replace('s', ''))
            total_seconds += seconds

    return total_seconds

def parse_memory(mem_str):
    """ Convert memory usage from 'X.Y GB' to float gigabytes. """
    return float(mem_str.split()[0])

def main(trace_file_path):
    # Read the trace file
    df = pd.read_csv(trace_file_path, delimiter='\t')
    
    # Parse duration, peak_vmem, and %cpu columns
    df['duration_seconds'] = df['duration'].apply(parse_duration)
    df['peak_vmem_gb'] = df['peak_vmem'].apply(parse_memory)
    df['%cpu'] = df['%cpu'].str.replace('%', '').astype(float)  # Convert %cpu to float
    
    # Calculate mean of duration in seconds, peak virtual memory in GB, and %cpu
    mean_duration = df['duration_seconds'].mean()
    mean_peak_vmem = df['peak_vmem_gb'].mean()
    mean_cpu = df['%cpu'].mean()
    
    # Print the results in a table format
    print(f"{'Column':<20}{'Mean Value'}")
    print(f"{'-'*31}")
    print(f"{'Duration (seconds)':<20}{mean_duration:.2f}")
    print(f"{'Peak VMEM (GB)':<20}{mean_peak_vmem:.2f}")
    print(f"{'CPU Usage (%)':<20}{mean_cpu:.2f}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process a Nextflow trace file and calculate mean values for duration, peak virtual memory, and CPU usage.')
    parser.add_argument('trace_file_path', type=str, help='Path to the trace file')
    args = parser.parse_args()
    main(args.trace_file_path)
