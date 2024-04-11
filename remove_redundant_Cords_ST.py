#!/usr/bin/env python3

def merge_overlapping_ranges(coords_file, output_file):
    seq_coords = {}
    with open(coords_file) as f:
        for line in f:
            line = line.strip()
            seq_name, coords = line.split(':')
            start, end = map(int, coords.split('-'))
            if seq_name in seq_coords:
                seq_coords[seq_name].append((start, end))
            else:
                seq_coords[seq_name] = [(start, end)]
    with open(output_file, 'w') as f:
        for seq_name, coords in seq_coords.items():
            coords.sort()
            merged_coords = []
            prev_start, prev_end = coords[0]
            for start, end in coords[1:]:
                if start <= prev_end:
                    prev_end = max(prev_end, end)
                else:
                    merged_coords.append((prev_start, prev_end))
                    prev_start, prev_end = start, end
            merged_coords.append((prev_start, prev_end))
            for start, end in merged_coords:
                f.write(f"{seq_name}:{start}-{end}\n")


if __name__ == "__main__":
    import sys
    coords_file = sys.argv[1]
    output_file = sys.argv[2]
    merge_overlapping_ranges(coords_file, output_file)
