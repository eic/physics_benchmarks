import sys

def filter_events(input_file,output_file, particle_id=80000, status=-3):
    with open(input_file, 'r') as f:
        lines = f.readlines()

    events = []
    current_event = []

    for line in lines:
        if line.startswith('E '):
            if current_event:
                events.append(current_event)
                current_event = []
        current_event.append(line)

    # Append the last event
    if current_event:
        events.append(current_event)

    filtered_events = []
    filtered_events.append(lines[0:4])
    for event in events:
        particle_found = False
        for line_index, line in enumerate(event):
            if line.startswith('P '):
                parts = line.split()
                if int(parts[3]) == particle_id and int(parts[2]) == status:
                    particle_found = True
                    parts[2] = '-2'  # Change the third index to -2
                    parts[-1] = '200'  # Change the last index to 200
                    event[line_index] = ' '.join(parts) + '\n'  # Join the parts back into a line
                    break
        if particle_found:
            filtered_events.append(event)

    with open(output_file,'w') as f:
        for event in filtered_events:
            for line in event:
                f.write(line)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python filter_hepmc.py input.hepmc")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    filter_events(input_file,output_file)
 
