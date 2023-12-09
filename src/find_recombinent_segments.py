
##### Find Recombinent Frame Code ######

def find_frame(cigar, pos):
    running_num = 0
    cigar_parsed = []
    for a in cigar:
        if a.isdigit():
            running_num *= 10
            running_num += int(a)
        else:
            cigar_parsed.append([running_num, a])
            running_num = 0
    if len(cigar_parsed) == 1:
        return None
    seq_len = 0
    for a in cigar_parsed:
        if a[1] == 'M' or a[1] == 'I':
            seq_len += a[0]
    left_side = None
    if cigar_parsed[0][1] == 'M' and (float(cigar_parsed[0][0]) / float(seq_len)) > .2:
        left_side = (float(cigar_parsed[0][0]) / float(seq_len))
    right_side = None
    if cigar_parsed[-1][1] == 'M' and (float(cigar_parsed[-1][0]) / float(seq_len)) > .2:
        right_side = (float(cigar_parsed[-1][0]) / float(seq_len))
    if (left_side is None and right_side is None):
        return None
    if right_side is None or (left_side is not None and left_side > right_side):
        return ('end', pos + cigar_parsed[0][0])
    left_mis = 0
    for a in range(len(cigar_parsed)-1):
        if cigar_parsed[a][1] == 'M' or cigar_parsed[a][1] == 'I':
            left_mis += cigar_parsed[a][0]
    return ('beg', pos + left_mis)

##### Reading in SAM File ##############

sam_file = open('data/align/align.sam')

##### Reaching First Line on Align #####

line = sam_file.readline()
while line.startswith('@'):
    line = sam_file.readline()

##### Iterating through Alignments #####

frame_mapping = {}

line_counter = 1
while line:
    alignment_stats = line.strip().split('\t')
    if alignment_stats[5] != '*':
        frame = find_frame(alignment_stats[5], int(alignment_stats[3]))
        if frame is not None:
            if frame[1] not in frame_mapping:
                frame_mapping[frame[1]] = [0,0]
            if frame[0] == 'beg':
                frame_mapping[frame[1]][0] += 1
            else:
                frame_mapping[frame[1]][1] += 1
    line = sam_file.readline()
    line_counter += 1

##### Find Potential Frames ############

sig_frames = {}
for k in frame_mapping:
    if frame_mapping[k][0] > 1500 or frame_mapping[k][1] > 1500:
        sig_frames[k] = frame_mapping[k]

positions = [int(p) for p in list(sig_frames.keys())]
positions.sort()

##### Find Frames ######################

frames = []

open_frame = False
prev_pos = 0
prev_count = 0
for position in positions:
    if open_frame and sig_frames[str(position)][1] > 1500:
        ave = ((float(prev_count) + float(sig_frames[str(position)][1])) / 2.0)
        portion = (ave / float(line_counter))
        frames.append([prev_pos, position, portion])
        open_frame = False
    if sig_frames[str(position)][0] > 1500:
        open_frame = True
        prev_pos = position
        prev_count = sig_frames[str(position)][0]

for frame in frames:
    out = 'There may be a variant call from position ' + str(frame[0])
    out += (' to position ' + str(frame[1]))
    print(out)

