import os

def countConfusionMatrix(out_dir, ext):
    for file in os.listdir(out_dir):
        cnt1 = 0
        cnt2 = 0
        if file.endswith(ext):
            filename = file[:-4]
            with open(os.path.join(out_dir, file)) as f:
                for line in f:
                    # chr22 indicates reference-originated reads
                    if line.startswith("chr22"):
                        # fp for contaminated file and tn for uncomtinaminated file
                        cnt1 += 1
                    else:
                        # tp for contaminated file and fn for uncomtinaminated file
                        cnt2 += 1
            if ext == ".contaminated.txt":
                with open(os.path.join(out_dir, filename + ".count.txt"), "w") as f:
                    f.write(f"fp: {cnt1}\ntp: {cnt2}\n")
            else: 
                with open(os.path.join(out_dir, filename + ".count.txt"), "w") as f:
                    f.write(f"tn: {cnt1}\nfn: {cnt2}\n")
    return


output_dir = "data/comparisons/Vecuum"
countConfusionMatrix(output_dir, ".contaminated.txt")
countConfusionMatrix(output_dir, ".noncontaminated.txt")