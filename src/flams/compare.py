run3 = {}
with open("log_ptm_counts.txt", "r") as file:
    # --- FIX: Skip the header line ("PTM|Count") ---
    next(file)
    for line in file:
        ptm, count = line.split("|")
        ptm = ptm.replace("__", " ")
        count = int(count.strip())
        if ptm == "Total":
            continue
        run3[ptm] = count

maria = {}
with open("maria_ptm_list.txt", "r") as filem:
    next(filem)
    for line in filem:
        # Note: file2 (ptm_list.txt) correctly uses the tab delimiter for its data.
        ptm, count = line.split("\t") 
        ptm = ptm.replace('"', "")
        count = int(count.strip())
        maria[ptm] = count

run3_set = set(run3.keys())
maria_set = set(maria.keys())

intersection = run3_set & maria_set
unique_run3 = run3_set - maria_set
unique_maria = maria_set - run3_set

print(f"intersection: {intersection}")
print(f"intersection count: {len(intersection)}")

print(f"unique run3: {unique_run3}")
print(f"unique run3 count: {len(unique_run3)}")

print(f"unique maria: {unique_maria}")
print(f"unique maria count: {len(unique_maria)}")

union = maria_set | run3_set

headers = ["PTM", "Maria", "Ieva"]

data = []
data.append(headers)

for entry in union:
    row = [entry, maria.get(entry, "N/A"), run3.get(entry, "N/A")]
    data.append(row)

import csv

with open("comaprison.tsv", "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")  # use tab instead of comma
    writer.writerows(data)
