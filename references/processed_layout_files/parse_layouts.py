print("WARNING: This script is defunct and should no longer be used!")
quit()

import re

# with open("/Users/jansauer/tmp/PROMISE_local/layouts/FDA_1443_edt_V170125_semicolon.txt", "r") as f:
#     dat = f.readlines()

with open("/Users/jansauer/tmp/PROMISE_local/layouts/KiStem_edt_V170125_semicolon.txt", "r") as f:
    dat = f.readlines()

headers = dat[0].split(";")
headers = ";".join(headers[1:4] + ["Row", "Col"] + [headers[-1].strip()])

output = []
for line in dat[1:]:
    linesplit = line.split(";")
    cat_id, compound, well, layout_id = linesplit[1:4] + [linesplit[-1].strip()]
    regex = re.match("([A-Z])([0-9]*)", well)
    row = regex.group(1)
    col = regex.group(2)
    output.append("%s;%s;%s;%s;%s;%s" % (cat_id, compound, well, row, col, layout_id))

layouts = [s.split(";")[-1] for s in output]

for layout_id in set(layouts):
    layout = [s + "\n" for s in output if s.split(";")[-1] == layout_id]
    with open("/Users/jansauer/tmp/PROMISE_local/layouts/L%s.txt" % layout_id, "w-") as f:
        f.write(headers+"\n")
        f.writelines(layout)
